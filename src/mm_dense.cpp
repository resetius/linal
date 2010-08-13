/*
 * Copyright 1993-2010 NVIDIA Corporation.  All rights reserved.
 *
 * NVIDIA Corporation and its licensors retain all intellectual property and 
 * proprietary rights in and to this software and related documentation. 
 * Any use, reproduction, disclosure, or distribution of this software 
 * and related documentation without an express license agreement from
 * NVIDIA Corporation is strictly prohibited.
 *
 * Please refer to the applicable NVIDIA end user license agreement (EULA) 
 * associated with this source code for terms and conditions that govern 
 * your use of this NVIDIA software.
 * 
 */

#include <omp.h>

////////////////////////////////////////////////////////////////////////////////
// export C interface
extern "C"
void computeGold( double*, const double*, const double*, unsigned int, unsigned int, unsigned int);

////////////////////////////////////////////////////////////////////////////////
//! Compute reference data set
//! C = A * B
//! @param C          reference data, computed but preallocated
//! @param A          matrix A as provided to device
//! @param B          matrix B as provided to device
//! @param hA         height of matrix A
//! @param wB         width of matrix B
////////////////////////////////////////////////////////////////////////////////
void
computeGold(double* C, const double* A, const double* B, unsigned int hA, unsigned int wA, unsigned int wB)
{
	for (unsigned int i = 0; i < hA; ++i) {
		for (unsigned int j = 0; j < wB; ++j) {
			double sum = 0;
			for (unsigned int k = 0; k < wA; ++k) {
				double a = A[i * wA + k];
				double b = B[k * wB + j];
				sum += a * b;
			}
			C[i * wB + j] = sum;
		}
	}
}

#define block_dim 128
static double As[block_dim][block_dim];
static double Bs[block_dim][block_dim];
static double Cs[block_dim][block_dim];

extern "C"
void
computeGold2(double * C, const double * A, const double * B, int n)
{
	int blocks = (n + block_dim - 1) / block_dim;

#pragma omp parallel
{

#pragma omp for
	for (int i = 0; i < n * n; ++i)
	{
		C[i] = 0;
	}

	for (int l = 0; l < blocks; ++l )
	{
		int fl = n * l;
		fl /= blocks;
		int ll = n * (l + 1);
		ll = ll / blocks - 1;
		int nll = ll - fl + 1;

		for (int m = 0; m < blocks; ++m)
		{
			int fm = n * m;
			fm /= blocks;
			int lm = n * (m + 1);
			lm = lm / blocks - 1;
			int nlm = lm - fm + 1;

			// C[l, m] = \sum A[l, k] * B[k, m]

			for (int k = 0; k < blocks; ++k)
			{
				int fk = n * k;
				fk /= blocks;
				int lk = n * (k + 1);
				lk = lk / blocks - 1;
				int nlk = lk - fk + 1;

#pragma omp for schedule(static)
				for (int i = 0; i < nll; ++i)
				{
					for (int j = 0; j < nlk; ++j)
					{
						As[i][j] = A[(i + fl) * n + (j + fk)];
					}
				}

#pragma omp for schedule(static)
				for (int i = 0; i < nlm; ++i)
				{
					for (int j = 0; j < nlk; ++j)
					{
						Bs[i][j] = B[(i + fk) * n + (j + fm)];
					}
				}

#pragma omp barrier

#pragma omp for schedule(static)
				for (int i = 0; i < nll; ++i)
				{
					for (int j = 0; j < nlm; ++j)
					{
						double s = 0.0;
						for (int k1 = 0; k1 < nlm; ++k1)
						{
							s += As[i][k1] * Bs[k1][j];
						}
	//					Cs[i][j] = s;
						C[(i + fl) * n + j + fm] += s;
					}
				}
#pragma omp barrier
				/*
#pragma omp for schedule(static)
				for (int i = 0; i < nll; ++i)
				{
					for (int j = 0; j < nlm; ++j)
					{
						C[(i + fl) * n + j + fm] += Cs[i][j];
					}
				}
				*/

#pragma omp barrier
			}
		}
	}
}

}

#include <time.h>

#if defined(_MSC_VER) || defined(_MSC_EXTENSIONS)
#define DELTA_EPOCH_IN_MICROSECS  11644473600000000Ui64
#else
#define DELTA_EPOCH_IN_MICROSECS  11644473600000000ULL
#endif

#ifdef _MSC_VER
#include <windows.h>

	struct timezone
	{
		int  tz_minuteswest; /* minutes W of Greenwich */
		int  tz_dsttime;     /* type of dst correction */
	};

	int gettimeofday (struct timeval *tv, struct timezone *tz)
	{
		FILETIME ft;
		unsigned __int64 tmpres = 0;
		static int tzflag;

		if (NULL != tv)
		{
			GetSystemTimeAsFileTime (&ft);

			tmpres |= ft.dwHighDateTime;
			tmpres <<= 32;
			tmpres |= ft.dwLowDateTime;

			/*converting file time to unix epoch*/
			tmpres /= 10;  /*convert into microseconds*/
			tmpres -= DELTA_EPOCH_IN_MICROSECS;
			tv->tv_sec = (long) (tmpres / 1000000UL);
			tv->tv_usec = (long) (tmpres % 1000000UL);
		}

		if (NULL != tz)
		{
			if (!tzflag)
			{
				_tzset();
				tzflag++;
			}
			tz->tz_minuteswest = _timezone / 60;
			tz->tz_dsttime = _daylight;
		}

		return 0;
	}
#else
#include <sys/time.h>
#endif

extern "C"
	double get_full_time()
	{
		struct timeval tv;
		gettimeofday (&tv, 0);
		return (double) (tv.tv_sec * 100.0 + tv.tv_usec / 10000.0);
	}


#ifdef TEST
#include <vector>
#include <stdio.h>

using namespace std;
#define N 4096L

int main()
{
	vector < double > A(N * N);
	vector < double > B(N * N);
	vector < double > C(N * N);
	double t  = get_full_time();
	omp_set_num_threads(4);
	computeGold2(&C[0], &A[0], &B[0], N);
	printf("%lf\n", (get_full_time() - t) / 100.0);
}
#endif

