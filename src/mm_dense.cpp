/* -*- charset: utf-8 -*- */
/*$Id$*/

/* Copyright (c) 2009 Alexey Ozeritsky (Алексей Озерицкий)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. The name of the author may not be used to endorse or promote products
 *    derived from this software without specific prior written permission
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
 * NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <string.h>
#include <math.h>
#include <stdio.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "timer.h"

namespace linal
{

void
mat_mult_mat_stupid(double * C, const double * A, const double * B, int n)
{
#pragma omp parallel for
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			double s = 0;
			for (int k = 0; k < n; ++k) {
				s += A[i * n + k] * B[k * n + j];
			}
			C[i * n + j] = s;
		}
	}
}

template < typename T >
void mat_mult_mat_(T * C, const T * A, const T * B, int n)
{
#define block_dim 64
	T As[block_dim][block_dim];
	T Bs[block_dim][block_dim];
	T Cs[block_dim][block_dim];

	int blocks = (n + block_dim - 1) / block_dim;
#undef block_dim

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

#pragma omp for schedule(static)
			for (int i = 0; i < nll; ++i)
			{
				memcpy(&Cs[i][0], &C[(i + fl) * n + fm], nlm * sizeof(T));
			}

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
					memcpy(&As[i][0], &A[(i + fl) * n + fk], nlk * sizeof(T));
				}

#pragma omp for schedule(static)
				for (int i = 0; i < nlm; ++i)
				{
					memcpy(&Bs[i][0], &B[(i + fk) * n + fm], nlk * sizeof(T));
				}

#pragma omp barrier

#pragma omp for schedule(static)
				for (int i = 0; i < nll; ++i)
				{
					for (int j = 0; j < nlm; ++j)
					{
						T s = 0.0;
						for (int k1 = 0; k1 < nlm; ++k1)
						{
							s += As[i][k1] * Bs[k1][j];
						}
						Cs[i][j] += s;
					}
				}
#pragma omp barrier
			}

#pragma omp for schedule(static)
			for (int i = 0; i < nll; ++i)
			{
				memcpy(&C[(i + fl) * n + fm], &Cs[i][0], nlm * sizeof(T));
			}
#pragma omp barrier
		}
	}
}

}

void mat_mult_mat(double * C, const double * A, const double * B, int n)
{
	mat_mult_mat_(C, A, B, n);
}

void mat_mult_mat(float * C, const float * A, const float * B, int n)
{
	mat_mult_mat_(C, A, B, n);
}

}

static int get_num_threads ()
{
#ifdef _OPENMP
	return omp_get_num_threads();
#else
	return 1;
#endif
}

static int get_my_id ()
{
#ifdef _OPENMP
	return omp_get_thread_num();
#else
	return 1;
#endif
}

void
mat_mult_mat_cannon(double * C, const double * A, const double * B, int n)
{
#define block_dim 64
#pragma omp parallel
{
	int blocks = (n + block_dim - 1) / block_dim;
	int id = get_my_id();
	int threads = get_num_threads();

	int first_row = n * id;
	first_row /= threads;
	int last_row = n * (id + 1);
	last_row = last_row / threads - 1;

#pragma omp for
	for (int i = 0; i < n * n; ++i) {
		C[i] = 0.0;
	}

//	for (int l = 0; l < threads; ++l) {
//		int k = (id + l) % threads;
//		int fk = n * k;
//		fk /= threads;
//		int lk = n * (k + 1);
//		lk = lk / threads - 1;
//		int nlk = lk - fk + 1;

	for (int l = 0; l < blocks; ++l) {
		int k = (id + l) % blocks;
		int fk = n * k;
		fk /= blocks;
		int lk = n * (k + 1);
		lk = lk / blocks - 1;
		int nlk = lk - fk + 1;

		for (int c_col = 0; c_col < n; c_col += block_dim)
		{
			int c_ncol = (c_col + block_dim <= n ? c_col + block_dim : n) - 1;

			for (int c_row = first_row; c_row <= last_row; c_row += block_dim)
			{
				int c_nrow = (c_row + block_dim - 1 <= last_row 
					? c_row + block_dim - 1 : last_row);

				for (int a_col = fk; a_col <= lk; a_col += block_dim)
				{
					int a_ncol = (a_col + block_dim - 1 <= lk 
						? a_col + block_dim - 1 : lk);

					int m = c_col;
					const double * pb = B + m;
					double * pc = C + m;

					for (; m <= c_ncol; m += 2, pb += 2, pc += 2)
					{
						for (int i = c_row; i <= c_nrow; i += 2)
						{
							double s00, s01, s10, s11;
							s00 = s01 = s10 = s11 = 0.0;

							int j = a_col;
							const double * pa = A + i * n + j;
							for (; j <= a_ncol; ++j, ++pa)
							{
								s00 += pa[0] * pb[j * n];
								s01 += pa[0] * pb[j * n + 1];
								s10 += pa[n] * pb[j * n];
								s11 += pa[n] * pb[j * n + 1];
							}

							pc[i * n]           += s00;
							pc[i * n + 1]       += s01;
							pc[(i + 1) * n]     += s10;
							pc[(i + 1) * n + 1] += s11;
						}
					}
				}
			}
		}

#pragma omp barrier
	}

}
#undef block_dim
}

void
mat_mult_mat1(double * C, const double * A, const double * B, int n)
{
#define block_dim 256
	int blocks = (n + block_dim - 1) / block_dim;

#pragma omp parallel
{
	double As[block_dim][block_dim];
	double Bs[block_dim][block_dim];

	int id  = get_my_id();
	int num = get_num_threads();

	memset(As, 0, block_dim * block_dim * sizeof(double));
	memset(Bs, 0, block_dim * block_dim * sizeof(double));

#pragma omp for
	for (int i = 0; i < n * n; ++i)
	{
		C[i] = 0;
	}

	for (int l = id; l < blocks; l += num )
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

				double * pb = &Bs[0][0];
				double * pa = &As[0][0];

				for (int i = 0; i < nlm; ++i)
				{
					const double * bb = &B[(i + fk) * n + fm];
					for (int j = 0; j < nlk; ++j)
					{
						*pb++ = *bb++;
					}
					
				}

				for (int i = 0; i < nll; ++i)
				{
					const double * aa = &A[(i + fl) * n + fk];
					for (int j = 0; j < nlk; ++j)
					{
						*pa++ = *aa++;
					}
				}

				for (int i = 0; i < nll; ++i)
				{
					pb = &Bs[0][0];
					double * pc = &C[(i + fl) * n + fm];
					for (int j = 0; j < nlm; ++j)
					{
						double s = 0.0;
						pa = &As[i][0];
						for (int k1 = 0; k1 < nlm; ++k1)
						{
							// A[i][k] * B[k][j]
							s += *pa++ * *pb++;
						}

						//C[(i + fl) * n + fm + j] += s;
						*pc ++ += s;
					}
				}
			}

#pragma omp barrier
		}
	}
}
#undef block_dim
}

#ifdef TEST
#include <vector>
#include <stdio.h>

using namespace std;
using namespace linal;

#define N 1024L
//#define CHECK

int main()
{
	int n = N;
	vector < double > A(n * n);
	vector < double > B(n * n);
	vector < double > C(n * n);
	vector < double > reference(n * n);

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			if (i > j) {
				A[i * n + j] = 1. / (1. + i + j);
				B[i * n + j] = (2. + i + j);
			} else {
				A[i * n + j] = -1. / (1. + i + j);
				B[i * n + j] = -(2. + i + j);
			}
		}
	}

	double t, seconds;
	t = get_full_time();
	//omp_set_num_threads(4);
	//mat_mult_mat(&C[0], &A[0], &B[0], n);
	//mat_mult_mat_cannon(&C[0], &A[0], &B[0], n);
	mat_mult_mat1(&C[0], &A[0], &B[0], n);

	seconds = (get_full_time() - t) / 100.0;
	printf("t=%lf, gflops=%lf\n", seconds, 
		1e-9 * 2.0 * (double)n * (double)n * (double)n / seconds);

	// check

#ifdef CHECK
	t = get_full_time();
	mat_mult_mat_stupid(&reference[0], &A[0], &B[0], n);
	seconds = (get_full_time() - t) / 100.0;
	printf("t=%lf, gflops=%lf\n", seconds, 
		1e-9 * 2.0 * (double)n * (double)n * (double)n / seconds);

	for (int i = 0; i < n * n; ++i) {
		if (fabs(reference[i] - C[i]) > 1e-7) {
			printf("error ! \n %.16le != %.16le\n", reference[i], C[i]);
			return -1;
		}
	}

	printf("test passed!\n");
#endif

	return 0;
}
#endif
