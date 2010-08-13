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

#ifdef _OPENMP
#include <omp.h>
#endif

#include "timer.h"

extern "C"
void
mat_mult_mat(double * C, const double * A, const double * B, int n)
{
#define block_dim 128
	double As[block_dim][block_dim];
	double Bs[block_dim][block_dim];
	double Cs[block_dim][block_dim];

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

#pragma omp for schedule(static)
			for (int i = 0; i < nll; ++i)
			{
				memcpy(&Cs[i][0], &C[(i + fl) * n + fm], nlm * sizeof(double));
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
					memcpy(&As[i][0], &A[(i + fl) * n + fk], nlk * sizeof(double));
				}

#pragma omp for schedule(static)
				for (int i = 0; i < nlm; ++i)
				{
					memcpy(&Bs[i][0], &B[(i + fk) * n + fm], nlk * sizeof(double));
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
						Cs[i][j] += s;
					}
				}
#pragma omp barrier
			}

#pragma omp for schedule(static)
			for (int i = 0; i < nll; ++i)
			{
				memcpy(&C[(i + fl) * n + fm], &Cs[i][0], nlm * sizeof(double));
			}
#pragma omp barrier
		}
	}
}

}

#ifdef TEST
#include <vector>
#include <stdio.h>

using namespace std;
#define N 2048L

int main()
{
	vector < double > A(N * N);
	vector < double > B(N * N);
	vector < double > C(N * N);
	double t  = get_full_time();
	//omp_set_num_threads(4);
	mat_mult_mat(&C[0], &A[0], &B[0], N);
	printf("%lf\n", (get_full_time() - t) / 100.0);
}
#endif
