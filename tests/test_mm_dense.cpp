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

#include "linal.h"
#include "timer.h"

#include <vector>
#include <stdio.h>

using namespace std;
using namespace linal;

#define N 4096L
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
	mat_mult_mat(&C[0], &A[0], &B[0], n);
	//mat_mult_mat_cannon(&C[0], &A[0], &B[0], n);
	//mat_mult_mat1(&C[0], &A[0], &B[0], n);

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

