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
//#define N 2L
//#define CHECK

static void
mat_mult_mat_stupid1(double * C, const double * A, const double * B, int n)
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
int do_all()
{
	int n = N;
	int iters = 1;

	linal_init();

	vector < T > Ah(n * n);
	vector < T > Bh(n * n);
	vector < T > Ch(n * n);
	ArrayDevice < T > A(n * n);
	ArrayDevice < T > B(n * n);
	ArrayDevice < T > C(n * n);
	vector < T > reference(n * n);

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			if (i > j) {
				Ah[i * n + j] = 1. / (1. + i + j);
				Bh[i * n + j] = (2. + i + j);
			} else {
				Ah[i * n + j] = -1. / (1. + i + j);
				Bh[i * n + j] = -(2. + i + j);
			}
		}
	}

	double t, seconds;
	vec_copy_from_host(&A[0], &Ah[0], n * n);
	vec_copy_from_host(&B[0], &Bh[0], n * n);
#ifdef GPGPU
	cudaThreadSynchronize();
#endif
	fprintf(stderr, "go!\n");

	t = get_full_time();
	//omp_set_num_threads(4);
	for (int i = 0; i < iters; ++i) {
		mat_mult_mat(&C[0], &A[0], &B[0], n);
	}
	//mat_mult_mat_cannon(&C[0], &A[0], &B[0], n);
	//mat_mult_mat1(&C[0], &A[0], &B[0], n);

#ifdef GPGPU
	cudaThreadSynchronize();
#endif
	seconds = (get_full_time() - t) / 100.0;

	vec_copy_from_device(&Ch[0], &C[0], n * n);
	//vec_copy_from_device(&Ch[0], &A[0], n * n);
	fprintf(stderr, "t=%lf, gflops=%lf\n", seconds, 
		1e-9 * 2.0 * (double)n * (double)n * (double)n / seconds * iters);

	// check

#ifdef CHECK
	t = get_full_time();
	mat_mult_mat_stupid1(&reference[0], &Ah[0], &Bh[0], n);
	seconds = (get_full_time() - t) / 100.0;
	fprintf(stderr, "t=%lf, gflops=%lf\n", seconds, 
		1e-9 * 2.0 * (double)n * (double)n * (double)n / seconds);

	for (int i = 0; i < n * n; ++i) {
		if (fabs(reference[i] - Ch[i]) > 1e-7) {
			fprintf(stderr, "error ! \n %.16le != %.16le\n", reference[i], Ch[i]);
			return -1;
		}
	}

	printf("test passed!\n");
#endif

	linal_shutdown();

	return 0;
}

int main()
{
	//do_all < double > ();
	do_all < float > ();
}
