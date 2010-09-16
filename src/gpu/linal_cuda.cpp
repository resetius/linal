/* -*- charset: utf-8 -*- */
/*$Id$*/

/* Copyright (c) 2009-2010 Alexey Ozeritsky
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

#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include <cublas.h>
#include <cuda_runtime_api.h>

#include "linal.h"
#include "linal_cuda.h"

namespace linal {

double vec_scalar2(const double * a, const double * b, int n)
{
	return cublasDdot(n, a, 1, b, 1);
}

/*
float vec_scalar2(const float * a, const float * b, int n)
{
	return cublasSdot(n, a, 1, b, 1);
}
*/

/* level 3*/
void mat_mult_mat(double*C, double const*A, double const*B, int n1)
{
	int n   = n1;
	double alpha = 0.0;
	double beta = 1.0;
	cublasDgemm ('N', 'N', n, n, n, beta, A, n, B, n, alpha, C, n);
/*
	int status = cublasGetError();
	if (status != CUBLAS_STATUS_SUCCESS)
	{
		exit(1);
	}
	*/
}

void mat_mult_mat(float*C, const float*A, const float*B, int n1)
{
	long n = n1;
	char c = 'T';
	float alpha = 0.0;
	float beta = 1.0;
	cublasSgemm (c, c, n, n, n, beta, A, n, B, n, alpha, C, n);
/*
	int status = cublasGetError();
	if (status != CUBLAS_STATUS_SUCCESS)
	{
		exit(1);
	}
	*/
}

int check_device_supports_double()
{
	double pi1 = M_PI;
	double pi2 = M_PI;
	cudaSetDoubleForDevice(&pi1);
	return pi1 == pi2;
}

void linal_init()
{
	if (cublasInit() != CUBLAS_STATUS_SUCCESS)
	{
		fprintf(stderr, "failed to initialize library\n");
		exit(1);
	} else {
		fprintf(stderr, "library initialized\n");
	}

	int device;
	struct cudaDeviceProp prop;
	cudaGetDevice(&device);
	cudaGetDeviceProperties(&prop, device);
	fprintf(stderr, "name: %s\n", prop.name);
	fprintf(stderr, "version: %d.%d\n", prop.major, prop.minor);
	fprintf(stderr, "totalGlobalMem: %lu\n", prop.totalGlobalMem);

	fprintf(stderr, "maxthreads: %d\n", prop.maxThreadsPerBlock);
	fprintf(stderr, "maxthreadsdim: %d %d %d\n", prop.maxThreadsDim[0], prop.maxThreadsDim[1], prop.maxThreadsDim[2]);
	fprintf(stderr, "maxgridsize: %d %d %d\n", prop.maxGridSize[0], prop.maxGridSize[1], prop.maxGridSize[2]);
}

void linal_shutdown()
{
	cublasShutdown();
}

void linal_sync()
{
	cudaThreadSynchronize();
}

} /* namespace phelm */

