/* -*- charset: utf-8 -*- */
/*$Id$*/

/* Copyright (c) 2010-2015 Alexey Ozeritsky
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

#include <stdio.h>
#include <stdlib.h>
#include "lu_solve.h"

#define N 10

using namespace linal;

static void mat_print(const double * A, int n)
{
	int i, j;
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < n; ++j)
		{
			fprintf(stderr, "%8.3lf ", A[i * n + j]);
		}
		fprintf(stderr, "\n");
	}
}

int test_lu(int argc, char * argv[])
{
	double * A = (double *)malloc(N * N * sizeof(double));
	double * L = (double *)malloc(N * N * sizeof(double));
	double * U = (double *)malloc(N * N * sizeof(double));
	int i, j, k;
	int n = N;

	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < n; ++j)
		{
			A[i * n + j] = 1.0 / (1.0 + i + j);
			if (i > j)
			{
				A[i * n + j] *= -1.0;
			}
		}
	}

	fprintf(stderr, "A:\n");
	mat_print(A, n);
	lu_build(L, U, A, n);
	fprintf(stderr, "L:\n");
	mat_print(L, n);
	fprintf(stderr, "U:\n");
	mat_print(U, n);

	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < n; ++j)
		{
			double s = 0;

			for (k = 0; k < n; ++k)
			{
				s += L[i * n + k] * U[k * n + j];
			}

			A[i * n + j] = s;
		}
	}
	fprintf(stderr, "LU:\n");
	mat_print(A, n);

	free(A);
	free(L);
	free(U);

	return 1;
}
