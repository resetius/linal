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

#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace linal
{

template < typename T >
void csr_add_matrix1_ (const int * oAp, T * oAx,
                       const int * Ap, const int * Ai, const T * Ax,
                       const T * x, int n)
{
#pragma omp parallel for
	for (int j = 0; j < n; ++j)
	{
		const T * a =  &Ax[Ap[j]];
		T * o = &oAx[oAp[j]];

		for (int i0 = Ap[j]; i0 < Ap[j + 1]; ++i0, ++a, ++o)
		{
			int i = Ai[i0];
			*o += *a * x[i];
		}
	}
}

void csr_add_matrix1 (const int * oAp, double * oAx,
                      const int * Ap, const int * Ai, const double * Ax,
                      const double * x, int n)
{
	csr_add_matrix1_ (oAp, oAx, Ap, Ai, Ax, x, n);
}

void csr_add_matrix1 (const int * oAp, float * oAx,
                      const int * Ap, const int * Ai, const float * Ax,
                      const float * x, int n)
{
	csr_add_matrix1_ (oAp, oAx, Ap, Ai, Ax, x, n);
}

template < typename T >
void csr_add_matrix2_ (const int * oAp, T * oAx,
                       const int * Ap, const int * Ai, const T * Ax,
                       const T * x, int n)
{
#pragma omp parallel for
	for (int j = 0; j < n; ++j)
	{
		const T * a =  &Ax[Ap[j]];
		T x1 = x[j];
		T * o = &oAx[oAp[j]];

		for (int i0 = Ap[j]; i0 < Ap[j + 1]; ++i0, ++a, ++o)
		{
			int i = Ai[i0];
			*o += *a * x1;
		}
	}
}

void csr_add_matrix2 (const int * oAp, double * oAx,
                      const int * Ap, const int * Ai, const double * Ax,
                      const double * x, int n)
{
	csr_add_matrix2_ (oAp, oAx, Ap, Ai, Ax, x, n);
}

void csr_add_matrix2 (const int * oAp, float * oAx,
                      const int * Ap, const int * Ai, const float * Ax,
                      const float * x, int n)
{
	csr_add_matrix2_ (oAp, oAx, Ap, Ai, Ax, x, n);
}

} /* namespace linal */

