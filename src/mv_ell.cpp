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
void ell_mult_vector_ (
    T * r,
    const int * Ai,
    const T * Ax,
    const T * x,
    int n,
    int cols,
    int stride)
{
#pragma omp parallel for
	for (int row = 0; row < n; ++row)
	{
		T sum = 0;
		for (int i0 = 0; i0 < cols; i0++)
		{
			const T A_ij = Ax[stride * i0 + row];
			const int col = Ai[stride * i0 + row];
			sum += A_ij * x[col];
		}
		r[row] = sum;
	}
}

void ell_mult_vector_r (double * r, const int * Ai,
                        const double * Ax, const double * x,
                        int n, int cols, int stride)
{
	ell_mult_vector_ (r, Ai, Ax, x, n, cols, stride);
}

void ell_mult_vector_r (float * r, const int * Ai,
                        const float * Ax, const float * x,
                        int n, int cols, int stride)
{
	ell_mult_vector_ (r, Ai, Ax, x, n, cols, stride);
}

} /* namespace linal */

