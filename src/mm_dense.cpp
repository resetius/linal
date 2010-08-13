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

#include "mm_dense.h"

namespace linal
{

static void set_num_threads (int threads)
{
#ifdef _OPENMP
	if (threads)
	{
		omp_set_num_threads (threads);
	}
#endif
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

template < typename T >
void mat_mult_vector_ (T * r, const T * A, const T * x, int n)
{
#pragma omp parallel
	{
		int procs  = get_num_threads();
		int blocks = procs;
		int id     = get_my_id();
		int f_row  = n * id;
		f_row /= procs;
		int l_row  = n * (id + 1);
		l_row = l_row / procs - 1;
		int m = id;

		for (int i = f_row; i <= l_row; ++i)
		{
			r[i] = 0;
		}

		for (int block = 0; block < blocks; ++block, m = (m + 1) % blocks)
		{
			int fm = n * m;
			fm /= blocks;
			int lm = n * (m + 1);
			lm = lm / blocks - 1;

			for (int i = f_row; i <= l_row; ++i)
			{
				T s = 0.0;
				const T * ax = &A[i * n + fm];
				const T * xx = &x[fm];

				for (int j = fm; j <= lm; ++j)
				{
					s += *ax++ * *xx++;
				}
				r[i] += s;
			}
		}
	}
}

void mat_mult_vector (double * r, const double * A, const double * x, int n)
{
	mat_mult_vector_ (r, A, x, n);
}

void mat_mult_vector (float * r, const float * A, const float * x, int n)
{
	mat_mult_vector_ (r, A, x, n);
}

template < typename T >
void mat_mult_vector_stupid_ (T * r, const T * A, const T * x, int n)
{
#pragma omp parallel for
	for (int i = 0; i < n; ++i)
	{
		T s = 0.0;
		for (int j = 0; j < n; ++j)
		{
			s += A[i * n + j] * x[j];
		}
		r[i] = s;
	}
}

void mat_mult_vector_stupid (double * r, const double * A, const double * x, int n)
{
	mat_mult_vector_stupid_ (r, A, x, n);
}

void mat_mult_vector_stupid (float * r, const float * A, const float * x, int n)
{
	mat_mult_vector_stupid_ (r, A, x, n);
}

} /* namespace */

