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

/**
 * r = k1 * a + k2 * b
 */
template < typename T >
void vec_sum1_ (T * r, const T * a, const T *b, T k1, T k2, int n)
{
	int i;
#pragma omp parallel for
	for (i = 0; i < n; ++i)
	{
		r[i] = k1 * a[i] + k2 * b[i];
	}
}

void vec_sum1 (double * r, const double * a, const double *b, double k1, double k2, int n)
{
	vec_sum1_ (r, a, b, k1, k2, n);
}

void vec_sum1 (float * r, const float * a, const float *b, float k1, float k2, int n)
{
	vec_sum1_ (r, a, b, k1, k2, n);
}

/**
 * r = a + k2 * b
 */
template < typename T >
void vec_sum2_ (T * r, const T * a, const T *b, T k2, int n)
{
	int i;
#pragma omp parallel for
	for (i = 0; i < n; ++i)
	{
		r[i] = a[i] + k2 * b[i];
	}
}

void vec_sum2 (double * r, const double * a, const double *b, double k2, int n)
{
	vec_sum2_ (r, a, b, k2, n);
}

void vec_sum2 (float * r, const float * a, const float *b, float k2, int n)
{
	vec_sum2_ (r, a, b, k2, n);
}

template < typename T >
void vec_sum_ (T * r, const T * a, const T *b, int n)
{
	int i;
#pragma omp parallel for
	for (i = 0; i < n; ++i)
	{
		r[i] = a[i] + b[i];
	}
}

void vec_sum (double * r, const double * a, const double *b, int n)
{
	vec_sum_ (r, a, b, n);
}

void vec_sum (float * r, const float * a, const float *b, int n)
{
	vec_sum_ (r, a, b, n);
}

/**
 * a = b * k
 */
template < typename T >
void vec_mult_scalar_ (T * a, const T * b, T k, int n)
{
	int i;
#pragma omp parallel for
	for (i = 0; i < n; ++i)
	{
		a[i] = b[i] * k;
	}
}

void vec_mult_scalar (double * a, const double * b, double k, int n)
{
	vec_mult_scalar_ (a, b, k, n);
}

void vec_mult_scalar (float * a, const float * b, float k, int n)
{
	vec_mult_scalar_ (a, b, k, n);
}

template < typename T >
void vec_mult_ (T * r, const T * a, const T *b, int n)
{
	int i;
#pragma omp parallel for
	for (i = 0; i < n; ++i)
	{
		r[i] = a[i] * b[i];
	}
}

void vec_mult (double * r, const double * a, const double *b, int n)
{
	vec_mult_ (r, a, b, n);
}

void vec_mult (float * r, const float * a, const float *b, int n)
{
	vec_mult_ (r, a, b, n);
}

/**
 * r = a - b
 */
template < typename T >
void vec_diff_ (T * r, const T * a, const T * b, int n)
{
	int i;
#pragma omp parallel for
	for (i = 0; i < n; ++i)
	{
		r[i] = a[i] - b[i];
	}
}

void vec_diff (double * r, const double * a, const double * b, int n)
{
	vec_diff_ (r, a, b, n);
}

void vec_diff (float * r, const float * a, const float * b, int n)
{
	vec_diff_ (r, a, b, n);
}

template < typename T >
T vec_scalar2_ (const T * a, const T * b, int n)
{
	T s = (T) 0.0;
	int i;
#pragma omp parallel for reduction(+:s)
	for (i = 0; i < n; ++i)
	{
		s = s + a[i] * b[i];
	}
	return s;
}

double vec_scalar2 (const double * a, const double * b, int n)
{
	return vec_scalar2_ (a, b, n);
}

float vec_scalar2 (const float * a, const float * b, int n)
{
	return vec_scalar2_ (a, b, n);
}

} /* namespace linal */

