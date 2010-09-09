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

extern "C"
{
	void dcopy_(long * n, const double * x, long *incx, const double * y, long *incy);
	void scopy_(long * n, const float * x, long *incx, const float * y, long *incy);

	void daxpy_(long *n, const double * alpha, const double * x, long *incx, double *y, long *incy);
	void saxpy_(long *n, const float * alpha, const float * x, long *incx, float *y, long *incy);

	void ddot_(long *n, const double *x, long *incx, const double *y, long *incy);
	void sdot_(long *n, const float *x, long *incx, const float *y, long *incy);
}

namespace linal
{

/**
 * r = k1 * a + k2 * b
 */
void vec_sum1 (double * r, const double * a, const double *b, double k1, double k2, int n)
{
}

void vec_sum1 (float * r, const float * a, const float *b, float k1, float k2, int n)
{
}

/**
 * r = a + k2 * b
 */
void vec_sum2 (double * r, const double * a, const double *b, double k2, int n)
{
}

void vec_sum2 (float * r, const float * a, const float *b, float k2, int n)
{
}

void vec_sum (double * r, const double * a, const double *b, int n)
{
}

void vec_sum (float * r, const float * a, const float *b, int n)
{
}

/**
 * a = b * k
 */
void vec_mult_scalar (double * a, const double * b, double k, int n)
{
}

void vec_mult_scalar (float * a, const float * b, float k, int n)
{
}

void vec_mult (double * r, const double * a, const double *b, int n)
{
}

void vec_mult (float * r, const float * a, const float *b, int n)
{
}

/**
 * r = a - b
 */
void vec_diff (double * r, const double * a, const double * b, int n)
{
	long n1 = n;
	long one = 1;
	double alpha = -1.0;
	if (r == a) {
		daxpy_(&n1, &alpha, b, &one, r, &one);
	} else {
		dcopy_(&n1, a, &one, r, &one);
		vec_diff(r, r, b, n);
	}
}

void vec_diff (float * r, const float * a, const float * b, int n)
{
	long n1 = n;
	long one = 1;
	float alpha = -1.0;

	if (r == a) {
		saxpy_(&n1, &alpha, b, &one, r, &one);
	} else {
		scopy_(&n1, a, &one, r, &one);
		vec_diff(r, r, b, n);
	}
}

double vec_scalar2 (const double * a, const double * b, int n)
{
	long one = 1;
	long n1  = n;
	return ddot_(&n1, a, &one, b, &one);
}

float vec_scalar2 (const float * a, const float * b, int n)
{
	long one = 1;
	long n1  = n;
	return sdot_(&n1, a, &one, b, &one);
}

} /* namespace linal */

