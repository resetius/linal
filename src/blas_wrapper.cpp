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

	double ddot_(long *n, const double *x, long *incx, const double *y, long *incy);
	float sdot_(long *n, const float *x, long *incx, const float *y, long *incy);

	void dscal_(long *n1, const double * k, double * a, long * one);
	void sscal_(long *n1, const float * k, float * a, long * one);

	void dgemv_(char *trans, long *, long *, const double *alpha, const double *A, long *, 
			const double * x, long *, const double *beta, double *y, long *one);
	void sgemv_(char *trans, long *, long *, const float *alpha, const float *A, long *, 
			const float * x, long *, const float *beta, float *y, long *one);

	void dgemm_ (char *, char *, long *, long *, long *, const double *beta, 
			const double *A, long *, const double *B, long *, const double *alpha, 
			double *C, long *);
	void sgemm_ (char *, char *, long *, long *, long *, const float *beta, 
			const float *A, long *, const float *B, long *, const float *alpha, 
			float *C, long *);
}

namespace linal
{

/**
 * r = a + k2 * b
 */
void vec_sum2 (double * r, const double * a, const double *b, double k2, int n)
{
	long n1 = n;
	long one = 1;
	double alpha  = k2;
	double alpha1 = 1.0;
	if (r == a) {
		daxpy_(&n1, &alpha, b, &one, r, &one);
	} else if (r == b) {
		dscal_(&n1, &k2, r, &one);
		daxpy_(&n1, &alpha1, a, &one, r, &one);
	} else {
		dcopy_(&n1, a, &one, r, &one);
		vec_sum2(r, r, b, k2, n);
	}
}

void vec_sum2 (float * r, const float * a, const float *b, float k2, int n)
{
	long n1 = n;
	long one = 1;
	float alpha  = k2;
	float alpha1 = 1.0;
	if (r == a) {
		saxpy_(&n1, &alpha, b, &one, r, &one);
	} else if (r == b) {
		sscal_(&n1, &k2, r, &one);
		saxpy_(&n1, &alpha1, a, &one, r, &one);
	} else {
		scopy_(&n1, a, &one, r, &one);
		vec_sum2(r, r, b, k2, n);
	}
}

/**
 * r = k1 * a + k2 * b
 */
void vec_sum (double * r, const double * a, const double *b, int n)
{
	vec_sum2(r, a, b, 1.0, n);
}

void vec_sum (float * r, const float * a, const float *b, int n)
{
	vec_sum2(r, a, b, 1.0, n);
}

/**
 * a = b * k
 */
void vec_mult_scalar (double * a, const double * b, double k, int n)
{
	long n1  = n;
	long one = 1;
	if (a == b) {
		dscal_(&n1, &k, a, &one);
	} else {
		dcopy_(&n1, b, &one, a, &one);
		vec_mult_scalar(a, a, k, n);
	}
}

void vec_mult_scalar (float * a, const float * b, float k, int n)
{
	long n1  = n;
	long one = 1;
	if (a == b) {
		sscal_(&n1, &k, a, &one);
	} else {
		scopy_(&n1, b, &one, a, &one);
		vec_mult_scalar(a, a, k, n);
	}
}

/**
 * r = k1 * a + k2 * b
 */
void vec_sum1 (double * r, const double * a, const double *b, double k1, double k2, int n)
{
	if (r == a && r == b) {
		long n1  = n;
		long one = 1;
		double k = k1 + k2;
		dscal_(&n1, &k, r, &one);
	} else if (r == b) {
		vec_sum1(r, b, a, k2, k1, n);
	} else {
		vec_mult_scalar(r, a, k1, n);
		vec_sum2(r, r, b, k2, n);
	}
}

void vec_sum1 (float * r, const float * a, const float *b, float k1, float k2, int n)
{
	if (r == a && r == b) {
		long n1  = n;
		long one = 1;
		float k  = k1 + k2;
		sscal_(&n1, &k, r, &one);
	} if (r == b) {
		vec_sum1(r, b, a, k2, k1, n);
	} else {
		vec_mult_scalar(r, a, k1, n);
		vec_sum2(r, r, b, k2, n);
	}
}

/**
 * r = a - b
 */
void vec_diff (double * r, const double * a, const double * b, int n)
{
	vec_sum2(r, a, b, -1.0, n);
}

void vec_diff (float * r, const float * a, const float * b, int n)
{
	vec_sum2(r, a, b, -1.0, n);
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

/* level 2*/

void mat_mult_vector(double*y, double const*A, double const*x, int n)
{
	char trans = 'T';
	long n1 = n, one = 1;
	double alpha = 1.0, beta = 0.0;
	dgemv_(&trans, &n1, &n1, &alpha, A, &n1, x, &one, &beta, y, &one);
}

void mat_mult_vector(float*y, float const*A, float const*x, int n)
{
	char trans = 'T';
	long n1 = n, one = 1;
	float alpha = 1.0, beta = 0.0;
	sgemv_(&trans, &n1, &n1, &alpha, A, &n1, x, &one, &beta, y, &one);
}

/* level 3*/
void mat_mult_mat(double*C, double const*A, double const*B, int n1)
{
	long n = n1;
	char c = 'T';
	double alpha = 0.0;
	double beta = 1.0;
	dgemm_ (&c, &c, &n, &n, &n, &beta, A, &n, B, &n, &alpha, C, &n);
}

void mat_mult_mat(float*C, const float*A, const float*B, int n1)
{
	long n = n1;
	char c = 'T';
	float alpha = 0.0;
	float beta = 1.0;
	sgemm_ (&c, &c, &n, &n, &n, &beta, A, &n, B, &n, &alpha, C, &n);
}

} /* namespace linal */

