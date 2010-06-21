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

#include "linal.h"

namespace linal
{

double vec_find_min (const double * v, int n)
        {
                double mn = v[0];
                for (int i = 0; i < n; ++i)
                {
                        if (v[i] < mn) mn = v[i];
                }
                return mn;
        }

        double vec_find_max (const double * v, int n)
        {
                double mx = v[0];
                for (int i = 0; i < n; ++i)
                {
                        if (v[i] > mx) mx = v[i];
                }
                return mx;
        }

void set_num_threads (int threads)
{
#ifdef _OPENMP
	if (threads)
	{
		omp_set_num_threads (threads);
	}
#endif
}

int get_num_threads ()
{
#ifdef _OPENMP
	return omp_get_num_threads();
#else
	return 1;
#endif
}

int get_my_id ()
{
#ifdef _OPENMP
	return omp_get_thread_num();
#else
	return 1;
#endif
}

void mat_transpose (double * out, const double * in, int n, int m)
{
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < m; ++j)
		{
			out[j * n + i] = in[i * m + j];
		}
	}
}

void mat_transpose1 (double * out, const double * in, double k, int n, int m)
{
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < m; ++j)
		{
			out[j * n + i] = k * in[i * m + j];
		}
	}
}

/**
 * Gauss
 * Copyright (c) 2009 Andrey Kornev, Andrey Ivanchikov, Alexey Ozeritsky
 * from manipro
 */
static void gauss_reverse (const double *A, const double *b, double *x, int n, int diagUnit)
{
	int j, k;

	for (k = n - 1; k >= 0; k--)
	{
		x[k] = b[k];
		for (j = k + 1; j < n; j++)
		{
			x[k] = x[k] - x[j] * A[k*n+j];
		}
		if (!diagUnit)
		{
			x[k] = x[k] / A[k*n+k];
		}
	}
}

int gauss (double *A, double *b, double *x, int n)
{
	int i, j, k;
	double p;
	int imax;
	double Eps = 1.e-15;

	for (k = 0; k < n; k++)
	{
		imax = k;

		for (i = k + 1; i < n; i++)
		{
			if (fabs (A[i*n+k]) > fabs (A[imax*n+k]) ) imax = i;
		}

		for (j = k; j < n; j++)
		{
			p = A[imax*n+j];
			A[imax*n+j] = A[k*n+j];
			A[k*n+j] = p;
		}
		p = b[imax];
		b[imax] = b[k];
		b[k] = p;

		p = A[k*n+k];

		if (fabs (p) < Eps)
		{
			printf ("Warning in %s %s : Near-null zero element\n", __FILE__, __FUNCTION__);
			return -1;
		}

		for (j = k; j < n; j++)
		{
			A[k*n+j] = A[k*n+j] / p;
		}
		b[k] = b[k] / p;

		for (i = k + 1; i < n; i++)
		{
			p = A[i*n+k];
			for (j = k; j < n; j++)
			{
				A[i*n+j] = A[i*n+j] - A[k*n+j] * p;
			}
			b[i] = b[i] - b[k] * p;
		}
	}

	gauss_reverse (A, b, x, n, true);

	return 0;
}

template < typename T >
void mat_print_ (FILE * f, const T * A, int n, int m, const char * fmt)
{
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			fprintf (f, fmt, A[i * n + j]);
		}
		fprintf (f, "\n");
	}
}

void mat_print (FILE * f, const float * A, int n, int m, const char * fmt)
{
	mat_print_ (f, A, n, m, fmt);
}

void mat_print (FILE * f, const double * A, int n, int m, const char * fmt)
{
	mat_print_ (f, A, n, m, fmt);
}

void mat_print (const char * fn, const float * A, int n, int m, const char * fmt)
{
	FILE * f = fopen (fn, "wb");
	if (!f) return;
	mat_print_ (f, A, n, m, fmt);
	fclose (f);
}

void mat_print (const char * fn, const double * A, int n, int m, const char * fmt)
{
	FILE * f = fopen (fn, "wb");
	if (!f) return;
	mat_print_ (f, A, n, m, fmt);
	fclose (f);
}

double vec_norm2 (const double * v, int n)
{
	return sqrt (vec_scalar2 (v, v, n) );
}

float vec_norm2 (const float * v, int n)
{
	float f = vec_scalar2 (v, v, n);
	return sqrtf (f);
}

template < typename T >
void vec_print_ (FILE * f, const T * A, int n, const char * fmt)
{
	for (int i = 0; i < n; ++i)
	{
		fprintf (f, fmt, A[i]);
	}
	fprintf (f, "\n");
}

void vec_print (FILE * f, const double * A, int n, const char * fmt)
{
	vec_print_ (f, A, n, fmt);
}

void vec_print (FILE * f, const float * A, int n, const char * fmt)
{
	vec_print_ (f, A, n, fmt);
}

void vec_print (const char * fn, const double * A, int n, const char * fmt)
{
	FILE * f = fopen (fn, "wb");
	if (!f) return;
	vec_print_ (f, A, n, fmt);
	fclose (f);
}

void vec_print (const char * fn, const float * A, int n, const char * fmt)
{
	FILE * f = fopen (fn, "wb");
	if (!f) return;
	vec_print_ (f, A, n, fmt);
	fclose (f);
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

template < typename T >
void csr_print_ (const int * Ap, const int * Ai,
                 const T * Ax, int n, FILE * f)
{
	int i, i0, j, k, i_old;
	const T * p = Ax;
	for (j = 0; j < n; ++j)
	{
		i_old = -1;
		for (i0 = Ap[j]; i0 < Ap[j + 1]; ++i0, ++p)
		{
			i = Ai[i0];
			for (k = i_old; k < i - 1; ++k)
			{
				fprintf (f, "%8.3le ", 0.0);
			}
			fprintf (f, "%8.3le ", (double) *p);
			i_old = i;
		}

		for (k = i_old + 1; k < n; ++k)
		{
			fprintf (f, "%8.3le ", 0.0);
		}
		fprintf (f, "\n");
	}
}

void sparse_print (const int * Ap, const int * Ai,
                   const double * Ax, int n, FILE * f)
{
	csr_print_ (Ap, Ai, Ax, n, f);
}

void sparse_print (const int * Ap, const int * Ai,
                   const float * Ax, int n, FILE * f)
{
	csr_print_ (Ap, Ai, Ax, n, f);
}

} /* namespace */
