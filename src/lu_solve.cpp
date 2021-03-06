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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "lu_solve.h"

namespace linal {

static void reverse_upper (double * x, const double *A, const double *b, int n)
{
	int j, k;

	for (k = n - 1; k >= 0; --k)
	{
		x[k] = b[k];
		for (j = k + 1; j < n; j++)
		{
			x[k] -= x[j] * A[k * n + j];
		}
		x[k] = x[k] / A[k * n + k];
	}
}

static void reverse_lower (double * x, const double *A, const double *b, int n)
{
	int j, k;

	for (k = 0; k < n; ++k)
	{
		x[k] = b[k];
		for (j = k - 1; j >= 0; --j)
		{
			x[k] -= x[j] * A[k * n + j];
		}
		x[k] = x[k] / A[k * n + k];
	}
}

void lu_solve(double * x, const double * L, const double * U, const double * b, int n)
{
	reverse_lower(x, L, b, n);
	reverse_upper(x, U, x, n);
}

void lu_build (double * L, double * U, const double * A, int n)
{
	int i, j, k;
	memcpy (U, A, n * n * sizeof (double) );
	memset (L, 0, n * n * sizeof (double) );
	for (i = 0; i < n; ++i)
	{
		L[i * n + i] = 1.0;
	}

	for (j = 0; j < n; j++)
	{
		double max = U[j * n + j];
		int r  = j;
		for (k = j; k < n; ++k)
		{
			double c = fabs (U[k * n + j]);
			if (max < c)
			{
				r   = k;
				max = c;
			}
		}

		if (r != j)
		{
			double temp;
			for (k = j; k < n; k++)
			{
				temp         = U[j * n + k];
				U[j * n + k] = U[r * n + k];
				U[r * n + k] = temp;
			}

			for (k = 0; k < n; k++)
			{
				temp         = L[j * n + k];
				L[j * n + k] = L[r * n + k];
				L[r * n + k] = temp;
			}
		}

		for (i = j + 1; i < n; i++)
		{
			double ba = U[i * n + j] / U[j * n + j];

			for (k = j; k < n; ++k)
			{
				U[i * n + k] -= U[j * n + k] * ba;
			}

			for (k = 0; k < n; ++k)
			{
				L[i * n + k] -= L[j * n + k] * ba;
			}
		}
	}

	for (k = 0; k < n; ++k)
	{
		for (i = 0; i < k; ++i) {
			double s = 0.0;
			for (j = k - 1; j >= 0; --j)
			{
				s -= L[j * n + i] * L[k * n + j];
			}
			L[k * n + i] = s / L[k * n + k];
		}
	}
}

} /* namespace linal */

