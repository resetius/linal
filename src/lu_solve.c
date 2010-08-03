#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

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

void build_lu (double * L, double * U, const double * A, int n)
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

#ifdef TEST
#define N 10

static void mat_print (const double * A, int n)
{
	int i, j;
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < n; ++j)
		{
			fprintf (stderr, "%8.3lf ", A[i * n + j]);
		}
		fprintf (stderr, "\n");
	}
}

int main()
{
	double * A = malloc (N * N * sizeof (double) );
	double * L = malloc (N * N * sizeof (double) );
	double * U = malloc (N * N * sizeof (double) );
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

	fprintf (stderr, "A:\n");
	mat_print (A, n);
	build_lu (L, U, A, n);
	fprintf (stderr, "L:\n");
	mat_print (L, n);
	fprintf (stderr, "U:\n");
	mat_print (U, n);

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
	fprintf (stderr, "LU:\n");
	mat_print (A, n);

	free (A);
	free (L);
	free (U);

	return 1;
}
#endif

