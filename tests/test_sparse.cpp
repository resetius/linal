/* -*- charset: utf-8 -*- */

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <vector>
#include <string>

#include "gmres.h"
#include "solver.h"

using namespace linal;
using namespace std;

void usage (const char * n)
{
	fprintf (stderr, "%s [--threads|-t n] \
[--task] [--dim] [--iters]\n", n);
	fprintf (stderr, "--threads n - sets the number of threads\n");
	fprintf (stderr, "--file file - sparse matrix file\n");
	fprintf (stderr, "--task - all\n"
	         "mult_sparse\n"
	         "invert_sparse\n");
	exit (-1);
}

bool check (float val)
{
	return fabs (val) < 1e-5;
}

bool check (double val)
{
	return fabs (val) < 1e-12;
}

template < typename T >
bool test_solve (int n, int iters)
{
	int i, j = 0;
	if (n <= 0)     n     = 320000;
	if (iters <= 0) iters = 1;

	bool ret = true;

	vector < T > b (n);
	vector < T > x1 (n);
	vector < T > x2 (n);
	vector < T > x3 (n);
	vector < T > v (n);
	fprintf (stderr, "n=%d, iters=%d\n", n, iters);

	/* матрицу записываем по строкам! */
	SparseSolver  < T, StoreELL < T , Allocator > , StoreELL < T , Allocator > > M1 (n);
	init_matrix (M1, n);

#ifdef UMFPACK
	UmfPackSolver < T, StoreCSR < T , Allocator > > M2 (n);
	init_matrix (M2, n);
#endif

#ifdef SUPERLU
	SuperLUSolver < T, StoreCSR < T , Allocator > > M3 (n);
	init_matrix (M1, n);
#endif

	for (i = 0; i < n; ++i)
	{
		b[i] = 1;
	}

	Timer t;

	t.restart();
	for (int k = 0; k < iters; ++k)
	{
		M1.solve (&x1[0], &b[0]);
	}
	fprintf (stdout, "gmres solve: %lf\n", t.elapsed() );
#ifdef UMFPACK
	t.restart();
	for (int k = 0; k < iters; ++k)
	{
		M2.solve (&x2[0], &b[0]);
	}
	fprintf (stdout, "umfpack solve: %lf\n", t.elapsed() );
#endif
#ifdef SUPERLU
	t.restart();
	for (int k = 0; k < iters; ++k)
	{
		M3.solve (&x3[0], &b[0]);
	}
	fprintf (stdout, "superlu solve: %lf\n", t.elapsed() );
#endif
	T nr;
#if defined(UMFPACK) && defined(SUPERLU)
	vec_diff (&v[0], &x2[0], &x3[0], (int) x1.size() );
	nr = vec_norm2 (&v[0], (int) x2.size() );
	fprintf (stderr, "%.16le\n", (double) nr);
	ret &= check (nr);
#endif
	vec_diff (&v[0], &x1[0], &x2[0], (int) x1.size() );
	nr = vec_norm2 (&v[0], (int) x2.size() );
	fprintf (stderr, "%.16le\n", (double) nr);
//	ret &= check(nr);

	return ret;
}

template < typename T >
struct MatrixLoader
{
};

template < typename M, typename V >
void mult_loop(V & x1, M & M1, V & b, int iters)
{
	Timer t;
	t.restart();
	for (int k = 0; k < iters; ++k)
	{
		M1.mult_vector (&x1[0], &b[0]);
	}
	fprintf (stdout, "M1 mult_vector: %lf\n", t.elapsed() );
}

template < typename T >
bool test_mult_sparse (FILE * f, int iters)
{
	int i, j = 0;
	if (iters <= 0) iters = 10000;

	bool ret = true;

	vector < T > b (n);
	vector < T > x1 (n);
	vector < T > x2 (n);
	vector < T > x3 (n);
	vector < T > v (n);
	fprintf (stderr, "n=%d, iters=%d\n", n, iters);

	for (i = 0; i < n; ++i)
	{
		b[i] = 1;
	}

	/* матрицу записываем по строкам! */
	SparseSolver  < T, StoreELL < T , Allocator > , StoreELL < T , Allocator > > M1 (n);
	init_matrix (M1, n);
	mult_loop(x1, M1, b, iters);

	return ret;
}

int main (int argc, char * argv[])
{
	bool result        = true;
	bool mult_sparse   = false;
	bool invert_sparse = false;

	int dim   = 0;
	int iters = 0;

	for (int i = 0; i < argc; ++i)
	{
		if (!strcmp (argv[i], "--threads") || !strcmp (argv[i], "-t") )
		{
			if (i == argc - 1)
			{
				usage (argv[0]);
			}

			int threads = atoi (argv[i + 1]);

			fprintf(stderr, "threads=%d\n", threads);
			set_num_threads (threads);
		}
		if (!strcmp (argv[i], "--task") )
		{
			if (i == argc - 1)
			{
				usage (argv[0]);
			}

			if (!strcmp (argv[i + 1], "mult_sparse") )
			{
				mult_sparse   = true;
			}
			if (!strcmp (argv[i + 1], "invert_sparse") )
			{
				invert_sparse = true;
			}
		}
		if (!strcmp (argv[i], "--iters") )
		{
			if (i == argc - 1)
			{
				usage (argv[0]);
			}

			iters = atoi (argv[i + 1]);
		}
		if (!strcmp (argv[i], "--help") || !strcmp (argv[i], "-h") )
		{
			usage (argv[0]);
		}
	}

	linal_init();

	try
	{
		int has_double = check_device_supports_double();
		fprintf (stderr, "has double: %d\n", has_double);

		Timer t;
		//burn (1.0);
		t.restart();
		if (use_double && has_double)
		{
			fprintf (stderr, "testing double:\n");

			if (invert_sparse)
			{
				t.restart();
				result &= test_solve < double > (dim, iters);
				fprintf (stderr, "test_solve < double > (): %lf, %d\n", t.elapsed(), (int) result);
			}

			if (mult_sparse)
			{
				t.restart();
				result &= test_mult < double > (dim, iters);
				fprintf (stderr, "test_mult < double > (): %lf, %d\n", t.elapsed(), (int) result);
			}
		}
		else
		{
			fprintf (stderr, "testing float:\n");

			if (invert_sparse)
			{
				t.restart();
				result &= test_solve < float > (dim, iters);
				fprintf (stderr, "test_solve < float > (): %lf, %d\n", t.elapsed(), (int) result);
			}

			if (mult_sparse)
			{
				t.restart();
				result &= test_mult < float > (dim, iters);
				fprintf (stderr, "test_mult < float > (): %lf, %d\n", t.elapsed(), (int) result);
			}

		}
		fprintf (stderr, "elapsed: %lf\n", t.elapsed() );
	}
	catch (const std::exception & e)
	{
		fprintf (stderr, "exception: %s\n", e.what() );
	}

	linal_shutdown();
	return (int) (!result);
}

