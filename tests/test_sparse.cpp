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
[--task] [--dim] [--iters] --file file \n", n);
	fprintf (stderr, "--threads n - sets the number of threads\n");
	fprintf (stderr, "--file file - sparse matrix file\n");
	fprintf (stderr, "--task - all\n"
	         "mult\n"
	         "invert\n");
	exit (-1);
}

template < typename Solver >
void do_test_mult_sparse(const Solver & solver, FILE * f, int iters)
{
	int n = solver.dim();
	int nz = solver.nonzero();
	vector < typename Solver::data_type > rp(n);
	vector < typename Solver::data_type > ans(n);

	fprintf(stderr, "mult: n:%d, nz:%d\n", n, nz);

	Timer t;
	for (int i = 0; i < iters; ++i) {
		solver.mult_vector(&ans[0], &rp[0]);
	}

	double seconds = t.elapsed();
	double gflops  = 1e-9 * 2.0 * (double)iters * (double)nz / seconds;
	fprintf(stderr, "mult: %.16lfs, %.16lfgflops\n", seconds, gflops);
}

template < typename Solver >
void do_test_invert_sparse(const Solver & solver, FILE * f, int iters)
{
	int n = solver.dim();
	int nz = solver.nonzero();
	vector < typename Solver::data_type > rp(n);
	vector < typename Solver::data_type > ans(n);

	fprintf(stderr, "invert: n:%d, nz:%d\n", n, nz);
	Timer t;
	for (int i = 0; i < 1; ++i) {
		rp[i] = 1.0;
		solver.solve(&ans[0], &rp[0]);
	}

#pragma omp parallel for
	for (int i = 1; i < n; ++i) {
		rp[i] = 1.0;
		solver.solve(&ans[0], &rp[0]);
	}

	double seconds = t.elapsed();
	fprintf(stderr, "invert: %.16lfs\n", seconds);
}

template < typename Store >
void do_test_sparse(FILE * f, int iters, bool mult, bool invert)
{
	Store store;
	store.restore(f);

	if (mult) {
		do_test_mult_sparse(make_sparse_solver(store, store), f, iters);
	}

	if (invert) {
#ifdef UMFPACK
		do_test_invert_sparse(make_umfpack_solver(store, store), f, iters);
#else
		do_test_invert_sparse(make_sparse_solver(store, store), f, iters);
#endif
	}
}

bool test_sparse (FILE * f, int iters, bool mult, bool invert)
{
	char tag[4];
	int size;

	fread(tag, 4, 1, f);
	fread(&size, 4, 1, f);
	fseek(f, 0, SEEK_SET);

	if (!memcmp(tag, "CSR ", 4)) {
		if (size == 4) {
			do_test_sparse < StoreCSR < float > > (f, iters, mult, invert);
		} else if (size == 8 && check_device_supports_double()) {
			do_test_sparse < StoreCSR < double > > (f, iters, mult, invert);
		} else {
			throw runtime_error("unsupported floating point format\n");
		}
	} else if (!memcmp(tag, "ELL ", 4)) {
		if (size == 4) {
			do_test_sparse < StoreELL < float > > (f, iters, mult, invert);
		} else if (size == 8 && check_device_supports_double()) {
			do_test_sparse < StoreELL < double > > (f, iters, mult, invert);
		} else {
			throw runtime_error("unsupported floating point format\n");
		}
	} else {
		throw runtime_error("unknown format\n");
	}

	return true;
}

int main (int argc, char * argv[])
{
	bool result = true;
	bool mult   = false;
	bool invert = false;
	FILE * f = 0;

	int dim   = 0;
	int iters = 1000;

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

			if (!strcmp (argv[i + 1], "mult") )
			{
				mult   = true;
			}
			if (!strcmp (argv[i + 1], "invert") )
			{
				invert = true;
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
		if (!strcmp (argv[i], "--file"))
		{
			if (i == argc - 1)
			{
				usage (argv[0]);
			}

			f = fopen (argv[i + 1], "rb");
			if (!f) {
				usage (argv[0]);
			}
		}
		if (!strcmp (argv[i], "--help") || !strcmp (argv[i], "-h") )
		{
			usage (argv[0]);
		}
	}

	if (!f) {
		usage(argv[0]);
	}

	fprintf(stderr, "mult: %d, invert: %d\n", (int)mult, (int)invert);
	linal_init();

//	try
//	{
		Timer t;
		//burn (1.0);
		t.restart();
		result &= test_sparse (f, iters, mult, invert);
		fprintf (stderr, "elapsed: %lf\n", t.elapsed() );
//	}
//	catch (const std::exception & e)
//	{
//		fprintf (stderr, "exception: %s\n", e.what() );
//	}

	linal_shutdown();
	if (f) fclose(f);
	return (int) (!result);
}

