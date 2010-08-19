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

template < typename T >
struct MatrixLoader
{
};

template < typename Store >
void do_test_mult_sparse(FILE * f, int iters)
{
	Store store;
	store.restore(f);
}

bool test_mult_sparse (FILE * f, int iters)
{
	char tag[4];
	int size;

	fread(tag, 4, 1, f);
	fread(&size, 4, 1, f);
	fseek(f, 0, SEEK_SET);

	if (!memcmp(tag, "CSR ", 4)) {
		if (size == 4) {
			do_test_mult_sparse < StoreCSR < float > > (f, iters);
		} else if (size == 8) {
			do_test_mult_sparse < StoreCSR < double > > (f, iters);
		} else {
			throw runtime_error("unknown format\n");
		}
	} else if (!memcmp(tag, "ELL ", 4)) {
		if (size == 4) {
			do_test_mult_sparse < StoreELL < float > > (f, iters);
		} else if (size == 8) {
			do_test_mult_sparse < StoreELL < double > > (f, iters);
		} else {
			throw runtime_error("unknown format\n");
		}
	} else {
		throw runtime_error("unknown format\n");
	}
}

int main (int argc, char * argv[])
{
	bool result        = true;
	bool mult_sparse   = false;
	bool invert_sparse = false;
	FILE * f = 0;

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

	linal_init();

	try
	{
		int has_double = check_device_supports_double();
		fprintf (stderr, "has double: %d\n", has_double);

		Timer t;
		//burn (1.0);
		t.restart();

		if (invert_sparse)
		{
			t.restart();
			fprintf (stderr, "invert_sparse (): %lf, %d\n", t.elapsed(), (int) result);
		}

		if (mult_sparse)
		{
			t.restart();
			result &= test_mult_sparse (f, iters);
			fprintf (stderr, "test_mult_sparse (): %lf, %d\n", t.elapsed(), (int) result);
		}
		fprintf (stderr, "elapsed: %lf\n", t.elapsed() );
	}
	catch (const std::exception & e)
	{
		fprintf (stderr, "exception: %s\n", e.what() );
	}

	linal_shutdown();
	fclose(f);
	return (int) (!result);
}

