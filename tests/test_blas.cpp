#include <vector>
#include <math.h>
#include "linal.h"

using namespace std;
using namespace linal;

void init_a_b(double * a, double * b, int n)
{
	for (int i = 0; i < n; ++i) {
		a[i] = 1;
		b[i] = 2;
	}
}

int main()
{
	int n = 10000;
	vector < double > a(n), b(n), c(n);

	init_a_b(&a[0], &b[0], n);
	vec_sum(&c[0], &a[0], &b[0], n);

	for (int i = 0; i < n; ++i) {
		if (fabs(c[i] - 3) > 1e-15) {
			fprintf(stderr, "vec_sum failed\n");
			return -1;
		}
	}

	init_a_b(&a[0], &b[0], n);
	vec_sum(&a[0], &a[0], &b[0], n);

	for (int i = 0; i < n; ++i) {
		if (fabs(a[i] - 3) > 1e-15) {
			fprintf(stderr, "vec_sum failed\n");
			return -1;
		}
	}

	init_a_b(&a[0], &b[0], n);
	vec_sum2(&c[0], &a[0], &b[0], -5, n);

	for (int i = 0; i < n; ++i) {
		if (fabs(c[i] - (-9)) > 1e-15) {
			fprintf(stderr, "vec_sum2 failed\n");
			return -1;
		}
	}

	init_a_b(&a[0], &b[0], n);
	vec_diff(&c[0], &a[0], &b[0], n);

	for (int i = 0; i < n; ++i) {
		if (fabs(c[i] - (-1)) > 1e-15) {
			fprintf(stderr, "vec_diff failed\n");
			return -1;
		}
	}

	init_a_b(&a[0], &b[0], n);
	vec_diff(&a[0], &a[0], &b[0], n);

	for (int i = 0; i < n; ++i) {
		if (fabs(a[i] - (-1)) > 1e-15) {
			fprintf(stderr, "vec_diff failed\n");
			return -1;
		}
	}

	init_a_b(&a[0], &b[0], n);
	vec_mult(&c[0], &a[0], &b[0], n);

	for (int i = 0; i < n; ++i) {
		if (fabs(c[i] - (2)) > 1e-15) {
			fprintf(stderr, "vec_mult failed\n");
			return -1;
		}
	}

	init_a_b(&a[0], &b[0], n);
	vec_sum1(&c[0], &a[0], &b[0], 2.0, -2.0, n);

	for (int i = 0; i < n; ++i) {
		if (fabs(c[i] - (-2)) > 1e-15) {
			fprintf(stderr, "vec_sum1 failed\n");
			return -1;
		}
	}

	init_a_b(&a[0], &b[0], n);
	vec_sum1(&a[0], &a[0], &b[0], 2.0, -2.0, n);

	for (int i = 0; i < n; ++i) {
		if (fabs(a[i] - (-2)) > 1e-15) {
			fprintf(stderr, "vec_sum1 failed\n");
			return -1;
		}
	}

	init_a_b(&a[0], &b[0], n);
	vec_mult_scalar(&c[0], &a[0], 2.0, n);

	for (int i = 0; i < n; ++i) {
		if (fabs(c[i] - (2)) > 1e-15) {
			fprintf(stderr, "vec_mult_scalar failed\n");
			return -1;
		}
	}

	init_a_b(&a[0], &b[0], n);
	vec_mult_scalar(&a[0], &a[0], 2.0, n);

	for (int i = 0; i < n; ++i) {
		if (fabs(a[i] - (2)) > 1e-15) {
			fprintf(stderr, "vec_mult_scalar failed\n");
			return -1;
		}
	}

	return 0;
}

