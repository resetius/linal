#ifndef MM_DENSE_H
#define MM_DENSE_H

namespace linal {
void
mat_mult_mat(double * C, const double * A, const double * B, int n);
void
mat_mult_mat(float * C, const float * A, const float * B, int n);
}

#endif /* MM_DENSE_H */

