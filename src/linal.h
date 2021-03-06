#ifndef PHELM_LA_H
#define PHELM_LA_H
/* -*- charset: utf-8 -*- */
/*$Id$*/
/**
 * @file
 * @author Alexey Ozeritsky <aozeritsky@gmail.com>
 *
 * @page License
 * @section LICENSE
 *
 * @verbatim
  Copyright (c) 2009-2010 Alexey Ozeritsky
  All rights reserved.
 
  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions
  are met:
  1. Redistributions of source code must retain the above copyright
     notice, this list of conditions and the following disclaimer.
  2. Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the following disclaimer in the
     documentation and/or other materials provided with the distribution.
  3. The name of the author may not be used to endorse or promote products
     derived from this software without specific prior written permission
 
  THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
  IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
  IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
  NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
  THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * @endverbatim
 *
 * @mainpage Linal Documentation
 * @section into_sec Introduction
 *
 * This library contains linear algebra subprograms such as
 * - vector operations 
 *   - vectors sum
 *   - dot products 
 *   - vector norm
 * - matrix vector operations
 *   - matrix vector product
 * - linear systems solver
 * - input/output for vectors and matrices
 *
 * The library supports the following matrices types
 * - dense
 * - sparse
 *   - CSR
 *   - ELL 
 * 
 * The linear system solver supports the following backends
 * - Gauss (for dense matrices)
 * - SuperLU (http://crd.lbl.gov/~xiaoye/SuperLU/)
 * - UMFPACK (http://www.cise.ufl.edu/research/sparse/umfpack/)
 * - GMRES (included)
 *
 * The vector operations subprograms have the following backends
 * - Simple
 * - OpenMP
 * - CUDA
 *
 * The library supports the following data types
 * - float
 * - double
 *
 */

#include <stdio.h>
#include <assert.h>
#include <string.h>

#ifdef GPGPU
#include <cuda_runtime_api.h>
#endif

#include "allocator.h"

namespace linal
{

/**
 * @defgroup la Linear Algebra functions and Classes.
 * @{
 */

/**
 * The Gauss method. Solves a linear system Ax=b.
 *
 * @param A - the matrix of the system.
 * @param b - the right part.
 * @param x - the answer.
 * @param n - dimension.
 * @return 0 on success.
 */
int gauss (double *A, double *b, double *x, int n);

/**
 * Print NxM matrix to file.
 * @param A - matrix
 * @param n - dimension
 * @param m - dimension
 * @param fmt - format
 */
void mat_print (FILE * f, const double * A, int n, int m, const char * fmt);
void mat_print (const char * f, const double * A, int n, int m, const char * fmt);

/**
 * Print NxM matrix to file.
 * @param A - matrix
 * @param n - dimension
 * @param m - dimension
 * @param fmt - format
 */
void mat_print (FILE * f, const float * A, int n, int m, const char * fmt);
void mat_print (const char * f, const float * A, int n, int m, const char * fmt);

/**
 * Load NxM matrix from file.
 */
void mat_load (FILE * f, std::vector < double > & A, int * n, int * m);
void mat_load (const char * f, std::vector < double > & A, int * n, int * m);
void mat_load (FILE * f, std::vector < float > & A, int * n, int * m);
void mat_load (const char * f, std::vector < float > & A, int * n, int * m);

/**
 * Print vector to file.
 * @param A - vector
 * @param n - dimension
 * @param fmt - format
 */
void vec_print (FILE * f, const double * A, int n, const char * fmt);
void vec_print (const char * f, const double * A, int n, const char * fmt);

/**
 * Print vector to file.
 * @param A - vector
 * @param n - dimension
 * @param fmt - format
 */
void vec_print (FILE * f, const float * A, int n, const char * fmt);
void vec_print (const char * f, const float * A, int n, const char * fmt);

/**
 * Finds minimal element of vector.
 */
double vec_find_min (const double * v, int n);

/**
 * Finds maximal element of vector.
 */
double vec_find_max (const double * v, int n);

/**
 * nxm matrix transpose.
 */
void mat_transpose (double * out, const double * in, int n, int m);
/**
 * nxm matrix transpose and multiply by k.
 */
void mat_transpose1 (double * out, const double * in, double k, int n, int m);

/**
 * @param n - число строк в матрице Ax
 */
void csr_add_matrix1 (const int * oAp, double * oAx,
                      const int * Ap, const int * Ai, const double * Ax,
                      const double * x, int n);

void csr_add_matrix1 (const int * oAp, float * oAx,
                      const int * Ap, const int * Ai, const float * Ax,
                      const float * x, int n);

void csr_add_matrix2 (const int * oAp, double * oAx,
                      const int * Ap, const int * Ai, const double * Ax,
                      const double * x, int n);

void csr_add_matrix2 (const int * oAp, float * oAx,
                      const int * Ap, const int * Ai, const float * Ax,
                      const float * x, int n);

/**
 * Print CSR sparse matrix to file.
 * @param A - input sparse matrix
 * @param n - dimension of sparse matrix
 * @param f - output file
 */
void csr_print (const int * Ap, const int * Ai,
                   const double * Ax, int n, FILE * f);

void csr_print (const int * Ap, const int * Ai,
                   const float * Ax, int n, FILE * f);

/**
 * Print ELL sparse matrix to file.
 * @param A - input sparse matrix
 * @param n - dimension of sparse matrix
 * @param f - output file
 */
void ell_print (const int * Ai, const double * Ax, int rows, int cols, int stride, FILE * f);

void ell_print (const int * Ai, const float * Ax, int rows, int cols, int stride, FILE * f);

/**
 * Linear combination of two vectors.
 * @param r - output vector
 * @param a - input vector
 * @param b - input vector
 * @param k1 - coefficient
 * @param k2 - coefficient
 * @param n - dimension of vectors
 * @return r = k1 * a + k2 * b
 */
void vec_sum1 (double * r, const double * a,
               const double *b, double k1, double k2, int n);

void vec_sum1 (float * r, const float * a,
               const float *b, float k1, float k2, int n);

/**
 * Linear combination of two vectors.
 * @param r - output vector
 * @param a - input vector
 * @param b - input vector
 * @param k2 - coefficient
 * @param n - dimension of vectors
 * @return r = a + k2 * b
 */
void vec_sum2 (double * r, const double * a,
               const double *b, double k2, int n);

void vec_sum2 (float * r, const float * a,
               const float *b, float k2, int n);

/**
 * Product of vector by number.
 * @param a - output vector
 * @param b - input vector
 * @param k - input number
 * @param n - dimension of vector
 * @return a = b * k
 */
void vec_mult_scalar (double * a, const double * b, double k, int n);
void vec_mult_scalar (float * a, const float * b, float k, int n);

/**
 * Sum of two vectors.
 * @param r - output vector
 * @param a - input vector
 * @param b - input vector
 * @param n - dimension of vectors
 * @return r = a + b
 */
void vec_sum (double * r, const double * a, const double *b, int n);
void vec_sum (float * r, const float * a, const float *b, int n);

/**
 * Vector norm.
  \f[
  \sqrt{\sum_{i=0}v_i^2}
  \f]
 * @param v - input vector
 * @param n - dimension of vector
 * @return vector norm
 */
double vec_norm2 (const double *v, int n);
float vec_norm2 (const float *v, int n);

/**
 * Inner product of two vectors.
   \f[
   \sum_{i=0}^{n}a_i b_i
   \f]
 * @param a - input vector
 * @param b - input vector
 * @param n - dimension of vectors
 * @return inner product of a and b
 */
double vec_scalar2 (const double * a, const double * b, int n);
float vec_scalar2 (const float * a, const float * b, int n);

/**
 * Element by element vector multiplication.
 * r = a * b
 * @param r - output vector
 * @param a - input vector
 * @param b - input vector
 * @param n - dimension of vectors
 */
void vec_mult (double * r, const double * a, const double * b, int n);
void vec_mult (float * r, const float * a, const float * b, int n);

/**
 * Difference of two vectors.
 * @param r - the output vector
 * @param a - the input vector
 * @param b - the input vector
 * @param n - the dimension of vectors
 * @return r = a - b
 */
void vec_diff (double * r, const double * a, const double * b, int n);
void vec_diff (float * r, const float * a, const float * b, int n);

/**
 * Copy vector a to vector b.
 * @param b - the output vector
 * @param a - the input vector
 * @param n - the dimension of vectors
 */
template < typename T >
void vec_copy (T * b, const T * a, int n)
{
#ifndef GPGPU
	memcpy (b, a, n * sizeof (T) );
#else
	cudaMemcpy (b, a, n * sizeof (T), cudaMemcpyDeviceToDevice);
#endif
}

template < typename T >
void vec_copy_from_host (T * b, const T * a, int n)
{
#ifndef GPGPU
	memcpy (b, a, n * sizeof (T) );
#else
	cudaMemcpy (b, a, n * sizeof (T), cudaMemcpyHostToDevice);
#endif
}

template < typename T >
void vec_copy_from_device (T * b, const T * a, int n)
{
#ifndef GPGPU
	memcpy (b, a, n * sizeof (T) );
#else
	cudaMemcpy (b, a, n * sizeof (T), cudaMemcpyDeviceToHost);
#endif
}

int check_device_supports_double();
void linal_init();
void linal_shutdown();
void linal_sync();

void set_num_threads (int threads);

/**
 * @}
 */

} /* namespace */

#include "mm_dense.h"
#include "mv_dense.h"
#include "mv_csr.h"
#include "mv_ell.h"
#include "pow.h"
#include "timer.h"
#include "fpe.h"
#include "array.h"
#include "array_pool.h"

#endif /* PHELM_LA_H */

