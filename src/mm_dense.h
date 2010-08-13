#ifndef PHELM_LA_H
#define PHELM_LA_H
/* -*- charset: utf-8 -*- */
/*$Id$*/
/**
 * @file
 * @author Alexey Ozeritsky <aozeritsky@gmail.com>
 *
 * Copyright (c) 2009-2010 Alexey Ozeritsky
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

namespace linal
{

/**
 * @defgroup la Linear Algebra functions and Classes.
 * @{
 */

/**
 * Product of NxN matrix and vector
 * @param r - output vector, r = Ax
 * @param A - intput matrix
 * @param x - input vector
 * @param n - dimension of matrix and vector
 */
void mat_mult_vector (double * r, const double * A, const double * x, int n);

/**
 * Product of NxN matrix and vector
 * @param r - output vector, r = Ax
 * @param A - intput matrix
 * @param x - input vector
 * @param n - dimension of matrix and vector
 */
void mat_mult_vector (float * r, const float * A, const float * x, int n);

/**
 * Product of NxN matrix and vector (stupid algorithm)
 * @param r - output vector, r = Ax
 * @param A - intput matrix
 * @param x - input vector
 * @param n - dimension of matrix and vector
 */
void mat_mult_vector_stupid (double * r, const double * A, const double * x, int n);

/**
 * Product of NxN matrix and vector (stupid algorithm)
 * @param r - output vector, r = Ax
 * @param A - intput matrix
 * @param x - input vector
 * @param n - dimension of matrix and vector
 */
void mat_mult_vector_stupid (float * r, const float * A, const float * x, int n);

/**
 * @}
 */

} /* namespace */

#endif /* PHELM_LA_H */

