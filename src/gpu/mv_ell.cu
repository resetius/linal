/* -*- charset: utf-8 -*- */
/*$Id$*/

/* Copyright (c) 2009-2010 Alexey Ozeritsky
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

#include "linal_cuda.h"
#include "texture.h"

register_texture(float, texX1);
register_texture(float, texAX);
register_texture(int, texAP);
register_texture(int, texAI);

namespace linal {

template < typename T, typename AIR, typename AXR, typename XR >
__global__ void ell_mult(
			   T * r, 
			   AIR Ai, 
			   AXR Ax,
			   XR x, 
			   int n,
			   int cols, 
			   int stride)
{
	int row = blockDim.x * blockIdx.x + threadIdx.x;
	if (row < n) {
		T sum = 0;

		for (int i0 = 0; i0 < cols; i0++){
			const T A_ij = Ax.get(stride * i0 + row);

			if (A_ij != 0) {
				const int col = Ai.get(stride * i0 + row);
				sum += A_ij * x.get(col);
			}
		}
	    r[row] = sum;
	}
}

__host__ void 
ell_mult_vector_r(float * r, const int * Ai, const float * Ax, 
	const float * x, int n, int cols, int stride)
{
	SPLAY2(n);

//	texture_reader(texX1) XR(x, n);
//	texture_reader(texAX) AXR(Ax, cols * stride);
//	texture_reader(texAI) AIR(Ai, cols * stride);

	simple_reader < float > XR(x);
	simple_reader < float > AXR(Ax);
	simple_reader < int > AIR(Ai);

	ell_mult<<<blocks, threads>>>(r, AIR, AXR, XR, n, cols, stride);
}

__host__ void 
ell_mult_vector_r(double * r, const int * Ai, const double * Ax, 
	const double * x, int n, int cols, int stride)
{
	SPLAY2(n);

	simple_reader < double > XR(x);
	simple_reader < double > AXR(Ax);
	simple_reader < int > AIR(Ai);
	
	ell_mult<<<blocks, threads>>>(r, AIR, AXR, XR, n, cols, stride);
}

}

