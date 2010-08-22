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

template < typename T, typename APR, typename AIR, typename XR, typename AXR >
__global__ void sparse_mult_vector_csr_(T * r, 
	APR Ap, 
	AIR Ai, 
	AXR Ax,
	XR x, 
	int n)
{
	int threads = gridDim.x  * blockDim.x;
	int i, i0, to, j;
	int start = blockDim.x * blockIdx.x + threadIdx.x;

	for (j = start; j < n; j += threads) {
		i0 = Ap.get(j);    
		to = Ap.get(j + 1);

		T rj = (T)0.0;

		for (; i0 < to; ++i0) {
			i   = Ai.get(i0);
			rj += Ax.get(i0) * x.get(i);
		}

		r[j] = rj;
	}
}

__host__ void csr_mult_vector_r(double * r, 
	const int * Ap, 
	const int * Ai, 
	const double * Ax,
	const double * x, 
	int n,
	int nz)
{
	SPLAY(n);
	simple_reader < double > XR(x);
	simple_reader < double > AXR(Ax);
	simple_reader < int > AIR(Ap);
	simple_reader < int > APR(Ap);

	sparse_mult_vector_csr_ <<< blocks, threads >>> (r, APR, AIR, AXR, XR, n);
}

__host__ void csr_mult_vector_r(float * r, 
	const int * Ap, 
	const int * Ai, 
	const float * Ax,
	const float * x, 
	int n,
	int nz)
{
	SPLAY(n);

	bool useTexture;

	useTexture = ((n + 1 < MAX_1DBUF_SIZE) && (nz < MAX_1DBUF_SIZE));

	if (n < 1000 || n > 10000) /* experimental bound */
	{
		useTexture = false;
	}

	if (useTexture) {
		texture_reader(texX1) XR(x, n);
		texture_reader(texAX) AXR(Ax, nz);
		texture_reader(texAI) AIR(Ai, nz);
		texture_reader(texAP) APR(Ap, n + 1);

		sparse_mult_vector_csr_ <<< blocks, threads >>> (r, APR, AIR, AXR, XR, n);
	} else {
		simple_reader < float > XR(x);
		simple_reader < float > AXR(Ax);
		simple_reader < int > AIR(Ai);
		simple_reader < int > APR(Ap);

		sparse_mult_vector_csr_ <<< blocks, threads >>> (r, APR, AIR, AXR, XR, n);
	}
}

}

