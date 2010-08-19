/* -*- charset: utf-8 -*- */
/*$Id$*/

/* Copyright (c) 2009 Alexey Ozeritsky (Алексей Озерицкий)
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

template < typename T, template < class > class Alloc >
void StoreCSR < T, Alloc > ::load (const std::vector < row_t > & A)
{
	int idx = 0;
	n_  = (int) A.size();
	nz_ = 0; // non-null elements
	for (uint i = 0; i < A.size(); ++i)
	{
		nz_ += (int) A[i].size();
	}

	Ax_.resize (nz_);
	Ai_.resize (nz_);
	Ap_.resize (n_ + 1);

	std::vector < T >   Ax (nz_);
	std::vector < int > Ai (nz_);
	std::vector < int > Ap (n_ + 1);

	Ap[0] = 0;

	for (uint i = 0; i < A.size(); ++i)
	{
		Ap[i + 1] = Ap[i] + (int) A[i].size();
		for (typename row_t::const_iterator it = A[i].begin();
		        it != A[i].end(); ++it)
		{
			Ax[idx] = it->second;
			Ai[idx] = it->first;
			idx += 1;
		}
	}

	vec_copy_from_host (&Ax_[0], &Ax[0], nz_);
	vec_copy_from_host (&Ai_[0], &Ai[0], nz_);
	vec_copy_from_host (&Ap_[0], &Ap[0], n_ + 1);
}

template < typename T, template < class > class Alloc >
void StoreCSR < T, Alloc > ::mult (T * r, const T * x) const
{
	csr_mult_vector_r (r, &Ap_[0], &Ai_[0], &Ax_[0], x, n_, nz_);
}

template < typename T, template < class > class Alloc >
void StoreCSR < T, Alloc > ::add_matrix1 (const my_type & A, const T * x)
{
	csr_add_matrix1 (&Ap_[0], &Ax_[0],
	                 &A.Ap_[0], &A.Ai_[0], &A.Ax_[0],
	                 x, A.n_);
}

template < typename T, template < class > class Alloc >
void StoreCSR < T, Alloc > ::add_matrix2 (const my_type & A, const T * x)
{
	csr_add_matrix2 (&Ap_[0], &Ax_[0],
	                 &A.Ap_[0], &A.Ai_[0], &A.Ax_[0],
	                 x, A.n_);
}

template < typename T, template < class > class Alloc >
void StoreCSR < T, Alloc > ::print(FILE * f)
{
	sparse_print(&Ap_[0], &Ai_[0], &Ax_[0], n_, f);
}

template < typename T, template < class > class Alloc >
void StoreCSR < T, Alloc > ::dump(FILE * f)
{
	int size = sizeof(T);
	fwrite("CSR ", 4, 1, f);
	fwrite(&size, 4, 1, f);
	fwrite(&n_,  sizeof(n_), 1, f);
	fwrite(&nz_, sizeof(nz_), 1, f);
	if (nz_) {
		fwrite(&Ap_[0], sizeof(int), Ap_.size(), f);
		fwrite(&Ai_[0], sizeof(int), nz_, f);
		fwrite(&Ax_[0], sizeof(T), nz_, f);
	}
}

template < typename T, template < class > class Alloc >
void StoreCSR < T, Alloc > ::restore(FILE * f)
{
	char tag[4];
	int size;
	off_t fsize;

	fseek(f, 0L, SEEK_END);
	fsize = ftell(f);
	fseek(f, 0L, SEEK_SET);

	if (fread(tag, 4, 1, f) != 1 || memcmp(tag, "CSR ", 4) != 0) {
		throw std::runtime_error("bad file format");
	}

	if (fread(&size, 4, 1, f) != 1 || size != sizeof(T)) {
		throw std::runtime_error("bad file format");
	}

	if (fread(&n_, sizeof(n_), 1, f) !=1 ||	fread(&nz_, sizeof(nz_), 1, f) != 1)
	{
		throw std::runtime_error("bad file format");
	}

	if (fsize != (n_ + 1) * sizeof(int) + nz_ * sizeof(int) + nz_ * sizeof(T) + 4 + 4 + 4 + 4)
	{
		throw std::runtime_error("bad file format");
	}

	if (nz_) {
		Ap_.resize(n_ + 1);
		Ai_.resize(nz_);
		Ax_.resize(nz_);
		fread(&Ap_[0], sizeof(int), n_ + 1, f);
		fread(&Ai_[0], sizeof(int), nz_, f);
		fread(&Ax_[0], sizeof(T), nz_, f);
	}
}

template < typename T, template < class > class Alloc >
void StoreELL < T, Alloc > ::load (const std::vector < row_t > & A)
{
	nz_ = 0; // non-null elements
	n_  = (int) A.size();
	for (uint i = 0; i < A.size(); ++i)
	{
		nz_ += (int) A[i].size();
	}
	cols_ = 0;
	for (uint i = 0; i < A.size(); ++i)
	{
		cols_ = std::max (cols_, (int) A[i].size() );
	}

	int stride_off = 1; //32;
	stride_   = stride_off * ( (n_ + stride_off - 1) / stride_off);
	int count = cols_ * stride_;

//	stride_   = stride_off * ((cols_ + stride_off - 1) / stride_off);
//	int count = n_ * stride_;

	Ax_.resize (count);
	Ai_.resize (count);
//	for (int i = 0; i < count; ++i) {
//		Ai_[i] = -1;
//	}

	std::vector < T > Ax (count);
	std::vector < int > Ai (count);

	for (uint i = 0; i < A.size(); ++i)
	{
		int idx = 0;
		for (typename row_t::const_iterator it = A[i].begin();
		        it != A[i].end(); ++it)
		{
			Ax[idx * stride_ + i] = it->second;
			Ai[idx * stride_ + i] = it->first;
			//Ax[i * stride_ + idx] = it->second;
			//Ai[i * stride_ + idx] = it->first;
			idx++;
		}
	}

	vec_copy_from_host (&Ax_[0], &Ax[0], count);
	vec_copy_from_host (&Ai_[0], &Ai[0], count);
}

template < typename T, template < class > class Alloc >
void StoreELL < T, Alloc > ::mult (T * r, const T * x) const
{
	ell_mult_vector_r (r, &Ai_[0], &Ax_[0], x, n_, cols_, stride_);
}

template < typename T, template < class > class Alloc >
void StoreELL < T, Alloc > ::dump (FILE * f)
{
	int size = sizeof(T);
	fwrite("ELL ", 4, 1, f);
	fwrite(&size, 4, 1, f);

	fwrite(&n_, sizeof(n_), 1, f);
	fwrite(&nz_, sizeof(nz_), 1, f);
	fwrite(&cols_, sizeof(cols_), 1, f);
	fwrite(&stride_, sizeof(stride_), 1, f);
	if (nz_) {
		fwrite(&Ai_[0], sizeof(int), cols_ * stride_, f);
		fwrite(&Ax_[0], sizeof(T), cols_ * stride_, f);
	}
}

template < typename T, template < class > class Alloc >
void StoreELL < T, Alloc > ::restore (FILE * f)
{
	char tag[4];
	int size;
	off_t fsize;

	fseek(f, 0L, SEEK_END);
	fsize = ftell(f);
	fseek(f, 0L, SEEK_SET);

	if (fread(tag, 4, 1, f) != 1 || memcmp(tag, "CSR ", 4) != 0) {
		throw std::runtime_error("bad file format");
	}

	if (fread(&size, 4, 1, f) != 1 || size != sizeof(T)) {
		throw std::runtime_error("bad file format");
	}

	if (fsize != 4 + 4 + sizeof(n_) + sizeof(nz_) + sizeof(cols_) + sizeof(stride_)
			+ cols_ * stride_ * sizeof(int) + cols_ * stride_ * sizeof(T))
	{
		throw std::runtime_error("bad file format");
	}

	fread(&n_, sizeof(n_), 1, f);
	fread(&nz_, sizeof(nz_), 1, f);
	fread(&cols_, sizeof(cols_), 1, f);
	fread(&stride_, sizeof(stride_), 1, f);
	if (nz_) {
		Ai_.resize(cols_ * stride_);
		Ax_.resize(cols_ * stride_);
		fread(&Ai_[0], sizeof(int), cols_ * stride_, f);
		fread(&Ax_[0], sizeof(T), cols_ * stride_, f);
	}
}
	
template < typename T >
void SimpleSolver < T > ::mult_vector (T * out, const T * in)
{
	mat_mult_vector (out, &A_[0], in, n_);
}

template < typename T >
void SimpleSolver < T > ::solve (T * x, const T * b)
{
	gauss (&A_[0], &b[0], &x[0], n_);
}

template < typename T >
void SimpleSolver < T > ::add (int i, int j, T a)
{
	A_[i * n_ + j] += a;
}

template < typename T >
void SimpleSolver < T > ::print(FILE * f)
{
	mat_print (f, &A_[0], n_, n_, "%23.16le ");
}

template < typename T, typename MultStore, typename InvStore  >
void SparseSolver < T, MultStore, InvStore > ::prepare()
{
	if (store_.mult.empty() )
	{
		store_.mult.load (A_);
	}

	if (store_.invert.empty() )
	{
		store_.invert.load (A_);
	}
}

template < typename T, typename MultStore, typename InvStore  >
void SparseSolver < T, MultStore, InvStore > ::add_matrix1 (my_type & A, const T * vec)
{
	prepare();
	A.prepare();
	store_.add_matrix1 (A.store_, vec);
}

template < typename T, typename MultStore, typename InvStore  >
void SparseSolver < T, MultStore, InvStore > ::add_matrix2 (my_type & A, const T * vec)
{
	prepare();
	A.prepare();
	store_.add_matrix2 (A.store_, vec);
}

template < typename T, typename MultStore, typename InvStore  >
void SparseSolver < T, MultStore, InvStore > ::mult_vector (T * out, const T * in)
{
	prepare();
	store_.mult.mult (out, in);
}

template < typename T, typename MultStore, typename InvStore  >
void SparseSolver < T, MultStore, InvStore > ::print(FILE * f)
{
	prepare();
	store_.print(f);
}

template < typename T, typename MultStore, typename InvStore  >
void SparseSolver < T, MultStore, InvStore > ::dump(FILE * f)
{
	prepare();
	store_.dump(f);
}

template < typename T, typename MultStore, typename InvStore  >
void SparseSolver < T, MultStore, InvStore > ::restore(FILE * f)
{
	prepare();
	store_.restore(f);
}

template < typename T, typename MultStore, typename InvStore  >
void SparseSolver < T, MultStore, InvStore > ::add (int i, int j, T a)
{
	A_[i][j] += a;
}

template < typename T, typename Invert >
void SparseSolver__Ax (T * r, const Invert * invert, const T * x, int n)
{
	invert->mult (r, x);
}

template < typename T, typename MultStore, typename InvStore  >
void SparseSolver < T, MultStore, InvStore > ::solve (T * x, const T * b)
{
	prepare();
	gmres (&x[0], &store_.invert, &b[0], SparseSolver__Ax < T, typename store_t::invert_t > ,
	       store_.invert.n_, 100, 1000);
}

template < typename MultStore, typename InvStore >
SparseSolver < typename MultStore::data_type, typename InvStore::data_type >
make_sparse_solver(MultStore & store1, InvStore & store2)
{
	SparseSolver < typename MultStore::data_type, typename InvStore::data_type > ret(0);
	ret.store_.mult   = store1;
	ret.store_.invert = store2; 
}

