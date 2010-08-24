#ifndef SOLVER_H
#define SOLVER_H
/* -*- charset: utf-8 -*- */
/*$Id$*/

/* Copyright (c) 2009-2010 Alexey Ozeritsky (Алексей Озерицкий)
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

/**
 * @file
 * @author Alexey Ozeritsky <aozeritsky@gmail.com>
 * @version $Revision$
 *
 * @section DESCRIPTION
 * The Matrix class and linear system solver.
 */

#include <assert.h>

#include <vector>
#include <stdexcept>
#include <map>

#ifdef SPARSE
#ifndef GMRES
#define UMFPACK
#endif
#endif

#undef SUPERLU

#include "allocator.h"
#include "array.h"
#include "gmres.h"

namespace linal
{

/**
 * @defgroup solver Linear equations solver.
 * Contains Linear equations solver class and
 * auxiliary routines.
 * @{
 */

/**
 *  Compressed Sparsed Row Format (CSR)
 *  Also known as Compressed Row Storage (CRS)
 */
template < typename T, template < class > class Alloc = Allocator >
struct StoreCSR
{
	typedef T data_type;
	typedef StoreCSR < T, Alloc> my_type;

	int n_;
	int nz_;
	Array < int, Alloc < int > > Ap_; // количества ненулевых в столбцах
	Array < int, Alloc < int > > Ai_; // индексы (номер строки) ненулевых элементов
	Array < T, Alloc < T > > Ax_;   // ненулевые элементы матрицы

	StoreCSR() : n_ (0), nz_ (0) {}
	StoreCSR (int n, int nz) : n_ (n), nz_ (nz), Ap_ (n_ + 1), Ai_ (nz_), Ax_ (nz_) {}

	template < template < class > class A >
	my_type & operator = (const StoreCSR < T, A > & o)
	{
		n_  = o.n_;
		nz_ = o.nz_;
		array_copy(Ap_, o.Ap_);
		array_copy(Ai_, o.Ai_);
		array_copy(Ax_, o.Ax_);
	}

	void resize (int n, int nz)
	{
		n_  = n;
		nz_ = nz;
		Ap_.resize (n_ + 1);
		Ai_.resize (nz_);
		Ax_.resize (nz_);
	}

	bool empty()
	{
		return n_ == 0;
	}

	typedef std::map < int , T > row_t;
	typedef std::vector < row_t > sparse_t;

	/**
	 * Fill CSR Matrix from unstructured data
	 */
	void load (const sparse_t & unstruct);

	void mult (T * r, const T * x) const;

	/**
	 * Прибавляет матрицу A размера mxm.
	 * Каждая строка матрицы A умножается поэлементно на vec.
	 * Структура строк матрицы A ДОЛЖНА совпадать со структурой
	 * главного минора размера mxm матрицы this.
	 */
	void add_matrix1 (const my_type & A, const T * vec);

	/**
	 * Прибавляет матрицу A размера mxm.
	 * Каждый столбец матрицы A умножается поэлементно на vec.
	 * Структура строк матрицы A ДОЛЖНА совпадать со структурой
	 * главного минора размера mxm матрицы this.
	 */
	void add_matrix2 (const my_type & A, const T * vec);

	void print(FILE * f) const;
	void dump(FILE * f) const;
	void restore(FILE * f);
};

template < typename T, template < class > class Alloc = Allocator >
struct StoreELL
{
	typedef T data_type;
	typedef StoreELL < T, Alloc > my_type;

	int n_;
	int nz_;
	int cols_;
	int stride_;

	Array < int, Alloc < int > > Ai_;
	Array < T, Alloc < T > > Ax_;

	StoreELL() : n_ (0), nz_ (0), cols_ (0), stride_ (0) {}

	StoreELL (int n, int nz, int cols) :
			n_ (n), nz_ (nz), cols_ (cols), stride_ (32 * ( (cols_ + 32 - 1) / 32) ),
			Ai_ (cols_ * stride_), Ax_ (cols_ * stride_)
	{
	}

	void resize (int n, int nz, int cols)
	{
		n_    = n;
		nz_   = nz;
		cols_ = cols;
		Ai_.resize (cols_ * stride_);
		Ax_.resize (cols_ * stride_);
	}

	bool empty() const
	{
		return n_ == 0;
	}

	typedef std::map < int , T > row_t;
	typedef std::vector < row_t > sparse_t;

	/**
	 * Fill ELL Matrix from unstructured data
	 */
	void load (const sparse_t & unstruct);

	void mult (T * r, const T * x) const;

	void add_matrix1 (const my_type & A, const T * vec)
	{
		assert (0);
	}

	void add_matrix2 (const my_type & A, const T * vec)
	{
		assert (0);
	}

	void print(FILE * f) const
	{
		// implement
	}

	void dump(FILE * f) const;
	void restore(FILE * f);
};

template < typename Store1, typename Store2 >
struct DoubleStore
{
	typedef Store1 mult_t;
	typedef Store2 invert_t;
	typedef DoubleStore < Store1, Store2 > my_type;

	Store1 mult;
	Store2 invert;

	void add_matrix1 (const my_type & A, const typename Store1::data_type * vec)
	{
		mult.add_matrix1 (A.mult, vec);
		invert.add_matrix1 (A.invert, vec);
	}

	void add_matrix2 (const my_type & A, const typename Store1::data_type * vec)
	{
		mult.add_matrix2 (A.mult, vec);
		invert.add_matrix2 (A.invert, vec);
	}

	void print(FILE * f)
	{
		mult.print(f);
	}

	void dump(FILE * f)
	{
		mult.dump(f);
	}

	void restore(FILE * f)
	{
		throw std::runtime_error("restore not supported\n");
	}
};

template < typename Store >
struct DoubleStore < Store, Store >
{
	typedef Store mult_t;
	typedef Store invert_t;
	typedef DoubleStore < Store, Store > my_type;

	Store both;
	Store & mult;
	Store & invert;
	DoubleStore () : mult (both), invert (both) {}
	DoubleStore (const DoubleStore < Store, Store > & r)
		: both(r.both), mult (both), invert (both)
	{
	}

	void add_matrix1 (const my_type & A, const typename Store::data_type * vec)
	{
		mult.add_matrix1 (A.mult, vec);
	}

	void add_matrix2 (const my_type & A, const typename Store::data_type * vec)
	{
		mult.add_matrix2 (A.mult, vec);
	}

	void print(FILE * f)
	{
		mult.print(f);
	}

	void dump(FILE * f)
	{
		mult.dump(f);
	}

	void restore(FILE * f)
	{
		mult.restore(f);
	}
};

/**
 * Solver class.
 * В солвере может быть два контейнера с разными аллокаторами для обращения
 * и для умножения.
 */
template < typename T, typename MultStore, typename InvStore = MultStore >
class SparseSolver
{
protected:
	typedef DoubleStore < MultStore, InvStore > store_t;
	mutable store_t store_;
	typedef std::map < int , T > row_t;
	typedef std::vector < row_t > sparse_t;

	mutable sparse_t A_;

public:
	typedef T data_type;
	typedef SparseSolver < T, MultStore, InvStore > my_type;

	SparseSolver (int n) : A_ (n)
	{
	}

	template < typename Store1, typename Store2 >
	SparseSolver (const Store1 & s1, const Store2 & s2)
	{
		if (s1.n_ != s2.n_) {
			throw std::runtime_error("uncompatible stores\n");
		}

		store_.mult   = s1;
		store_.invert = s2;
	}

	int dim() const
	{
		return store_.mult.n_;
	}

	int nonzero() const
	{
		return store_.mult.nz_;
	}

	/**
	 *  Add a number to element (i, j) (A[i][j] += a).
	 *  @param i - index
	 *  @param j - index
	 *  @param a - value
	 */
	void add (int i, int j, T a);

	/**
	 * Solve equation Ax = b.
	 * That function uses GMRES
	 * @param x - answer
	 * @param b - right part
	 */
	void solve (T * x, const T * b) const;

	/**
	 * Product of matrix by vector (out = A in).
	 * @param out - result
	 * @param in  - input vector
	 */
	void mult_vector (T * out, const T * in) const;

	/**
	 * print matrix to file.
	 */
	void print(FILE * f) const;
	void dump(FILE * f) const;
	void restore(FILE * f);

	/**
	 * fill permanent storage from temporary storage
	 */
	void prepare() const;

	void add_matrix1 (my_type & A, const T * vec);
	void add_matrix2 (my_type & A, const T * vec);
};

#ifdef UMFPACK
// implements UmfPackSolver
#include "impl/solver_umfpack.h"
#endif

#ifdef SUPERLU
// implements SuperLUSolver
#include "impl/solver_superlu.h"
#endif

/**
 * Matrix class.
 */
template < typename T >
class SimpleSolver
{
	int n_;
	Array < T, Allocator < T > > A_;  // матрица

public:
	typedef T data_type;

	SimpleSolver (int n) : n_ (n), A_ (n * n) {}

	/**
	 *  Add a number to element (i, j) (A[i][j] += a).
	 *  @param i - index
	 *  @param j - index
	 *  @param a - value
	 */
	void add (int i, int j, T a);

	/**
	 * Solve equation Ax = b.
	 * That function uses Gauss
	 * @param x - answer
	 * @param b - right part
	 */
	void solve (T * x, const T * b);

	/**
	 * Product of matrix by vector (out = A in).
	 * @param out - result
	 * @param in  - input vector
	 */
	void mult_vector (T * out, const T * in);

	/**
	 * print matrix to file.
	 */
	void print(FILE * f);
	void dump(FILE * f);
	void restore(FILE * f);

	/**
	 * Do nothing. For compatibility with sparse solvers.
	 */
	void prepare() {}
};


#if defined(SUPERLU) && !defined(GPGPU)

template < typename T >
class Solver: public SuperLUSolver < T, StoreCSR < T , Allocator >  >
{
	typedef SuperLUSolver < T, StoreCSR < T , Allocator >  > base;
public:
	Solver (int n) : base (n) {}
};

#elif defined(UMFPACK) && !defined(GPGPU)

template < typename T >
class Solver: public UmfPackSolver < T, StoreCSR < T , Allocator >  >
{
	typedef UmfPackSolver < T, StoreCSR < T , Allocator >  > base;
public:
	Solver (int n) : base (n) {}
};

#else
template < typename T >
class Solver: public SparseSolver < T, StoreELL < T , Allocator > , StoreELL < T , Allocator > >
{
	typedef SparseSolver < T, StoreELL < T , Allocator > , StoreELL < T , Allocator > >  base;
public:
	Solver (int n) : base (n) {}
};
#endif

/** @} */ /* solver */

#include "impl/solver_impl.h"

}

#endif /* SOLVER_H */

