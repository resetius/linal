#ifndef PHELM_ARRAY_H
#define PHELM_ARRAY_H
/* -*- charset: utf-8 -*- */
/*$Id$*/
/**
 * Copyright (c) 2009-2010 Alexey Ozeritsky
 * All rights reserved.
 
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

template < typename T, typename Alloc >
class Array
{
	typedef Array < T, Alloc > my_type;

	size_t size_;
	T * data_;
	Alloc alloc_;

public:
	Array() : size_ (0), data_ (0) {}
	Array (size_t size) : size_ (size), data_ (0)
	{
		data_ = alloc_.allocate (size_);
	}

	~Array()
	{
		if (data_)
		{
			alloc_.deallocate (data_, size_);
		}
	}

	Array (const my_type & other) : size_ (0), data_ (0)
	{
		operator = (other);
	}

	Array & operator = (const my_type & other)
	{
		if (data_)
		{
			alloc_.deallocate (data_, size_);
			data_ = 0;
		}
		size_ = other.size_;
		if (size_ > 0)
		{
			data_ = alloc_.allocate (size_);
			vec_copy (data_, other.data_, (int) size_);
		}
		return *this;
	}

	void resize (size_t size)
	{
		if (size > size_)
		{
			T * p = alloc_.allocate (size);
			if (data_)
			{
				vec_copy (p, data_, (int) size);
				alloc_.deallocate (data_, size_);
			}
			data_ = p;
			size_ = size;
		}
	}

	size_t size() const
	{
		return size_;
	}
	bool empty() const
	{
		return size_ == 0;
	}

	T & operator [] (int i)
	{
		assert (i < (int) size_);
		return data_[i];
	}

	const T & operator [] (int i) const
	{
		assert (i < (int) size_);
		return data_[i];
	}
};

template < typename T > struct ArrayDevice;

template < typename T >
struct ArrayHost: public Array < T, std::allocator < T > >
{
	ArrayHost() : Array < T, std::allocator < T > > () {}
	ArrayHost (size_t size) : Array < T, std::allocator < T > > (size) {}

	ArrayHost<T> & operator = (const ArrayDevice < T > & other);
};

template < typename T >
struct ArrayDevice: public Array < T, Allocator < T > >
{
	ArrayDevice() : Array < T,  Allocator < T > > () {}
	ArrayDevice (size_t size) : Array < T,  Allocator < T > > (size) {}

	ArrayDevice<T> & operator = (const ArrayHost < T > & other);
};

template < typename T >
ArrayHost<T> & ArrayHost<T>::operator = (const ArrayDevice < T > & other)
{
	resize(other.size());
	vec_copy_from_device(&(*this)[0], &other[0], this->size());
}

template < typename T >
ArrayDevice<T> & ArrayDevice<T>::operator = (const ArrayHost < T > & other)
{
	resize(other.size());
	vec_copy_from_host(&(*this)[0], &other[0], this->size());
}

template < typename T, template<class> class Alloc1, template<class> class Alloc2 >
struct ArrayCopier
{
};

template < typename T, template<class> class Alloc1 >
struct ArrayCopier < T, Alloc1, Alloc1 >
{
	void copy(Array < T, Alloc1 < T > > & ar1, const Array < T, Alloc1 < T > > & ar2)
	{
		ar1 = ar2;
	}
};

template < typename T >
struct ArrayCopier < T, std::allocator, Allocator >
{
	void copy(Array < T, std::allocator < T > > & ar1, const Array < T, Allocator < T > > & ar2)
	{
		ar1.resize(ar2.size());
		vec_copy_from_device(&ar1[0], &ar2[0], ar1.size());
	}
};

template < typename T >
struct ArrayCopier < T, Allocator, std::allocator >
{
	void copy(Array < T, Allocator < T > > & ar1, const Array < T, std::allocator < T > > & ar2)
	{
		ar1.resize(ar2.size());
		vec_copy_from_host(&ar1[0], &ar2[0], ar1.size());
	}
};

template < typename T, template<class> class Alloc1, template<class> class Alloc2 >
void array_copy(Array < T, Alloc1 < T > > & ar1, const Array < T, Alloc2 < T > > & ar2)
{
	ArrayCopier < T, Alloc1, Alloc2 > copier;
	copier.copy(ar1, ar2);
}

/**
 * @ingroup misc
 * array.
 */
typedef Array < double, Allocator < double > > vec;

} /* namespace */

#endif /* PHELM_ARRAY_H */

