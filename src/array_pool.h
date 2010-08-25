#ifndef ARRAY_POOL_H
#define ARRAY_POOL_H

#include <list>

template < typename T >
class ArrayPool
{
	typedef std::list < T > allocated_t;
	allocated_t allocated;
	typename allocated_t::iterator it;
	typename allocated_t::iterator cp;
	int n;

public:
	class Checkpoint
	{
		ArrayPool < T > & pool;

	public:
		Checkpoint(ArrayPool < T > & what): pool(what)
		{
			pool.checkpoint();
		}

		~Checkpoint()
		{
			pool.rollback();
		}
	};

	ArrayPool(int elems): it(allocated.begin()), n(elems) {}

	T & create() 
	{
		if (it == allocated.end())
		{
			if (cp == allocated.end()) {
				allocated.push_back(T(n));
				cp = allocated.end();
				--cp;
			} else {
				allocated.push_back(T(n));
			}
			it = allocated.end();
			return allocated.back();
		}

		return *it++;
	}

	void checkpoint()
	{
		cp = it;
	}

	void rollback()
	{
		it = cp;
	}
};

#endif /* ARRAY_POOL_H */

