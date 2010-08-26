#ifndef ARRAY_POOL_H
#define ARRAY_POOL_H

#include <list>

template < typename T >
class ArrayPool
{
	typedef std::list < T > allocated_t;
	typedef std::list < typename allocated_t::iterator > checkpoints_t;
	allocated_t allocated;
	checkpoints_t checkpoints;
	typename allocated_t::iterator it;
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
			typename allocated_t::iterator prev;
			allocated.push_back(T(n));
			it = allocated.end();
			prev = it; --prev;
			for (typename checkpoints_t::iterator i = checkpoints.begin();
			     i != checkpoints.end(); ++i)
			{
				if (*i == allocated.end()) {
					*i = prev;
				}
			}
			return allocated.back();
		}

		return *it++;
	}

	void checkpoint()
	{
		checkpoints.push_back(it);
	}

	void rollback()
	{
		it = checkpoints.back();
		checkpoints.pop_back();
	}
};

#endif /* ARRAY_POOL_H */

