#pragma once
#include "vector.hpp"
#include <boost/noncopyable.hpp>

namespace dmc
{
	template <class Scalar>
	class tree_node : boost::noncopyable
	{
	public:
		typedef Scalar scalar_type;
		typedef vector<scalar_type, 3> vector_type;

		virtual ~tree_node()
		{
		}
	};
}
