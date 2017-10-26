#pragma once
#include "vector.hpp"

namespace dmc
{
	template <class Scalar>
	class object
	{
	public:
		typedef Scalar scalar_type;
		typedef vector<scalar_type, 3> vector_type;

		virtual scalar_type value(const vector_type& p) const = 0;

		virtual vector_type grad(const vector_type& p) const = 0;
	};
}
