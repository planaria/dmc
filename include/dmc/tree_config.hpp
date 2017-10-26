#pragma once
#include <cstddef>

namespace dmc
{
	template <class Scalar>
	struct tree_config
	{
		typedef Scalar scalar_type;

		scalar_type grid_width = static_cast<scalar_type>(1.0);
		scalar_type tolerance = static_cast<scalar_type>(0.1);
		std::size_t maximum_depth = static_cast<std::size_t>(-1);
		scalar_type nominal_weight = static_cast<scalar_type>(0.1);
	};
}
