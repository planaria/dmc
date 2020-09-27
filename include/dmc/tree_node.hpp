#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>
#include <boost/noncopyable.hpp>

namespace dmc
{
	template <class Scalar>
	class tree_node : boost::noncopyable
	{
	public:
		typedef Scalar scalar_type;
		typedef Eigen::Matrix<scalar_type, 3, 1> vector_type;

		virtual ~tree_node()
		{
		}
	};
}
