#pragma once
#include "tree_node.hpp"
#include "vertex.hpp"

namespace dmc
{
	template <class Scalar>
	class leaf_tree_node : public tree_node<Scalar>
	{
	public:
		typedef tree_node<Scalar> base_type;

		using typename base_type::scalar_type;
		typedef dmc::vertex<scalar_type> vertex_type;

		explicit leaf_tree_node(const vertex_type& vertex)
			: vertex_(vertex)
		{
		}

		const vertex_type& vertex() const
		{
			return vertex_;
		}

	private:
		vertex_type vertex_;
	};
}
