#pragma once
#include "tree_node.hpp"
#include <memory>

namespace dmc
{
	template <class Scalar>
	class branch_tree_node : public tree_node<Scalar>
	{
	public:
		typedef tree_node<Scalar> base_type;

		explicit branch_tree_node(std::array<std::unique_ptr<base_type>, 8>&& children)
			: children_(std::move(children))
		{
		}

		const std::array<std::unique_ptr<base_type>, 8>& children() const
		{
			return children_;
		}

	private:
		std::array<std::unique_ptr<base_type>, 8> children_;
	};
}
