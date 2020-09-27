#pragma once
#include "branch_tree_node.hpp"
#include "leaf_tree_node.hpp"
#include "marching_cubes.hpp"
#include "object.hpp"
#include "tree_config.hpp"
#include "tree_node.hpp"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <memory>
#include <vector>
#include <random>

namespace dmc
{
	template <class Scalar>
	class tree
	{
	public:
		typedef Scalar scalar_type;
		typedef Eigen::Matrix<scalar_type, 3, 1> vector_type;
		typedef Eigen::Matrix<std::size_t, 3, 1> vector_size_type;
		typedef dual<scalar_type, 3> dual_type;
		typedef object<scalar_type> object_type;
		typedef tree_config<scalar_type> config_type;
		typedef tree_node<scalar_type> node_type;
		typedef branch_tree_node<scalar_type> branch_node_type;
		typedef leaf_tree_node<scalar_type> leaf_node_type;
		typedef vertex<scalar_type> vertex_type;

		explicit tree(const vector_type& minimum, const vector_type& maximum, const config_type& config = config_type())
			: minimum_(minimum)
			, config_(config)
		{
			vector_type v = maximum - minimum;

			for (int i = 0; i < 3; ++i)
			{
				size_[i] = static_cast<std::size_t>(std::max(static_cast<scalar_type>(1.0), std::ceil(v[i] / config_.grid_width)));
			}

			children_.resize(size_.prod());
		}

		void generate(const object_type& obj, const std::function<void(double)>& progress_receiver)
		{
			generate(obj, &progress_receiver);
		}

		void generate(const object_type& obj, const std::function<void(double)>* progress_receiver = nullptr)
		{
			std::size_t progress = 0;

			if (progress_receiver)
			{
				(*progress_receiver)(0.0);
			}

			auto total_size = size_.x() * size_.y() * size_.z();

			std::vector<std::size_t> indices(total_size);

			for (std::size_t i = 0; i < total_size; ++i)
			{
				indices[i] = i;
			}

			std::mt19937 gen;
			std::shuffle(indices.begin(), indices.end(), gen);

#pragma omp parallel for schedule(dynamic)
			for (std::ptrdiff_t i = 0; i < static_cast<std::ptrdiff_t>(total_size); ++i)
			{
				auto j = indices[i];
				auto ix = j % size_.x();
				auto iy = j / size_.x() % size_.y();
				auto iz = j / size_.x() / size_.y();
				vector_type minimum = minimum_ + vector_size_type(ix, iy, iz).cast<scalar_type>() * config_.grid_width;
				vector_type maximum = minimum_ + vector_size_type(ix + 1, iy + 1, iz + 1).cast<scalar_type>() * config_.grid_width;

				children_[index(ix, iy, iz)] = generate_impl(obj, minimum, maximum, 0);

				if (progress_receiver)
				{
#pragma omp critical
					{
						(*progress_receiver)(static_cast<double>(++progress) / static_cast<double>(total_size));
					}
				}
			}
		}

		template <class Receiver>
		void enumerate(Receiver receiver)
		{
			for (std::size_t iz = 0; iz < size_.z(); ++iz)
			{
				for (std::size_t iy = 0; iy < size_.y(); ++iy)
				{
					for (std::size_t ix = 0; ix < size_.x(); ++ix)
					{
						enumerate_impl_c(*children_[index(ix, iy, iz)], receiver);

						if (ix != size_.x() - 1)
							enumerate_impl_f_x(*children_[index(ix, iy, iz)], *children_[index(ix + 1, iy, iz)], receiver);

						if (iy != size_.y() - 1)
							enumerate_impl_f_y(*children_[index(ix, iy, iz)], *children_[index(ix, iy + 1, iz)], receiver);

						if (iz != size_.z() - 1)
							enumerate_impl_f_z(*children_[index(ix, iy, iz)], *children_[index(ix, iy, iz + 1)], receiver);

						if (ix != size_.x() - 1 && iy != size_.y() - 1)
						{
							enumerate_impl_e_xy(
								*children_[index(ix, iy, iz)],
								*children_[index(ix + 1, iy, iz)],
								*children_[index(ix, iy + 1, iz)],
								*children_[index(ix + 1, iy + 1, iz)],
								receiver);
						}

						if (iy != size_.y() - 1 && iz != size_.z() - 1)
						{
							enumerate_impl_e_yz(
								*children_[index(ix, iy, iz)],
								*children_[index(ix, iy + 1, iz)],
								*children_[index(ix, iy, iz + 1)],
								*children_[index(ix, iy + 1, iz + 1)],
								receiver);
						}

						if (ix != size_.x() - 1 && iz != size_.z() - 1)
						{
							enumerate_impl_e_xz(
								*children_[index(ix, iy, iz)],
								*children_[index(ix + 1, iy, iz)],
								*children_[index(ix, iy, iz + 1)],
								*children_[index(ix + 1, iy, iz + 1)],
								receiver);
						}

						if (ix != size_.x() - 1 && iy != size_.y() - 1 && iz != size_.z() - 1)
						{
							enumerate_impl_v(
								*children_[index(ix, iy, iz)],
								*children_[index(ix + 1, iy, iz)],
								*children_[index(ix, iy + 1, iz)],
								*children_[index(ix + 1, iy + 1, iz)],
								*children_[index(ix, iy, iz + 1)],
								*children_[index(ix + 1, iy, iz + 1)],
								*children_[index(ix, iy + 1, iz + 1)],
								*children_[index(ix + 1, iy + 1, iz + 1)],
								receiver);
						}
					}
				}
			}
		}

	private:
		std::size_t index(std::size_t ix, std::size_t iy, std::size_t iz) const
		{
			return iz * size_.y() * size_.x() + iy * size_.x() + ix;
		}

		std::unique_ptr<node_type> generate_impl(const object_type& obj, const vector_type& minimum, const vector_type& maximum, std::size_t depth) const
		{
			std::array<vector_type, 8> points =
				{{
					{minimum.x(), minimum.y(), minimum.z()},
					{maximum.x(), minimum.y(), minimum.z()},
					{minimum.x(), maximum.y(), minimum.z()},
					{maximum.x(), maximum.y(), minimum.z()},
					{minimum.x(), minimum.y(), maximum.z()},
					{maximum.x(), minimum.y(), maximum.z()},
					{minimum.x(), maximum.y(), maximum.z()},
					{maximum.x(), maximum.y(), maximum.z()},
				}};

			std::array<dual_type, 8> values;

			auto sanitize = [](scalar_type x) {
				return std::isnan(x) ? 0.0 : x;
			};

			std::transform(points.begin(), points.end(), values.begin(), [&](const auto& p) {
				auto d = obj.value_grad(p);

				for (int i = 0; i < d.grad().rows(); ++i)
				{
					d.grad()(i) = sanitize(d.grad()(i));
				}

				return d;
			});

			Eigen::Matrix<scalar_type, 11, 4> a;

			for (int i = 0; i < 8; ++i)
			{
				a(i, 0) = values[i].grad()(0);
				a(i, 1) = values[i].grad()(1);
				a(i, 2) = values[i].grad()(2);
				a(i, 3) = static_cast<scalar_type>(-1.0);
			}

			a(8, 0) = config_.nominal_weight;
			a(8, 1) = 0.0;
			a(8, 2) = 0.0;
			a(8, 3) = 0.0;
			a(9, 0) = 0.0;
			a(9, 1) = config_.nominal_weight;
			a(9, 2) = 0.0;
			a(9, 3) = 0.0;
			a(10, 0) = 0.0;
			a(10, 1) = 0.0;
			a(10, 2) = config_.nominal_weight;
			a(10, 3) = 0.0;

			Eigen::Matrix<scalar_type, 11, 1> b;

			vector_type medium = (minimum + maximum) * static_cast<scalar_type>(0.5);

			for (int i = 0; i < 8; ++i)
			{
				b(i) = values[i].grad().dot(points[i] - medium) - values[i].value();
			}

			b(8) = 0.0;
			b(9) = 0.0;
			b(10) = 0.0;

			Eigen::Matrix<scalar_type, 4, 1> x = a.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV).solve(b);

			vector_type center = vector_type(x(0), x(1), x(2)) + medium;
			auto offset = obj.value(center);

			auto error = scalar_type();

			for (int i = 0; i < 8; ++i)
			{
				error += squared(offset - values[i].value() - values[i].grad().dot(center - points[i]));
			}

			if (depth >= config_.maximum_depth || error < squared(config_.tolerance))
			{
				return std::make_unique<leaf_node_type>(vertex_type(center, offset));
			}
			else
			{
				std::array<std::unique_ptr<node_type>, 8> nodes =
					{{
						generate_impl(obj, {minimum.x(), minimum.y(), minimum.z()}, {medium.x(), medium.y(), medium.z()}, depth + 1),
						generate_impl(obj, {medium.x(), minimum.y(), minimum.z()}, {maximum.x(), medium.y(), medium.z()}, depth + 1),
						generate_impl(obj, {minimum.x(), medium.y(), minimum.z()}, {medium.x(), maximum.y(), medium.z()}, depth + 1),
						generate_impl(obj, {medium.x(), medium.y(), minimum.z()}, {maximum.x(), maximum.y(), medium.z()}, depth + 1),
						generate_impl(obj, {minimum.x(), minimum.y(), medium.z()}, {medium.x(), medium.y(), maximum.z()}, depth + 1),
						generate_impl(obj, {medium.x(), minimum.y(), medium.z()}, {maximum.x(), medium.y(), maximum.z()}, depth + 1),
						generate_impl(obj, {minimum.x(), medium.y(), medium.z()}, {medium.x(), maximum.y(), maximum.z()}, depth + 1),
						generate_impl(obj, {medium.x(), medium.y(), medium.z()}, {maximum.x(), maximum.y(), maximum.z()}, depth + 1),
					}};

				return std::make_unique<branch_node_type>(std::move(nodes));
			}
		}

		template <class Receiver>
		void enumerate_impl_c(const node_type& n, Receiver receiver)
		{
			if (auto b = dynamic_cast<const branch_node_type*>(&n))
			{
				enumerate_impl_c(*b->children()[0], receiver);
				enumerate_impl_c(*b->children()[1], receiver);
				enumerate_impl_c(*b->children()[2], receiver);
				enumerate_impl_c(*b->children()[3], receiver);
				enumerate_impl_c(*b->children()[4], receiver);
				enumerate_impl_c(*b->children()[5], receiver);
				enumerate_impl_c(*b->children()[6], receiver);
				enumerate_impl_c(*b->children()[7], receiver);

				enumerate_impl_f_x(*b->children()[0], *b->children()[1], receiver);
				enumerate_impl_f_x(*b->children()[2], *b->children()[3], receiver);
				enumerate_impl_f_x(*b->children()[4], *b->children()[5], receiver);
				enumerate_impl_f_x(*b->children()[6], *b->children()[7], receiver);

				enumerate_impl_f_y(*b->children()[0], *b->children()[2], receiver);
				enumerate_impl_f_y(*b->children()[1], *b->children()[3], receiver);
				enumerate_impl_f_y(*b->children()[4], *b->children()[6], receiver);
				enumerate_impl_f_y(*b->children()[5], *b->children()[7], receiver);

				enumerate_impl_f_z(*b->children()[0], *b->children()[4], receiver);
				enumerate_impl_f_z(*b->children()[1], *b->children()[5], receiver);
				enumerate_impl_f_z(*b->children()[2], *b->children()[6], receiver);
				enumerate_impl_f_z(*b->children()[3], *b->children()[7], receiver);

				enumerate_impl_e_xy(*b->children()[0], *b->children()[1], *b->children()[2], *b->children()[3], receiver);
				enumerate_impl_e_xy(*b->children()[4], *b->children()[5], *b->children()[6], *b->children()[7], receiver);

				enumerate_impl_e_yz(*b->children()[0], *b->children()[2], *b->children()[4], *b->children()[6], receiver);
				enumerate_impl_e_yz(*b->children()[1], *b->children()[3], *b->children()[5], *b->children()[7], receiver);

				enumerate_impl_e_xz(*b->children()[0], *b->children()[1], *b->children()[4], *b->children()[5], receiver);
				enumerate_impl_e_xz(*b->children()[2], *b->children()[3], *b->children()[6], *b->children()[7], receiver);

				enumerate_impl_v(
					*b->children()[0],
					*b->children()[1],
					*b->children()[2],
					*b->children()[3],
					*b->children()[4],
					*b->children()[5],
					*b->children()[6],
					*b->children()[7],
					receiver);
			}
		}

		template <class Receiver>
		void enumerate_impl_f_x(const node_type& n1, const node_type& n2, Receiver receiver)
		{
			auto b1 = dynamic_cast<const branch_node_type*>(&n1);
			auto b2 = dynamic_cast<const branch_node_type*>(&n2);

			if (b1 || b2)
			{
				enumerate_impl_f_x(b1 ? *b1->children()[1] : n1, b2 ? *b2->children()[0] : n2, receiver);
				enumerate_impl_f_x(b1 ? *b1->children()[3] : n1, b2 ? *b2->children()[2] : n2, receiver);
				enumerate_impl_f_x(b1 ? *b1->children()[5] : n1, b2 ? *b2->children()[4] : n2, receiver);
				enumerate_impl_f_x(b1 ? *b1->children()[7] : n1, b2 ? *b2->children()[6] : n2, receiver);

				enumerate_impl_e_xy(
					b1 ? *b1->children()[1] : n1,
					b2 ? *b2->children()[0] : n2,
					b1 ? *b1->children()[3] : n1,
					b2 ? *b2->children()[2] : n2,
					receiver);

				enumerate_impl_e_xy(
					b1 ? *b1->children()[5] : n1,
					b2 ? *b2->children()[4] : n2,
					b1 ? *b1->children()[7] : n1,
					b2 ? *b2->children()[6] : n2,
					receiver);

				enumerate_impl_e_xz(
					b1 ? *b1->children()[1] : n1,
					b2 ? *b2->children()[0] : n2,
					b1 ? *b1->children()[5] : n1,
					b2 ? *b2->children()[4] : n2,
					receiver);

				enumerate_impl_e_xz(
					b1 ? *b1->children()[3] : n1,
					b2 ? *b2->children()[2] : n2,
					b1 ? *b1->children()[7] : n1,
					b2 ? *b2->children()[6] : n2,
					receiver);

				enumerate_impl_v(
					b1 ? *b1->children()[1] : n1,
					b2 ? *b2->children()[0] : n2,
					b1 ? *b1->children()[3] : n1,
					b2 ? *b2->children()[2] : n2,
					b1 ? *b1->children()[5] : n1,
					b2 ? *b2->children()[4] : n2,
					b1 ? *b1->children()[7] : n1,
					b2 ? *b2->children()[6] : n2,
					receiver);
			}
		}

		template <class Receiver>
		void enumerate_impl_f_y(const node_type& n1, const node_type& n2, Receiver receiver)
		{
			auto b1 = dynamic_cast<const branch_node_type*>(&n1);
			auto b2 = dynamic_cast<const branch_node_type*>(&n2);

			if (b1 || b2)
			{
				enumerate_impl_f_y(b1 ? *b1->children()[2] : n1, b2 ? *b2->children()[0] : n2, receiver);
				enumerate_impl_f_y(b1 ? *b1->children()[3] : n1, b2 ? *b2->children()[1] : n2, receiver);
				enumerate_impl_f_y(b1 ? *b1->children()[6] : n1, b2 ? *b2->children()[4] : n2, receiver);
				enumerate_impl_f_y(b1 ? *b1->children()[7] : n1, b2 ? *b2->children()[5] : n2, receiver);

				enumerate_impl_e_xy(
					b1 ? *b1->children()[2] : n1,
					b1 ? *b1->children()[3] : n1,
					b2 ? *b2->children()[0] : n2,
					b2 ? *b2->children()[1] : n2,
					receiver);

				enumerate_impl_e_xy(
					b1 ? *b1->children()[6] : n1,
					b1 ? *b1->children()[7] : n1,
					b2 ? *b2->children()[4] : n2,
					b2 ? *b2->children()[5] : n2,
					receiver);

				enumerate_impl_e_yz(
					b1 ? *b1->children()[3] : n1,
					b2 ? *b2->children()[1] : n2,
					b1 ? *b1->children()[7] : n1,
					b2 ? *b2->children()[5] : n2,
					receiver);

				enumerate_impl_e_yz(
					b1 ? *b1->children()[2] : n1,
					b2 ? *b2->children()[0] : n2,
					b1 ? *b1->children()[6] : n1,
					b2 ? *b2->children()[4] : n2,
					receiver);

				enumerate_impl_v(
					b1 ? *b1->children()[2] : n1,
					b1 ? *b1->children()[3] : n1,
					b2 ? *b2->children()[0] : n2,
					b2 ? *b2->children()[1] : n2,
					b1 ? *b1->children()[6] : n1,
					b1 ? *b1->children()[7] : n1,
					b2 ? *b2->children()[4] : n2,
					b2 ? *b2->children()[5] : n2,
					receiver);
			}
		}

		template <class Receiver>
		void enumerate_impl_f_z(const node_type& n1, const node_type& n2, Receiver receiver)
		{
			auto b1 = dynamic_cast<const branch_node_type*>(&n1);
			auto b2 = dynamic_cast<const branch_node_type*>(&n2);

			if (b1 || b2)
			{
				enumerate_impl_f_z(b1 ? *b1->children()[4] : n1, b2 ? *b2->children()[0] : n2, receiver);
				enumerate_impl_f_z(b1 ? *b1->children()[5] : n1, b2 ? *b2->children()[1] : n2, receiver);
				enumerate_impl_f_z(b1 ? *b1->children()[6] : n1, b2 ? *b2->children()[2] : n2, receiver);
				enumerate_impl_f_z(b1 ? *b1->children()[7] : n1, b2 ? *b2->children()[3] : n2, receiver);

				enumerate_impl_e_xz(
					b1 ? *b1->children()[4] : n1,
					b1 ? *b1->children()[5] : n1,
					b2 ? *b2->children()[0] : n2,
					b2 ? *b2->children()[1] : n2,
					receiver);

				enumerate_impl_e_xz(
					b1 ? *b1->children()[6] : n1,
					b1 ? *b1->children()[7] : n1,
					b2 ? *b2->children()[2] : n2,
					b2 ? *b2->children()[3] : n2,
					receiver);

				enumerate_impl_e_yz(
					b1 ? *b1->children()[4] : n1,
					b1 ? *b1->children()[6] : n1,
					b2 ? *b2->children()[0] : n2,
					b2 ? *b2->children()[2] : n2,
					receiver);

				enumerate_impl_e_yz(
					b1 ? *b1->children()[5] : n1,
					b1 ? *b1->children()[7] : n1,
					b2 ? *b2->children()[1] : n2,
					b2 ? *b2->children()[3] : n2,
					receiver);

				enumerate_impl_v(
					b1 ? *b1->children()[4] : n1,
					b1 ? *b1->children()[5] : n1,
					b1 ? *b1->children()[6] : n1,
					b1 ? *b1->children()[7] : n1,
					b2 ? *b2->children()[0] : n2,
					b2 ? *b2->children()[1] : n2,
					b2 ? *b2->children()[2] : n2,
					b2 ? *b2->children()[3] : n2,
					receiver);
			}
		}

		template <class Receiver>
		void enumerate_impl_e_xy(
			const node_type& n1,
			const node_type& n2,
			const node_type& n3,
			const node_type& n4,
			Receiver receiver)
		{
			auto b1 = dynamic_cast<const branch_node_type*>(&n1);
			auto b2 = dynamic_cast<const branch_node_type*>(&n2);
			auto b3 = dynamic_cast<const branch_node_type*>(&n3);
			auto b4 = dynamic_cast<const branch_node_type*>(&n4);

			if (b1 || b2 || b3 || b4)
			{
				enumerate_impl_e_xy(
					b1 ? *b1->children()[3] : n1,
					b2 ? *b2->children()[2] : n2,
					b3 ? *b3->children()[1] : n3,
					b4 ? *b4->children()[0] : n4,
					receiver);

				enumerate_impl_e_xy(
					b1 ? *b1->children()[7] : n1,
					b2 ? *b2->children()[6] : n2,
					b3 ? *b3->children()[5] : n3,
					b4 ? *b4->children()[4] : n4,
					receiver);

				enumerate_impl_v(
					b1 ? *b1->children()[3] : n1,
					b2 ? *b2->children()[2] : n2,
					b3 ? *b3->children()[1] : n3,
					b4 ? *b4->children()[0] : n4,
					b1 ? *b1->children()[7] : n1,
					b2 ? *b2->children()[6] : n2,
					b3 ? *b3->children()[5] : n3,
					b4 ? *b4->children()[4] : n4,
					receiver);
			}
		}

		template <class Receiver>
		void enumerate_impl_e_yz(
			const node_type& n1,
			const node_type& n2,
			const node_type& n3,
			const node_type& n4,
			Receiver receiver)
		{
			auto b1 = dynamic_cast<const branch_node_type*>(&n1);
			auto b2 = dynamic_cast<const branch_node_type*>(&n2);
			auto b3 = dynamic_cast<const branch_node_type*>(&n3);
			auto b4 = dynamic_cast<const branch_node_type*>(&n4);

			if (b1 || b2 || b3 || b4)
			{
				enumerate_impl_e_yz(
					b1 ? *b1->children()[6] : n1,
					b2 ? *b2->children()[4] : n2,
					b3 ? *b3->children()[2] : n3,
					b4 ? *b4->children()[0] : n4,
					receiver);

				enumerate_impl_e_yz(
					b1 ? *b1->children()[7] : n1,
					b2 ? *b2->children()[5] : n2,
					b3 ? *b3->children()[3] : n3,
					b4 ? *b4->children()[1] : n4,
					receiver);

				enumerate_impl_v(
					b1 ? *b1->children()[6] : n1,
					b1 ? *b1->children()[7] : n1,
					b2 ? *b2->children()[4] : n2,
					b2 ? *b2->children()[5] : n2,
					b3 ? *b3->children()[2] : n3,
					b3 ? *b3->children()[3] : n3,
					b4 ? *b4->children()[0] : n4,
					b4 ? *b4->children()[1] : n4,
					receiver);
			}
		}

		template <class Receiver>
		void enumerate_impl_e_xz(
			const node_type& n1,
			const node_type& n2,
			const node_type& n3,
			const node_type& n4,
			Receiver receiver)
		{
			auto b1 = dynamic_cast<const branch_node_type*>(&n1);
			auto b2 = dynamic_cast<const branch_node_type*>(&n2);
			auto b3 = dynamic_cast<const branch_node_type*>(&n3);
			auto b4 = dynamic_cast<const branch_node_type*>(&n4);

			if (b1 || b2 || b3 || b4)
			{
				enumerate_impl_e_xz(
					b1 ? *b1->children()[5] : n1,
					b2 ? *b2->children()[4] : n2,
					b3 ? *b3->children()[1] : n3,
					b4 ? *b4->children()[0] : n4,
					receiver);

				enumerate_impl_e_xz(
					b1 ? *b1->children()[7] : n1,
					b2 ? *b2->children()[6] : n2,
					b3 ? *b3->children()[3] : n3,
					b4 ? *b4->children()[2] : n4,
					receiver);

				enumerate_impl_v(
					b1 ? *b1->children()[5] : n1,
					b2 ? *b2->children()[4] : n2,
					b1 ? *b1->children()[7] : n1,
					b2 ? *b2->children()[6] : n2,
					b3 ? *b3->children()[1] : n3,
					b4 ? *b4->children()[0] : n4,
					b3 ? *b3->children()[3] : n3,
					b4 ? *b4->children()[2] : n4,
					receiver);
			}
		}

		template <class Receiver>
		void enumerate_impl_v(
			const node_type& n1,
			const node_type& n2,
			const node_type& n3,
			const node_type& n4,
			const node_type& n5,
			const node_type& n6,
			const node_type& n7,
			const node_type& n8,
			Receiver receiver)
		{
			auto b1 = dynamic_cast<const branch_node_type*>(&n1);
			auto b2 = dynamic_cast<const branch_node_type*>(&n2);
			auto b3 = dynamic_cast<const branch_node_type*>(&n3);
			auto b4 = dynamic_cast<const branch_node_type*>(&n4);
			auto b5 = dynamic_cast<const branch_node_type*>(&n5);
			auto b6 = dynamic_cast<const branch_node_type*>(&n6);
			auto b7 = dynamic_cast<const branch_node_type*>(&n7);
			auto b8 = dynamic_cast<const branch_node_type*>(&n8);

			if (b1 || b2 || b3 || b4 || b5 || b6 || b7 || b8)
			{
				enumerate_impl_v(
					b1 ? *b1->children()[7] : n1,
					b2 ? *b2->children()[6] : n2,
					b3 ? *b3->children()[5] : n3,
					b4 ? *b4->children()[4] : n4,
					b5 ? *b5->children()[3] : n5,
					b6 ? *b6->children()[2] : n6,
					b7 ? *b7->children()[1] : n7,
					b8 ? *b8->children()[0] : n8,
					receiver);
			}
			else
			{
				auto l1 = static_cast<const leaf_node_type*>(&n1);
				auto l2 = static_cast<const leaf_node_type*>(&n2);
				auto l3 = static_cast<const leaf_node_type*>(&n3);
				auto l4 = static_cast<const leaf_node_type*>(&n4);
				auto l5 = static_cast<const leaf_node_type*>(&n5);
				auto l6 = static_cast<const leaf_node_type*>(&n6);
				auto l7 = static_cast<const leaf_node_type*>(&n7);
				auto l8 = static_cast<const leaf_node_type*>(&n8);

				std::array<const vertex_type*, 8> vertices =
					{{
						&l1->vertex(),
						&l2->vertex(),
						&l3->vertex(),
						&l4->vertex(),
						&l5->vertex(),
						&l6->vertex(),
						&l7->vertex(),
						&l8->vertex(),
					}};

				marching_cubes<scalar_type>(vertices, receiver);
			}
		}

		vector_type minimum_;
		config_type config_;

		std::size_t grid_size_;
		scalar_type grid_width_;

		vector_size_type size_;
		std::vector<std::unique_ptr<node_type>> children_;
	};
}
