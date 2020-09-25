#pragma once
#include "vector.hpp"
#include "dual.hpp"

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

	template <class Scalar, class Derived>
	class dual_object : public object<Scalar>
	{
	private:
		typedef object<Scalar> base_type;

	public:
		using typename base_type::scalar_type;
		using typename base_type::vector_type;

		virtual scalar_type value(const vector_type& p) const override
		{
			return derived().templated_value(p);
		}

		virtual vector_type grad(const vector_type& p) const override
		{
			typedef dual<scalar_type> dual_type;
			typedef vector<dual_type, 3> dual_vector_type;

			dual_vector_type dual_x(dual_type(p.x(), static_cast<scalar_type>(1.0)), dual_type(p.y(), static_cast<scalar_type>(0.0)), dual_type(p.z(), static_cast<scalar_type>(0.0)));
			dual_vector_type dual_y(dual_type(p.x(), static_cast<scalar_type>(0.0)), dual_type(p.y(), static_cast<scalar_type>(1.0)), dual_type(p.z(), static_cast<scalar_type>(0.0)));
			dual_vector_type dual_z(dual_type(p.x(), static_cast<scalar_type>(0.0)), dual_type(p.y(), static_cast<scalar_type>(0.0)), dual_type(p.z(), static_cast<scalar_type>(1.0)));

			auto px = derived().templated_value(dual_x);
			auto py = derived().templated_value(dual_y);
			auto pz = derived().templated_value(dual_z);

			return vector_type(px.grad(), py.grad(), pz.grad());
		}

	private:
		const Derived& derived() const
		{
			return static_cast<const Derived&>(*this);
		}
	};
}
