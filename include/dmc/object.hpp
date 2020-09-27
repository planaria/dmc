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
		typedef dual<scalar_type, 3> dual_type;

		virtual dual_type value_grad(const vector_type& p) const = 0;

		virtual scalar_type value(const vector_type& p) const = 0;
	};

	template <class Scalar, class Derived>
	class dual_object : public object<Scalar>
	{
	private:
		typedef object<Scalar> base_type;

	public:
		using typename base_type::scalar_type;
		using typename base_type::vector_type;
		using typename base_type::dual_type;

		virtual dual_type value_grad(const vector_type& p) const override
		{
			typedef dual<scalar_type, 3> dual_type;
			typedef vector<dual_type, 3> dual_vector_type;

			dual_type dx(p.x(), vector_type(static_cast<scalar_type>(1.0), static_cast<scalar_type>(0.0), static_cast<scalar_type>(0.0)));
			dual_type dy(p.y(), vector_type(static_cast<scalar_type>(0.0), static_cast<scalar_type>(1.0), static_cast<scalar_type>(0.0)));
			dual_type dz(p.z(), vector_type(static_cast<scalar_type>(0.0), static_cast<scalar_type>(0.0), static_cast<scalar_type>(1.0)));
			dual_vector_type dp(dx, dy, dz);

			return derived().templated_value(dp);
		}

		virtual scalar_type value(const vector_type& p) const override
		{
			return derived().templated_value(p);
		}

	private:
		const Derived& derived() const
		{
			return static_cast<const Derived&>(*this);
		}
	};
}
