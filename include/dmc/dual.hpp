#pragma once
#include <boost/operators.hpp>

namespace dmc
{
	template <class Scalar>
	class dual
		: boost::addable<dual<Scalar>
		, boost::addable<dual<Scalar>, Scalar
		, boost::subtractable<dual<Scalar>
		, boost::subtractable<dual<Scalar>, Scalar
		, boost::multipliable<dual<Scalar>
		, boost::multipliable<dual<Scalar>, Scalar
		, boost::dividable<dual<Scalar>
		, boost::dividable<dual<Scalar>, Scalar
		, boost::equality_comparable<dual<Scalar>
		, boost::equality_comparable<dual<Scalar>, Scalar
		, boost::less_than_comparable<dual<Scalar>
		, boost::less_than_comparable<dual<Scalar>, Scalar
		>>>>>>>>>>>>
	{
	public:
		typedef Scalar scalar_type;

		dual() = default;

		dual(scalar_type value)
			: value_(value)
		{
		}

		dual(scalar_type value, scalar_type	grad)
			: value_(value)
			, grad_(grad)
		{
		}

		scalar_type value() const
		{
			return value_;
		}

		scalar_type grad() const
		{
			return grad_;
		}

		const dual& operator+() const
		{
			return *this;
		}

		dual operator-() const
		{
			return dual(-value_, -grad_);
		}

		dual& operator+=(const dual& rhs)
		{
			value_ += rhs.value_;
			grad_ += rhs.grad_;
			return *this;
		}

		dual& operator-=(const dual& rhs)
		{
			value_ -= rhs.value_;
			grad_ -= rhs.grad_;
			return *this;
		}

		dual& operator*=(const dual& rhs)
		{
			grad_ = value_ * rhs.grad_ + rhs.value_ * grad_;
			value_ *= rhs.value_;
			return *this;
		}

		dual& operator/=(const dual& rhs)
		{
			grad_ = (grad_ * rhs.value_ - value_ * rhs.grad_) / (rhs.value_ * rhs.value_);
			value_ /= rhs.value_;
			return *this;
		}

		friend bool operator==(const dual& lhs, const dual& rhs)
		{
			return lhs.value_ == rhs.value_;
		}

		friend bool operator<(const dual& lhs, const dual& rhs)
		{
			return lhs.value_ < rhs.value_;
		}

		friend dual abs(const dual& d)
		{
			return d.value() < static_cast<scalar_type>(0.0) ? -d : d;
		}

		friend dual sqrt(const dual& d)
		{
			auto value = std::sqrt(d.value());
			auto grad = value * static_cast<scalar_type>(0.5) * d.grad() / d.value();
			return dual(value, grad);
		}

	private:
		scalar_type value_{};
		scalar_type grad_{};
	};
}
