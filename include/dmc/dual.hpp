#pragma once
#include "vector.hpp"
#include <boost/operators.hpp>

namespace dmc
{
	template <class Scalar, int Dimension>
	class dual
		: boost::addable<dual<Scalar, Dimension>
		, boost::addable<dual<Scalar, Dimension>, Scalar
		, boost::subtractable<dual<Scalar, Dimension>
		, boost::subtractable<dual<Scalar, Dimension>, Scalar
		, boost::multipliable<dual<Scalar, Dimension>
		, boost::multipliable<dual<Scalar, Dimension>, Scalar
		, boost::dividable<dual<Scalar, Dimension>
		, boost::dividable<dual<Scalar, Dimension>, Scalar
		, boost::equality_comparable<dual<Scalar, Dimension>
		, boost::equality_comparable<dual<Scalar, Dimension>, Scalar
		, boost::less_than_comparable<dual<Scalar, Dimension>
		, boost::less_than_comparable<dual<Scalar, Dimension>, Scalar
		>>>>>>>>>>>>
	{
	public:
		typedef Scalar scalar_type;
		static const int dimension = Dimension;
		typedef vector<scalar_type, dimension> vector_type;

		dual() = default;

		dual(scalar_type value)
			: value_(value)
		{
		}

		dual(scalar_type value, const vector_type& grad)
			: value_(value)
			, grad_(grad)
		{
		}

		scalar_type& value()
		{
			return value_;
		}

		const scalar_type& value() const
		{
			return value_;
		}

		vector_type& grad()
		{
			return grad_;
		}

		const vector_type& grad() const
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

		friend dual sin(const dual& d)
		{
			auto value = std::sin(d.value());
			auto grad = d.grad() * std::cos(d.value());
			return dual(value, grad);
		}

		friend dual cos(const dual& d)
		{
			auto value = std::cos(d.value());
			auto grad = -d.grad() * std::sin(d.value());
			return dual(value, grad);
		}

	private:
		scalar_type value_{};
		vector_type grad_;
	};
}
