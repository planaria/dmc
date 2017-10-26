#pragma once
#include "math_util.hpp"
#include <boost/assert.hpp>
#include <boost/operators.hpp>
#include <cmath>
#include <limits>
#include <numeric>

namespace dmc
{
	template <class Scalar, int Dimension>
	class vector;

	template <class Derived, class Scalar, int Dimension>
	class vector_base
		: boost::addable<
			  Derived,
			  boost::subtractable<
				  Derived,
				  boost::multipliable<
					  Derived,
					  Scalar,
					  boost::dividable<
						  Derived,
						  Scalar,
						  boost::equality_comparable<
							  Derived>>>>>
	{
	public:
		typedef Scalar scalar_type;
		static const int dimension = Dimension;

		scalar_type* data()
		{
			return &values_[0];
		}

		const scalar_type* data() const
		{
			return &values_[0];
		}

		scalar_type& operator[](int index)
		{
			BOOST_ASSERT(0 <= index && index < dimension);
			return values_[index];
		}

		scalar_type operator[](int index) const
		{
			BOOST_ASSERT(0 <= index && index < dimension);
			return values_[index];
		}

		Derived operator-() const
		{
			return map([](auto x) { return -x; });
		}

		const Derived& operator+() const
		{
			return derived();
		}

		Derived& operator+=(const Derived& rhs)
		{
			for (int i = 0; i < dimension; ++i)
				(*this)[i] += rhs[i];
			return derived();
		}

		Derived& operator-=(const Derived& rhs)
		{
			for (int i = 0; i < dimension; ++i)
				(*this)[i] -= rhs[i];
			return derived();
		}

		Derived& operator*=(scalar_type rhs)
		{
			for (int i = 0; i < dimension; ++i)
				(*this)[i] *= rhs;
			return derived();
		}

		Derived& operator/=(scalar_type rhs)
		{
			for (int i = 0; i < dimension; ++i)
				(*this)[i] /= rhs;
			return derived();
		}

		template <class F>
		vector<typename std::result_of<F(scalar_type)>::type, dimension> map(F f) const
		{
			vector<typename std::result_of<F(scalar_type)>::type, dimension> result;
			for (int i = 0; i < dimension; ++i)
				result[i] = f((*this)[i]);
			return result;
		}

		template <class T>
		auto cast() const
		{
			return map([](auto x) { return static_cast<T>(x); });
		}

		auto sign() const
		{
			return map([](auto x) { return dmc::sign(x); });
		}

		template <class T, class F>
		auto reduce(T t, F f) const
		{
			for (int i = 0; i < dimension; ++i)
				t = f(t, (*this)[i]);
			return t;
		}

		template <class F>
		auto reduce(F f) const
		{
			auto t = (*this)[0];
			for (int i = 1; i < dimension; ++i)
				t = f(t, (*this)[i]);
			return t;
		}

		auto sum() const
		{
			return reduce([](auto x, auto y) {
				return x + y;
			});
		}

		auto product() const
		{
			return reduce([](auto x, auto y) {
				return x * y;
			});
		}

		auto abs() const
		{
			return map([](auto x) { using std::abs; return abs(x); });
		}

		auto norm_l1() const
		{
			return abs().sum();
		}

		auto squared() const
		{
			return map([](auto x) { return squared(x); });
		}

		auto norm_l2_sq() const
		{
			return squared().sum();
		}

		auto norm_l2() const
		{
			using std::sqrt;
			return sqrt(norm_l2_sq());
		}

		auto max() const
		{
			return reduce([](auto x, auto y) { using std::max; return max(x, y); });
		}

		auto min() const
		{
			return reduce([](auto x, auto y) { using std::min; return min(x, y); });
		}

		bool try_normalize()
		{
			auto n = norm_l2();
			if (n < std::numeric_limits<scalar_type>::epsilon())
				return false;

			*this /= n;
			return true;
		}

		Derived clamp(const Derived& minimum, const Derived& maximum)
		{
			Derived result;
			for (int i = 0; i < dimension; ++i)
				result[i] = std::max(minimum[i], std::min(maximum[i], (*this)[i]));
			return result;
		}

		static Derived all(scalar_type s)
		{
			Derived result;
			for (int i = 0; i < dimension; ++i)
				result[i] = s;
			return result;
		}

		friend bool operator==(const Derived& lhs, const Derived& rhs)
		{
			for (int i = 0; i < dimension; ++i)
				if (lhs[i] != rhs[i])
					return false;

			return true;
		}

	private:
		Derived& derived()
		{
			return static_cast<Derived&>(*this);
		}

		const Derived& derived() const
		{
			return static_cast<const Derived&>(*this);
		}

		scalar_type values_[dimension] = {};
	};
}
