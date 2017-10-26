#pragma once
#include "vector_base.hpp"
#include <iostream>

namespace dmc
{
	template <class Scalar, int Dimension>
	class vector : public vector_base<vector<Scalar, Dimension>, Scalar, Dimension>
	{
	};

	template <class Scalar>
	class vector<Scalar, 2> : public vector_base<vector<Scalar, 2>, Scalar, 2>
	{
	private:
		typedef vector_base<vector<Scalar, 2>, Scalar, 2> base_type;

	public:
		using base_type::dimension;
		using typename base_type::scalar_type;

		vector()
		{
		}

		vector(scalar_type x_, scalar_type y_)
		{
			x() = x_;
			y() = y_;
		}

		scalar_type x() const
		{
			return (*this)[0];
		}

		scalar_type& x()
		{
			return (*this)[0];
		}

		scalar_type y() const
		{
			return (*this)[1];
		}

		scalar_type& y()
		{
			return (*this)[1];
		}
	};

	template <class Scalar>
	class vector<Scalar, 3> : public vector_base<vector<Scalar, 3>, Scalar, 3>
	{
	private:
		typedef vector_base<vector<Scalar, 3>, Scalar, 3> base_type;

	public:
		using base_type::dimension;
		using typename base_type::scalar_type;

		vector()
		{
		}

		vector(scalar_type x_, scalar_type y_, scalar_type z_)
		{
			x() = x_;
			y() = y_;
			z() = z_;
		}

		scalar_type& x()
		{
			return (*this)[0];
		}

		scalar_type x() const
		{
			return (*this)[0];
		}

		scalar_type& y()
		{
			return (*this)[1];
		}

		scalar_type y() const
		{
			return (*this)[1];
		}

		scalar_type& z()
		{
			return (*this)[2];
		}

		scalar_type z() const
		{
			return (*this)[2];
		}
	};

	template <class Scalar>
	class vector<Scalar, 4> : public vector_base<vector<Scalar, 4>, Scalar, 4>
	{
	private:
		typedef vector_base<vector<Scalar, 4>, Scalar, 4> base_type;

	public:
		using base_type::dimension;
		using typename base_type::scalar_type;

		vector()
		{
		}

		vector(scalar_type x_, scalar_type y_, scalar_type z_, scalar_type w_)
		{
			x() = x_;
			y() = y_;
			z() = z_;
			w() = w_;
		}

		vector(const vector<scalar_type, 3>& xyz_, scalar_type w_)
		{
			xyz() = xyz_;
			w() = w_;
		}

		scalar_type& x()
		{
			return (*this)[0];
		}

		scalar_type x() const
		{
			return (*this)[0];
		}

		scalar_type& y()
		{
			return (*this)[1];
		}

		scalar_type y() const
		{
			return (*this)[1];
		}

		scalar_type& z()
		{
			return (*this)[2];
		}

		scalar_type z() const
		{
			return (*this)[2];
		}

		scalar_type& w()
		{
			return (*this)[3];
		}

		scalar_type w() const
		{
			return (*this)[3];
		}

		vector<Scalar, 3>& xyz()
		{
			return *reinterpret_cast<vector<Scalar, 3>*>(this->data());
		}

		const vector<Scalar, 3>& xyz() const
		{
			return *reinterpret_cast<const vector<Scalar, 3>*>(this->data());
		}
	};

	template <class Scalar1, class Scalar2, int Dimension, class F>
	auto combine(const vector<Scalar1, Dimension>& v1, const vector<Scalar2, Dimension>& v2, F f)
	{
		typedef typename std::result_of<F(Scalar1, Scalar2)>::type result_scalar_type;

		vector<result_scalar_type, Dimension> result;

		for (int i = 0; i < Dimension; ++i)
			result[i] = f(v1[i], v2[i]);

		return result;
	}

	template <class Scalar, int Dimension>
	vector<Scalar, Dimension> minimum(const vector<Scalar, Dimension>& v1, const vector<Scalar, Dimension>& v2)
	{
		return combine(v1, v2, [](auto x, auto y) { return std::min(x, y); });
	}

	template <class Scalar, int Dimension>
	vector<Scalar, Dimension> maximum(const vector<Scalar, Dimension>& v1, const vector<Scalar, Dimension>& v2)
	{
		return combine(v1, v2, [](auto x, auto y) { return std::max(x, y); });
	}

	template <class Scalar, int Dimension>
	Scalar dot_product(const vector<Scalar, Dimension>& v1, const vector<Scalar, Dimension>& v2)
	{
		Scalar result = Scalar();

		for (int i = 0; i < Dimension; ++i)
			result += v1[i] * v2[i];

		return result;
	}

	template <class Scalar>
	Scalar cross_product(const vector<Scalar, 2>& v1, const vector<Scalar, 2>& v2)
	{
		return v1.x() * v2.y() - v1.y() * v2.x();
	}

	template <class Scalar>
	vector<Scalar, 3> cross_product(const vector<Scalar, 3>& v1, const vector<Scalar, 3>& v2)
	{
		return vector<Scalar, 3>(
			v1.y() * v2.z() - v1.z() * v2.y(),
			v1.z() * v2.x() - v1.x() * v2.z(),
			v1.x() * v2.y() - v1.y() * v2.x());
	}

	template <class Scalar, int Dimension>
	std::ostream& operator<<(std::ostream& os, const vector<Scalar, Dimension>& v)
	{
		for (int i = 0; i < Dimension; ++i)
		{
			if (i != 0)
				os << ", ";

			os << v[i];
		}

		return os;
	}

	typedef vector<double, 3> vector3d;
	typedef vector<float, 3> vector3f;
}
