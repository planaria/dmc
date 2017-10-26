#pragma once
#include "vector.hpp"

namespace dmc
{
	template <class Vector>
	class triangle
	{
	public:
		typedef Vector vector_type;

		triangle()
		{
		}

		triangle(const vector_type& p1, const vector_type& p2, const vector_type& p3)
			: p1_(p1)
			, p2_(p2)
			, p3_(p3)
		{
		}

		const vector_type& p1() const
		{
			return p1_;
		}

		const vector_type& p2() const
		{
			return p2_;
		}

		const vector_type& p3() const
		{
			return p3_;
		}

	private:
		vector_type p1_;
		vector_type p2_;
		vector_type p3_;
	};

	template <class Vector>
	auto make_triangle(const Vector& p1, const Vector& p2, const Vector& p3)
	{
		return triangle<Vector>(p1, p2, p3);
	}

	typedef triangle<vector3d> triangle3d;
	typedef triangle<vector3f> triangle3f;
}
