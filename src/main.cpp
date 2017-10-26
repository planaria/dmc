#include <dmc/dmc.hpp>
#include <fstream>

template <class Scalar>
struct test_object : dmc::object<Scalar>
{
public:
	typedef dmc::object<Scalar> base_type;

	using typename base_type::scalar_type;
	using typename base_type::vector_type;

	explicit test_object(scalar_type radius)
		: radius_(radius)
	{
	}

	virtual scalar_type value(const vector_type& p) const override
	{
		auto abs_p = p.map([](auto x) { return std::abs(x); });

		if (abs_p.x() < abs_p.y())
		{
			if (abs_p.y() < abs_p.z())
				return radius_ - abs_p.z();
			else
				return radius_ - abs_p.y();
		}
		else
		{
			if (abs_p.x() < abs_p.z())
				return radius_ - abs_p.z();
			else
				return radius_ - abs_p.x();
		}
	}

	virtual vector_type grad(const vector_type& p) const override
	{
		auto eps = 1.0e-6;

		return vector_type(
				   value(p + vector_type(eps, 0.0, 0.0)) - value(p - vector_type(eps, 0.0, 0.0)),
				   value(p + vector_type(0.0, eps, 0.0)) - value(p - vector_type(0.0, eps, 0.0)),
				   value(p + vector_type(0.0, 0.0, eps)) - value(p - vector_type(0.0, 0.0, eps))) /
			   (2.0 * eps);
	}

private:
	scalar_type radius_;
};

int main(int /*argc*/, char* /*argv*/ [])
{
	dmc::tree<double> t({-3.0, -3.0, -3.0}, {3.0, 3.0, 3.0});
	t.generate(test_object<double>(1.5f));

	std::vector<dmc::triangle3d> triangles;

	t.enumerate([&](const auto& t) {
		triangles.push_back(t);
	});

	std::ofstream os("a.stl", std::ios::binary);

	write_stl(os, triangles);
}
