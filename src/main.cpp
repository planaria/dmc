#include <dmc/dmc.hpp>
#include <iostream>
#include <fstream>
#include <iomanip>

template <class Scalar>
struct test_object : dmc::dual_object<Scalar, test_object<Scalar>>
{
public:
	typedef dmc::object<Scalar> base_type;

	using typename base_type::scalar_type;
	using typename base_type::vector_type;

	explicit test_object(scalar_type radius)
		: radius_(radius)
	{
	}

	template <class T>
	T templated_value(const dmc::vector<T, 3>& p) const
	{
		auto cube1 = 1.0 - p.norm_l_inf();
		auto cube2 = 1.0 - (p - dmc::vector<T, 3>(0.5, 0.5, 0.5)).norm_l_inf();

		return std::min(cube1, -cube2);
	}

private:
	scalar_type radius_;
};

int main(int /*argc*/, char* /*argv*/ [])
{
	dmc::tree_config<double> config;
	config.tolerance = 0.001;

	dmc::tree<double> t({-3.0, -3.0, -3.0}, {3.0, 3.0, 3.0}, config);

	t.generate(test_object<double>(1.5f), [](double progress) {
		std::cout << std::fixed << std::setprecision(3) << progress << std::endl;
	});

	std::vector<dmc::triangle3d> triangles;

	t.enumerate([&](const auto& t) {
		triangles.push_back(t);
	});

	std::ofstream os("a.stl", std::ios::binary);

	write_stl(os, triangles);
}
