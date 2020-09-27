#include <dmc/dmc.hpp>
#include <fstream>

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

	template <class Vector>
	auto templated_value(const Vector& p) const
	{
		return 1.0 - p.norm_l_inf();
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
