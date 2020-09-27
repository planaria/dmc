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
	T templated_value(const Eigen::Matrix<T, 3, 1>& p) const
	{
		auto cube1 = 1.0 - p.cwiseAbs().maxCoeff();
		auto cube2 = 1.0 - (p - Eigen::Matrix<T, 3, 1>(0.5, 0.5, 0.5)).cwiseAbs().maxCoeff();

		return std::min(cube1, -cube2);
	}

private:
	scalar_type radius_;
};

int main(int /*argc*/, char* /*argv*/ [])
{
	Eigen::initParallel();

	dmc::tree_config<double> config;
	config.grid_width = 0.1;
	config.tolerance = 0.001;

	dmc::tree<double> t({-3.0, -3.0, -3.0}, {3.0, 3.0, 3.0}, config);

	double last_progress = 0.0;

	t.generate(test_object<double>(1.5f), [&](double progress) {
		if (progress > last_progress + 0.01)
		{
			std::cout << std::fixed << std::setprecision(3) << progress << std::endl;
			last_progress = progress;
		}
	});

	std::vector<dmc::triangle<Eigen::Vector3d>> triangles;

	t.enumerate([&](const auto& t) {
		triangles.push_back(t);
	});

	std::ofstream os("a.stl", std::ios::binary);

	write_stl(os, triangles);
}
