#include <dmc/dmc.hpp>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>

template <class Scalar>
struct test_object : dmc::dual_object<Scalar, test_object<Scalar>>
{
private:
	typedef dmc::object<Scalar> base_type;

public:
	using typename base_type::scalar_type;
	using typename base_type::vector_type;

	template <class T>
	T templated_value(const Eigen::Matrix<T, 3, 1>& p) const
	{
		auto cube1 = 1.0 - p.cwiseAbs().maxCoeff();
		auto cube2 = 1.0 - (p - Eigen::Matrix<T, 3, 1>(0.5, 0.5, 0.5)).cwiseAbs().maxCoeff();

		return std::min(cube1, -cube2);
	}
};

int main(int /*argc*/, char* /*argv*/ [])
{
	dmc::tree_config<double> config;
	config.grid_width = 0.1;
	config.tolerance = 0.001;

	auto start_time = std::chrono::high_resolution_clock::now();

	dmc::tree<double> t({-3.0, -3.0, -3.0}, {3.0, 3.0, 3.0}, config);

	double last_progress = 0.0;

	test_object<double> obj;

	t.generate(obj, [&](double progress) {
		if (progress > last_progress + 0.01)
		{
			std::cout << std::fixed << std::setprecision(3) << progress << std::endl;
			last_progress = progress;
		}
	});

	auto generated_time = std::chrono::high_resolution_clock::now();

	auto triangles = t.enumerate();

	auto enumerated_time = std::chrono::high_resolution_clock::now();

	std::cout << "Generation: " << std::fixed << std::setprecision(2) << std::chrono::duration<double>(generated_time - start_time).count() << "s" << std::endl;
	std::cout << "Enumeration: " << std::fixed << std::setprecision(2) << std::chrono::duration<double>(enumerated_time - generated_time).count() << "s" << std::endl;

	std::ofstream os("a.stl", std::ios::binary);
	write_stl(os, triangles);
}
