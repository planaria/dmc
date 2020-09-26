#pragma once

namespace dmc
{
	template <class Scalar>
	Scalar clamp(Scalar x, Scalar minimum, Scalar maximum)
	{
		return std::min(std::max(x, minimum), maximum);
	}

	template <class Scalar>
	Scalar smooth_union(Scalar x1, Scalar x2, Scalar k)
	{
		auto h = clamp<Scalar>(0.5 + 0.5 * (x2 - x1) / k, 0.0, 1.0);
		return lerp(x2, x1, h) - k * h * (1.0 - h);
	}
}
