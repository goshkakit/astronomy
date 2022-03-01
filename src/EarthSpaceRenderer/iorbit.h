#pragma once

namespace Scene
{
	struct vec3d_t
	{
		double x, y, z;

		vec3d_t(double _x = 0.0, double _y = 0.0, double _z = 0.0) : x(_x), y(_y), z(_z) {}
		vec3d_t(const double _v[]) : x(_v[0]), y(_v[1]), z(_v[2]) {}
	};

	class IOrbit
	{
	public:
		/*!
		\return new position, 10^3 km
		\param[in] time new time moment, min
		*/
		virtual struct Scene::vec3d_t GetNewPosition(double time) = 0;
	};
}
