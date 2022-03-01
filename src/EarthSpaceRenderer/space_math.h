#pragma once

namespace Space
{
	struct vec3d
	{
		union
		{
			struct { double x, y, z; };
			struct { double v[3]; };
		};

		vec3d(double _x = 0.0, double _y = 0.0, double _z = 0.0) : x(_x), y(_y), z(_z) {}
		vec3d(const double _v[]) { for (int i = 0; i < 3; ++i) { v[i] = _v[i]; } }

		double sqr() const;
		double abs() const;
		struct vec3d normalized() const;

		bool isfinite() const;
		bool isinf() const;
		bool isnan() const;
		bool isnormal() const;
	};

	struct m3x3d
	{
		union
		{
			struct { double xx, xy, xz, yx, yy, yz, zx, zy, zz; };
			struct { double v[9]; };
			struct { double m[3][3]; };
		};

		m3x3d(double _xx = 0.0, double _xy = 0.0, double _xz = 0.0
			, double _yx = 0.0, double _yy = 0.0, double _yz = 0.0
			, double _zx = 0.0, double _zy = 0.0, double _zz = 0.0)
			: xx(_xx), xy(_xy), xz(_xz)
			, yx(_yx), yy(_yy), yz(_yz)
			, zx(_zx), zy(_zy), zz(_zz) {}
		m3x3d(double const _v[]) { for (int i = 0; i < 9; ++i) { v[i] = _v[i]; } }

		double det() const;
		double norm1() const;
		double normesq() const;
		double norme() const;
		double normi() const;
		struct m3x3d transposed() const;
		struct m3x3d inverse() const;
	};

	struct m3x3d munit();
	struct m3x3d mrotX(double fi);
	struct m3x3d mrotY(double fi);
	struct m3x3d mrotZ(double fi);
	struct m3x3d mrotV(double fi, const struct vec3d &v);
}

struct Space::vec3d operator*(double a, const struct Space::vec3d &v);
struct Space::vec3d operator*(const struct Space::vec3d &v, double a);
struct Space::vec3d operator/(const struct Space::vec3d &v, double a);

double operator*(const struct Space::vec3d &v, const struct Space::vec3d &w);

struct Space::vec3d operator+(const struct Space::vec3d &v, const struct Space::vec3d &w);
struct Space::vec3d operator-(const struct Space::vec3d &v, const struct Space::vec3d &w);
struct Space::vec3d operator+(const struct Space::vec3d &v);
struct Space::vec3d operator-(const struct Space::vec3d &v);

struct Space::m3x3d operator*(double a, const struct Space::m3x3d &m);
struct Space::m3x3d operator*(const struct Space::m3x3d &m, double a);
struct Space::m3x3d operator/(const struct Space::m3x3d &m, double a);

struct Space::vec3d operator*(const struct Space::vec3d &v, const struct Space::m3x3d &m);
struct Space::vec3d operator*(const struct Space::m3x3d &m, const struct Space::vec3d &v);
struct Space::m3x3d operator*(const struct Space::m3x3d &m, const struct Space::m3x3d &n);

struct Space::m3x3d operator+(const struct Space::m3x3d &m, const struct Space::m3x3d &n);
struct Space::m3x3d operator-(const struct Space::m3x3d &m, const struct Space::m3x3d &n);
struct Space::m3x3d operator+(const struct Space::m3x3d &m);
struct Space::m3x3d operator-(const struct Space::m3x3d &m);
