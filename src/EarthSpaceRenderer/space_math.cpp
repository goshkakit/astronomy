#include "space_math.h"

#include <math.h>

struct Space::vec3d operator*(double a, const struct Space::vec3d &v)
{
	return Space::vec3d(a*v.x, a*v.y, a*v.z);
}

struct Space::vec3d operator*(const struct Space::vec3d &v, double a)
{
	return Space::vec3d(v.x*a, v.y*a, v.z*a);
}

struct Space::vec3d operator/(const struct Space::vec3d &v, double a)
{
	return Space::vec3d(v.x / a, v.y / a, v.z / a);
}

double operator*(const struct Space::vec3d &v, const struct Space::vec3d &w)
{
	return (v.x*w.x) + (v.y*w.y) + (v.z*w.z);
}

struct Space::vec3d operator+(const struct Space::vec3d &v, const struct Space::vec3d &w)
{
	return Space::vec3d(v.x + w.x, v.y + w.y, v.z + w.z);
}

struct Space::vec3d operator-(const struct Space::vec3d &v, const struct Space::vec3d &w)
{
	return Space::vec3d(v.x - w.x, v.y - w.y, v.z - w.z);
}

struct Space::vec3d operator+(const struct Space::vec3d &v)
{
	return Space::vec3d(v.x, v.y, v.z);
}

struct Space::vec3d operator-(const struct Space::vec3d &v)
{
	return Space::vec3d(-v.x, -v.y, -v.z);
}

double Space::vec3d::sqr() const
{
	return (x*x) + (y*y) + (z*z);
}

double Space::vec3d::abs() const
{
	return sqrt(this->sqr());
}

struct Space::vec3d Space::vec3d::normalized() const
{
	double a = this->abs();

	return Space::vec3d(x / a, y / a, z / a);
}

bool Space::vec3d::isfinite() const
{
	return ::isfinite(x) && ::isfinite(y) && ::isfinite(z);
}

bool Space::vec3d::isinf() const
{
	return ::isinf(x) || ::isinf(y) || ::isinf(z);
}

bool Space::vec3d::isnormal() const
{
	return ::isnormal(x) && ::isnormal(y) && ::isnormal(z);
}

bool Space::vec3d::isnan() const
{
	return ::isnan(x) || ::isnan(y) || ::isnan(z);
}

double Space::m3x3d::det() const
{
	return (xx*yy*zz) - (xx*yz*zy) - (xy*yx*zz) + (xy*yz*zx) + (xz*yx*zy) - (xz*yy*zx);
}

double Space::m3x3d::norm1() const
{
	return fmax(fmax(fabs(xx) + fabs(yx) + fabs(zx), fabs(xy) + fabs(yy) + fabs(zy)), fabs(xz) + fabs(yz) + fabs(zz));
}

double Space::m3x3d::normesq() const
{
	return (xx*xx) + (xy*xy) + (xz*xz) + (yx*yx) + (yy*yy) + (yz*yz) + (zx*zx) + (zy*zy) + (zz*zz);
}

double Space::m3x3d::norme() const
{
	return sqrt(normesq());
}

double Space::m3x3d::normi() const
{
	return fmax(fmax(fmax(fmax(fabs(xx), fabs(xy)), fmax(fabs(xz), fabs(yx))), fmax(fmax(fabs(yy), fabs(yz)), fmax(fabs(zx), fabs(zy)))), fabs(zz));
}

struct Space::m3x3d Space::m3x3d::transposed() const
{
	return Space::m3x3d(
		xx, yx, zx,
		xy, yy, zy,
		xz, yz, zz
	);
}

struct Space::m3x3d Space::m3x3d::inverse() const
{
	double d = det();

	return Space::m3x3d(
		(yy*zz - yz*zy) / d, (xz*zy - xy*zz) / d, (xy*yz - xz*yy) / d,
		(yz*zx - yx*zz) / d, (xx*zz - xz*zx) / d, (xz*yx - xx*yz) / d,
		(yx*zy - yy*zx) / d, (xy*zx - xx*zy) / d, (xx*yy - xy*yx) / d
	);
}

struct Space::m3x3d operator*(double a, const struct Space::m3x3d &m)
{
	return Space::m3x3d(
		a*m.xx, a*m.xy, a*m.xz,
		a*m.yx, a*m.yy, a*m.yz,
		a*m.zx, a*m.zy, a*m.zz
	);
}

struct Space::m3x3d operator*(const struct Space::m3x3d &m, double a)
{
	return Space::m3x3d(
		m.xx*a, m.xy*a, m.xz*a,
		m.yx*a, m.yy*a, m.yz*a,
		m.zx*a, m.zy*a, m.zz*a
	);
}

struct Space::m3x3d operator/(const struct Space::m3x3d &m, double a)
{
	return Space::m3x3d(
		m.xx / a, m.xy / a, m.xz / a,
		m.yx / a, m.yy / a, m.yz / a,
		m.zx / a, m.zy / a, m.zz / a
	);
}

struct Space::vec3d operator*(const struct Space::vec3d &v, const struct Space::m3x3d &m)
{
	return Space::vec3d(
		v.x*m.xx + v.y*m.yx + v.z*m.zx,
		v.x*m.xy + v.y*m.yy + v.z*m.zy,
		v.x*m.xz + v.y*m.yz + v.z*m.zz
	);
}

struct Space::vec3d operator*(const struct Space::m3x3d &m, const struct Space::vec3d &v)
{
	return Space::vec3d(
		m.xx*v.x + m.xy*v.y + m.xz*v.z,
		m.yx*v.x + m.yy*v.y + m.yz*v.z,
		m.zx*v.x + m.zy*v.y + m.zz*v.z
	);
}

struct Space::m3x3d operator*(const struct Space::m3x3d &m, const struct Space::m3x3d &n)
{
	return Space::m3x3d(
		m.xx*n.xx + m.xy*n.yx + m.xz*n.zx, m.xx*n.xy + m.xy*n.yy + m.xz*n.zy, m.xx*n.xz + m.xy*n.yz + m.xz*n.zz,
		m.yx*n.xx + m.yy*n.yx + m.yz*n.zx, m.yx*n.xy + m.yy*n.yy + m.yz*n.zy, m.yx*n.xz + m.yy*n.yz + m.yz*n.zz,
		m.zx*n.xx + m.zy*n.yx + m.zz*n.zx, m.zx*n.xy + m.zy*n.yy + m.zz*n.zy, m.zx*n.xz + m.zy*n.yz + m.zz*n.zz
	);
}

struct Space::m3x3d operator+(const struct Space::m3x3d &m, const struct Space::m3x3d &n)
{
	return Space::m3x3d(
		m.xx + n.xx, m.xy + n.xy, m.xz + n.xz,
		m.yx + n.yx, m.yy + n.yy, m.yz + n.yz,
		m.zx + n.zx, m.zy + n.zy, m.zz + n.zz
	);
}

struct Space::m3x3d operator-(const struct Space::m3x3d &m, const struct Space::m3x3d &n)
{
	return Space::m3x3d(
		m.xx - n.xx, m.xy - n.xy, m.xz - n.xz,
		m.yx - n.yx, m.yy - n.yy, m.yz - n.yz,
		m.zx - n.zx, m.zy - n.zy, m.zz - n.zz
	);
}

struct Space::m3x3d operator+(const struct Space::m3x3d &m)
{
	return Space::m3x3d(
		m.xx, m.xy, m.xz,
		m.yx, m.yy, m.yz,
		m.zx, m.zy, m.zz
	);
}

struct Space::m3x3d operator-(const struct Space::m3x3d &m)
{
	return Space::m3x3d(
		-m.xx, -m.xy, -m.xz,
		-m.yx, -m.yy, -m.yz,
		-m.zx, -m.zy, -m.zz
	);
}

struct Space::m3x3d Space::munit()
{
	return Space::m3x3d(
		1., 0., 0.,
		0., 1., 0.,
		0., 0., 1.
	);
}

struct Space::m3x3d Space::mrotX(double fi)
{
	double s = sin(fi);
	double c = cos(fi);

	return Space::m3x3d(
		1., 0., 0.,
		0., c, -s,
		0., s, c
	);
}

struct Space::m3x3d Space::mrotY(double fi)
{
	double s = sin(fi);
	double c = cos(fi);

	return Space::m3x3d(
		c, 0., s,
		0., 1., 0.,
		-s, 0., c
	);
}

struct Space::m3x3d Space::mrotZ(double fi)
{
	double s = sin(fi);
	double c = cos(fi);

	return Space::m3x3d(
		c, -s, 0.,
		s, c, 0.,
		0., 0., 1.
	);
}

struct Space::m3x3d Space::mrotV(double fi, const struct Space::vec3d &v)
{
	struct Space::vec3d n = v.normalized();

	double s = sin(fi);
	double c = cos(fi);

	return Space::m3x3d(
		c + (1. - c)*n.x*n.x, (1. - c)*n.x*n.y - s*n.z, (1. - c)*n.x*n.z + s*n.y,
		(1. - c)*n.y*n.x + s*n.z, c + (1. - c)*n.y*n.y, (1. - c)*n.y*n.z - s*n.x,
		(1. - c)*n.z*n.x - s*n.y, (1. - c)*n.z*n.y + s*n.x, c + (1. - c)*n.z*n.z
	);
}
