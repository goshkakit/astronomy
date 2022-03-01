#include "cearth.h"
#include "space_force.h"

Space::CEarth::CEarth(double jd)
{
	setTime(jd);
}

Space::CEarth::CEarth(const Time::long_jd &ljd)
{
	setTime(ljd);
}

Space::CEarth::~CEarth()
{

}

void Space::CEarth::setTime(double jd)
{
	setTime(Time::long_jd::fromJD(jd));
}

void Space::CEarth::setTime(const Time::long_jd &ljd)
{
	fljd = ljd;

	Force()->set_time(date_msk(), time_msk(), &fajd, &fjdelt, &fjt);

	Force()->iers_update_matrix(fjt, fori.v, fajd, fjdelt);
}
