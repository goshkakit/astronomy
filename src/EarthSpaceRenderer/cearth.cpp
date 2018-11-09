#include "cearth.h"
#include "space_force.h"

#include "common/DataConverter.h"

Space::CEarth::CEarth(double jd)
{
	setTime(jd);
}

Space::CEarth::~CEarth()
{

}

void Space::CEarth::setTime(double jd)
{
	DataConverter dc;

	fjd = jd;
	fdate_utc = dc.JDtoYYYYMMDD(fjd);
	ftime_utc = dc.SECtoHHMMSS(fdate_utc, fjd);

	double fjd_msk = jd + 0.125;
	fdate_msk = dc.JDtoYYYYMMDD(fjd_msk);
	ftime_msk = dc.SECtoHHMMSS(fdate_msk, fjd_msk);

	Force->set_time(fdate_msk, ftime_msk, &fajd, &fjdelt, &fjt);

	Force->iers_update_matrix(fjt, fori.v, fajd, fjdelt);
}
