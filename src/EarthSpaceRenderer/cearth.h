#pragma once

#include "space_defines.h"
#include "space_math.h"

namespace Space
{
	class CEarth
	{
	public:
		CEarth(double jd = 2400000.);
		virtual ~CEarth();

		void setTime(double jd);

		double jd() const {
			return fjd;
		}
		double date_utc() const {
			return fdate_utc;
		}
		double time_utc() const {
			return ftime_utc;
		}
		double date_msk() const {
			return fdate_msk;
		}
		double time_msk() const {
			return ftime_msk;
		}
		double ajd() const {
			return fajd;
		}
		double jdelt() const {
			return fjdelt;
		}
		double jt() const {
			return fjt;
		}

		const struct m3x3d & orientation() const {
			return fori;
		}
		const struct m3x3d rotator() const {
			return fori.transposed();
		}

	private:
		double fjd, fdate_utc, ftime_utc;
		double fdate_msk, ftime_msk;
		double fajd, fjdelt, fjt;
		m3x3d fori;
	};
}
