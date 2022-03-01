#pragma once

#include "space_defines.h"
#include "space_math.h"
#include "common/long_jd.h"

namespace Space
{
	class CEarth
	{
	public:
		CEarth(double jd = 2400000.);
		CEarth(const Time::long_jd &ljd);
		virtual ~CEarth();

		void setTime(double jd);
		void setTime(const Time::long_jd &ljd);

		const Time::long_jd & ljd() const {
			return fljd;
		}
		double jd() const {
			return fljd.JD();
		}
		double date_utc() const {
			return fljd.date();
		}
		double time_utc() const {
			return fljd.time();
		}
		double date_msk() const {
			return fljd.date(10800.);
		}
		double time_msk() const {
			return fljd.time(10800.);
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
		Time::long_jd fljd;
		double fajd, fjdelt, fjt;
		m3x3d fori;
	};
}
