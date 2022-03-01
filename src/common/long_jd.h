#pragma once

#ifndef _LONG_JD_H
#define _LONG_JD_H

namespace Time
{
	typedef long long long_t;

	class long_jd_base
	{
	public:
		static const size_t BYTES = (sizeof(long_t));
		static const size_t BITS = (BYTES * 8);
		static const size_t BITS_SIG = (1);
		static const size_t BITS_TIME = (46);
		static const size_t BITS_DATE = (BITS - BITS_SIG - BITS_TIME);
		static const long_t STEP_DATE_LL = (1LL << BITS_TIME); // 2^46
		static const double SCALE_DAYS; // day (1/2^46 ~ 1.421e-14)
		static const double SCALE_SECONDS; // s (86400/2^46 ~ 1.228e-9)
		static const double SCALE_MILLIS; // ms (86400*10^3/2^46 ~ 1.228e-6)
		static const double SCALE_MICROS; // us (86400*10^6/2^46 ~ 1.228e-3)
		static const double SCALE_NANOS; // ns (86400*10^9/2^46 ~ 1.228)
		static const long_t CAP_DAYS = (1LL << BITS_DATE); // 2^17
		static const long_t MAX_DAY = (CAP_DAYS - 1); // 2^17 - 1
		static const long_t MASK_TIME = (STEP_DATE_LL - 1); // 2^46 - 1
		static const long_t MASK_DATE = (MAX_DAY << BITS_TIME);
		static const long_t MASK_SIG = (1LL << (BITS - 1)); // 2^63
		static const long_t BASE_DATE_LL = (2400000LL);
		static const double BASE_DATE; // 2400000.0
		static const long_t UNIX_SHIFT = ((4881175LL - (BASE_DATE_LL * 2)) * 43200LL); // (2440587.5 - 2400000) * 86400
		static const double UNIX_SCALE; // (2^46/86400. ~ 8.145e+8)

		// initial accuracy ~ 1 us
		static long_t longTfromRJD(double _RJD);
		// initial accuracy ~ 0.5 ns
		static long_t longTfromJDT(double _JDN, double _jtime);
		// initial accuracy ~ 0.5 ns (use longTfromJDT() internally)
		static long_t longTfromDTP(int year, int month, int day, int hour = 12, int minute = 0, double second = 0., double local_offset_s = 0.);
		// initial accuracy ~ 0.5 ns (use longTfromJDT() internally)
		static long_t longTfromDTF(double _date, double _time, double local_offset_s = 0.);
		// initial accuracy ~ 0.5 ns (use longTfromJDT() internally)
		static long_t longTfromUXT(long long _UXT);
	};

	class long_jd : public long_jd_base
	{
	public:
		// initial accuracy < 0.1 ms
		static long_jd fromJD(double _JD) {
			return long_jd(longTfromRJD(_JD - BASE_DATE));
		}
		// initial accuracy ~ 1 us
		static long_jd fromRJD(double _RJD) {
			return long_jd(longTfromRJD(_RJD));
		}
		// initial accuracy ~ 1 us
		static long_jd fromMJD(double _MJD) {
			return long_jd(longTfromRJD(_MJD + 0.5));
		}
		// initial accuracy ~ 0.5 ns
		static long_jd fromJDT(double _JDN, double _jtime) {
			return long_jd(longTfromJDT(_JDN, _jtime));
		}
		// initial accuracy ~ 0.5 ns
		static long_jd fromDTP(int year, int month, int day, int hour = 12, int minute = 0, double second = 0., double local_offset_s = 0.) {
			return long_jd(longTfromDTP(year, month, day, hour, minute, second, local_offset_s));
		}
		// initial accuracy ~ 0.5 ns
		static long_jd fromDTF(double _date, double _time, double local_offset_s = 0.) {
			return long_jd(longTfromDTF(_date, _time, local_offset_s));
		}
		// initial accuracy ~ 0.5 ns
		static long_jd fromUnixTime(long long _UXT) {
			return long_jd(longTfromUXT(_UXT));
		}

		double JD() const {
			return (SCALE_DAYS * lt + BASE_DATE);
		}
		double RJD() const {
			return (SCALE_DAYS * lt);
		}
		double MJD() const {
			return (SCALE_DAYS * lt - 0.5);
		}
		int JDN() const {
			return (int)((lt >> BITS_TIME) + BASE_DATE_LL);
		}
		double jtime() const {
			return (SCALE_DAYS * (lt & MASK_TIME));
		}
		long_t nanos() const;
		long_t micros() const;
		double millis() const {
			return (SCALE_MILLIS * lt);
		}
		double seconds() const {
			return (SCALE_SECONDS * lt);
		}
		double days() const {
			return (SCALE_DAYS * lt);
		}
		void date(int *year, int *month, int *day, double local_offset_s = 0.) const;
		void time(int *hour, int *minute, double *second, double local_offset_s = 0.) const;
		double date(double local_offset_s = 0.) const;
		double time(double local_offset_s = 0.) const;
		long_t unix_time() const;

		long_jd & addNanos(long long _nanos);
		long_jd & addNanos(double _nanos) {
			lt += longTfromRJD(_nanos / 86400.e9); return *this;
		}
		long_jd & addMicros(double _micros) {
			lt += longTfromRJD(_micros / 86400.e6); return *this;
		}
		long_jd & addMillis(double _millis) {
			lt += longTfromRJD(_millis / 86400.e3); return *this;
		}
		long_jd & addSeconds(double _seconds) {
			lt += longTfromRJD(_seconds / 86400.); return *this;
		}
		long_jd & addDays(double _days) {
			lt += longTfromRJD(_days); return *this;
		}
		long_jd & addDays(int _days) {
			lt += ((long_t)_days) << BITS_TIME; return *this;
		}

		bool operator==(const long_jd &ljd) const {
			return (lt == ljd.lt);
		}
		bool operator!=(const long_jd &ljd) const {
			return (lt != ljd.lt);
		}
		bool operator>=(const long_jd &ljd) const {
			return (lt >= ljd.lt);
		}
		bool operator<=(const long_jd &ljd) const {
			return (lt <= ljd.lt);
		}
		bool operator>(const long_jd &ljd) const {
			return (lt > ljd.lt);
		}
		bool operator<(const long_jd &ljd) const {
			return (lt < ljd.lt);
		}

		long_jd operator+(const long_jd &_ljd) const {
			return long_jd(lt + _ljd.lt);
		}
		long_jd operator-(const long_jd &_ljd) const {
			return long_jd(lt - _ljd.lt);
		}
		long_jd operator+() const {
			return long_jd(+lt);
		}
		long_jd operator-() const {
			return long_jd(-lt);
		}
		long_jd operator*(long long m) const {
			return long_jd(lt * m);
		}
		long_jd operator/(long long m) const {
			return long_jd(lt / m);
		}
		long_jd operator*(double m) const;
		long_jd operator/(double m) const;

		long_jd & operator+=(const long_jd &_ljd) {
			lt += _ljd.lt; return *this;
		}
		long_jd & operator-=(const long_jd &_ljd) {
			lt -= _ljd.lt; return *this;
		}
		long_jd & operator*=(long long m) {
			lt *= m; return *this;
		}
		long_jd & operator/=(long long m) {
			lt /= m; return *this;
		}
		long_jd & operator*=(double m);
		long_jd & operator/=(double m);

		long_jd() : lt(0LL) {}
		long_jd(const long_jd &_ljd) : lt(_ljd.lt) {}
		long_jd & operator=(const long_jd &_ljd) {
			lt = _ljd.lt; return *this;
		}
		explicit long_jd(long_t _lt) : lt(_lt) {}
		explicit long_jd(double _JD) : lt(longTfromRJD(_JD - BASE_DATE)) {}
		long_jd(double _JDN, double _jtime) : lt(longTfromJDT(_JDN, _jtime)) {}

		long_t lt;
	};
}

Time::long_jd operator*(long long m, const Time::long_jd &ljd);
Time::long_jd operator*(double m, const Time::long_jd &ljd);

namespace Tests
{
	namespace Time
	{
		struct tReferenceValues
		{
			long long unix_time;
			int JDN;
			double JD, RJD, MJD, jtime, date_utc, time_utc, date_msk, time_msk, date_sfr, time_sfr;
		};

		extern const struct tReferenceValues ref1, ref2, ref3;

		bool compareToRef(const ::Time::long_jd &val, const struct tReferenceValues &ref, double thr);

		bool testLongJD(const struct tReferenceValues &ref);

		bool testLongJDops(const struct tReferenceValues &ref);

		bool testLongJDmods(const struct tReferenceValues &ref);

		bool testLongJD();
	}
}

#endif//_LONG_JD_H
