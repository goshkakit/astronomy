#include "long_jd.h"

#include <math.h>

const double Time::long_jd_base::SCALE_DAYS = (1. / STEP_DATE_LL);
const double Time::long_jd_base::SCALE_SECONDS = (86400. / STEP_DATE_LL);
const double Time::long_jd_base::SCALE_MILLIS = (86400.e3 / STEP_DATE_LL);
const double Time::long_jd_base::SCALE_MICROS = (86400.e6 / STEP_DATE_LL);
const double Time::long_jd_base::SCALE_NANOS = (86400.e9 / STEP_DATE_LL);
const double Time::long_jd_base::BASE_DATE = ((double)long_jd_base::BASE_DATE_LL);
const double Time::long_jd_base::UNIX_SCALE = (STEP_DATE_LL / 86400.);

Time::long_t Time::long_jd_base::longTfromRJD(double _RJD)
{
	return (long_t)round(_RJD * STEP_DATE_LL);
}

Time::long_t Time::long_jd_base::longTfromJDT(double _JDN, double _jtime)
{
	double JDN_d = floor(_JDN);
	long_t JDN_i = (long_t)JDN_d;
	double JDN_f = _JDN - JDN_d;

	_jtime += JDN_f;

	double jt_d = floor(_jtime);
	long_t jt_i = (long_t)jt_d;
	double jt_f = _jtime - jt_d;

	long_t JDN_l = JDN_i + jt_i - BASE_DATE_LL;
	long_t JT_l = (long_t)round(jt_f * STEP_DATE_LL);

	JDN_l += (JT_l & STEP_DATE_LL) ? 1LL : 0LL;

	long_t _sig = (JDN_l & MASK_SIG);
	long_t _day = (JDN_l & MAX_DAY) << (BITS_TIME);
	long_t _time = (JT_l & MASK_TIME);

	return (_sig | _day | _time);
}

Time::long_t Time::long_jd_base::longTfromDTP(int year, int month, int day, int hour /*= 12*/, int minute /*= 0*/, double second /*= 0.*/, double local_offset_s /*= 0.*/)
{
	int a = (14 - month) / 12;
	int y = year + 4800 - a;
	int m = month + 12 * a - 3;

	int _JDN = day + (153 * m + 2) / 5 + 365 * y + y / 4 - y / 100 + y / 400 - 32045;
	double _jtime = (second / 86400.) + (minute / 1440.) + (hour / 24.) - 0.5 - (local_offset_s / 86400.);

	return longTfromJDT(_JDN, _jtime);
}

Time::long_t Time::long_jd_base::longTfromDTF(double _date, double _time, double local_offset_s /*= 0.*/)
{
	int date_i = (int)floor(_date);

	int year = date_i / 10000;
	date_i %= 10000;
	int month = date_i / 100;
	date_i %= 100;
	int day = date_i;

	double time_d = floor(_time);
	double time_f = _time - time_d;
	int time_i = (int)time_d;

	int hour = time_i / 10000;
	time_i %= 10000;
	int minute = time_i / 100;
	time_i %= 100;
	double second = time_f + time_i;

	return longTfromDTP(year, month, day, hour, minute, second, local_offset_s);
}

Time::long_t Time::long_jd_base::longTfromUXT(long_t _UXT)
{
	long_t UXT_l = _UXT + UNIX_SHIFT;
	double UXT_date = (double)(UXT_l / 86400) + BASE_DATE;
	double UXT_time = (double)(UXT_l % 86400) / 86400.;

	return longTfromJDT(UXT_date, UXT_time);
}

void Time::long_jd::date(int *year, int *month, int *day, double local_offset_s /*= 0.*/) const
{
	long_t lt_l = (lt + (long_t)round(UNIX_SCALE*local_offset_s + 0.5*STEP_DATE_LL));
	int _JDN = (int)((lt_l >> BITS_TIME) + BASE_DATE_LL);

	int a = _JDN + 32044;
	int b = (4 * a + 3) / 146097;
	int c = a - ((146097 * b) / 4);
	int d = (4 * c + 3) / 1461;
	int e = c - ((1461 * d) / 4);
	int m = (5 * e + 2) / 153;
	int n = m / 10;

	if (year) *year = 100 * b + d - 4800 + n;
	if (month) *month = m + 3 - 12 * n;
	if (day) *day = e - ((153 * m + 2) / 5) + 1;
}

void Time::long_jd::time(int *hour, int *minute, double *second, double local_offset_s /*= 0.*/) const
{
	long_t lt_l = (lt + (long_t)round(UNIX_SCALE*local_offset_s + 0.5*STEP_DATE_LL));

	double _time = SCALE_SECONDS * (lt_l & MASK_TIME);
	double time_d = floor(_time);
	double time_f = _time - time_d;
	int time_i = (int)time_d;

	if (hour) *hour = time_i / 3600;
	time_i %= 3600;
	if (minute) *minute = time_i / 60;
	time_i %= 60;
	if (second) *second = time_f + time_i;
}

double Time::long_jd::date(double local_offset_s /*= 0.*/) const
{
	int year, month, day;
	date(&year, &month, &day, local_offset_s);
	return ((double)day + 100.*month + 10000.*year);
}

double Time::long_jd::time(double local_offset_s /*= 0.*/) const
{
	int hour, minute;
	double second;
	time(&hour, &minute, &second, local_offset_s);
	return (second + 100.*minute + 10000.*hour);
}

Time::long_t Time::long_jd::unix_time() const
{
	return ((long_t)round(SCALE_SECONDS * lt) - UNIX_SHIFT);
}

Time::long_t Time::long_jd::nanos() const
{
    long_t l_days = lt >> BITS_TIME;
    long_t l_time = lt & MASK_TIME;

    double days_ns = 86400.e9 * l_days;
    double time_ns = SCALE_NANOS * l_time;

    double days_ns_d = floor(days_ns);
    long_t days_ns_l = (long_t)days_ns_d;
    double days_ns_f = days_ns - days_ns_d;

    long_t time_ns_l = (long_t)round(time_ns + days_ns_f);

    return (days_ns_l + time_ns_l);
}

Time::long_t Time::long_jd::micros() const
{
    long_t l_days = lt >> BITS_TIME;
    long_t l_time = lt & MASK_TIME;

    double days_us = 86400.e6 * l_days;
    double time_us = SCALE_MICROS * l_time;

    double days_us_d = floor(days_us);
    long_t days_us_l = (long_t)days_us_d;
    double days_us_f = days_us - days_us_d;

    long_t time_us_l = (long_t)round(time_us + days_us_f);

    return (days_us_l + time_us_l);
}

Time::long_jd & Time::long_jd::addNanos(long_t _nanos)
{
    double l_date = (double)(_nanos / 86400000000000LL) + BASE_DATE;
    double l_time = (double)(_nanos % 86400000000000LL) / 86400.e9;

    lt += longTfromJDT(l_date, l_time);
    return *this;
}

Time::long_jd Time::long_jd::operator*(double m) const
{
	return long_jd((long_t)round(lt * m));
}

Time::long_jd Time::long_jd::operator/(double m) const
{
	return long_jd((long_t)round(lt / m));
}

Time::long_jd & Time::long_jd::operator*=(double m)
{
	lt = (long_t)round(lt * m); return *this;
}

Time::long_jd & Time::long_jd::operator/=(double m)
{
	lt = (long_t)round(lt / m); return *this;
}

Time::long_jd operator*(long long m, const Time::long_jd &ljd)
{
	return Time::long_jd(m * ljd.lt);
}

Time::long_jd operator*(double m, const Time::long_jd &ljd)
{
	return Time::long_jd((Time::long_t)round(m * ljd.lt));
}

const Tests::Time::tReferenceValues Tests::Time::ref1 = {
	1542119304LL,
	2458436,
	2458436.1030555555555555555556,
	58436.103055555555555555555556,
	58435.603055555555555555555556,
	.10305555555555555555555555556,
	20181113., 142824.,
	20181113., 172824.,
	20181113.,  62824.,
};

const Tests::Time::tReferenceValues Tests::Time::ref2 = {
	1470531493LL,
	2457607,
	2457607.5404282407407407407407,
	57607.540428240740740740740741,
	57607.040428240740740740740741,
	.54042824074074074074074074074,
	20160807.,   5813.,
	20160807.,  35813.,
	20160806., 165813.,
};

const Tests::Time::tReferenceValues Tests::Time::ref3 = {
	7817806798LL,
	2531071,
	2531071.3749768518518518518519,
	131071.37497685185185185185185,
	131070.87497685185185185185185,
	.37497685185185185185185185185,
	22170926., 205958.,
	22170926., 235958.,
	22170926., 125958.,
};

#include <stdio.h>

template<class T> static T _abs_t(T v)
{
	return (v < 0) ? (-v) : (v);
}

bool Tests::Time::compareToRef(const ::Time::long_jd &val, const struct tReferenceValues &ref, double thr)
{
	long long ll;
	int i;
	double d, md = 0.;

	printf("Test: long_jd::unix_time()\n");
	if ((ll = _abs_t(val.unix_time() - ref.unix_time)) > 1LL)
	{
		printf("FAILED: %lld > %lld\n", ll, 1LL);
		return false;
	}

	printf("Test: long_jd::JDN()\n");
	if ((i = _abs_t(val.JDN() - ref.JDN)) > 0)
	{
		printf("FAILED: %d > %d\n", i, 0);
		return false;
	}

	printf("Test: long_jd::JD()\n");
	if ((d = (fabs(val.JD() - ref.JD))*86400) > thr)
	{
		printf("FAILED: %e > %e\n", d, thr);
		return false;
	}
	md = fmax(md, d);

	printf("Test: long_jd::RJD()\n");
	if ((d = (fabs(val.RJD() - ref.RJD))*86400) > thr)
	{
		printf("FAILED: %e > %e\n", d, thr);
		return false;
	}
	md = fmax(md, d);

	printf("Test: long_jd::MJD()\n");
	if ((d = (fabs(val.MJD() - ref.MJD))*86400) > thr)
	{
		printf("FAILED: %e > %e\n", d, thr);
		return false;
	}
	md = fmax(md, d);

	printf("Test: long_jd::jtime()\n");
	if ((d = (fabs(val.jtime() - ref.jtime))*86400) > thr)
	{
		printf("FAILED: %e > %e\n", d, thr);
		return false;
	}
	md = fmax(md, d);

	printf("Test: long_jd::date(UTC)\n");
	if ((d = (fabs(val.date() - ref.date_utc))*86400) > thr)
	{
		printf("FAILED: %e > %e\n", d, thr);
		return false;
	}
	md = fmax(md, d);

	printf("Test: long_jd::time(UTC)\n");
	if ((d = fabs(val.time() - ref.time_utc)) > thr)
	{
		printf("FAILED: %e > %e\n", d, thr);
		return false;
	}
	md = fmax(md, d);

	printf("Test: long_jd::date(GMT+3)\n");
	if ((d = (fabs(val.date(3. * 3600.) - ref.date_msk))*86400) > thr)
	{
		printf("FAILED: %e > %e\n", d, thr);
		return false;
	}
	md = fmax(md, d);

	printf("Test: long_jd::time(GMT+3)\n");
	if ((d = fabs(val.time(3. * 3600.) - ref.time_msk)) > thr)
	{
		printf("FAILED: %e > %e\n", d, thr);
		return false;
	}
	md = fmax(md, d);

	printf("Test: long_jd::date(GMT-8)\n");
	if ((d = (fabs(val.date(-8. * 3600.) - ref.date_sfr))*86400) > thr)
	{
		printf("FAILED: %e > %e\n", d, thr);
		return false;
	}
	md = fmax(md, d);

	printf("Test: long_jd::time(GMT-8)\n");
	if ((d = fabs(val.time(-8. * 3600.) - ref.time_sfr)) > thr)
	{
		printf("FAILED: %e > %e\n", d, thr);
		return false;
	}
	md = fmax(md, d);

	printf("Max delta: %.2e s\n", md);
	return true;
}

bool Tests::Time::testLongJD(const struct tReferenceValues &ref)
{
	printf("\nTest: long_jd::fromUnixTime()\n");
	if (!compareToRef(::Time::long_jd::fromUnixTime(ref.unix_time), ref, 1e-9))
		return false;

	printf("\nTest: long_jd::fromJD()\n");
	if (!compareToRef(::Time::long_jd::fromJD(ref.JD), ref, 1e-4))
		return false;

	printf("\nTest: long_jd::fromRJD()\n");
	if (!compareToRef(::Time::long_jd::fromRJD(ref.RJD), ref, 2e-6))
		return false;

	printf("\nTest: long_jd::fromMJD()\n");
	if (!compareToRef(::Time::long_jd::fromMJD(ref.MJD), ref, 2e-6))
		return false;

	printf("\nTest: long_jd::fromJDT()\n");
	if (!compareToRef(::Time::long_jd::fromJDT(ref.JDN, ref.jtime), ref, 1e-9))
		return false;

	printf("\nTest: long_jd::fromDTF(UTC)\n");
	if (!compareToRef(::Time::long_jd::fromDTF(ref.date_utc, ref.time_utc), ref, 1e-9))
		return false;

	printf("\nTest: long_jd::fromDTF(GMT+3)\n");
	if (!compareToRef(::Time::long_jd::fromDTF(ref.date_msk, ref.time_msk, 10800), ref, 1e-9))
		return false;

	printf("\nTest: long_jd::fromDTF(GMT-8)\n");
	if (!compareToRef(::Time::long_jd::fromDTF(ref.date_sfr, ref.time_sfr, -28800), ref, 1e-9))
		return false;

	return true;
}

bool Tests::Time::testLongJDops(const struct tReferenceValues &ref)
{
    const ::Time::long_jd t0 = ::Time::long_jd::fromUnixTime(ref.unix_time);

    ::Time::long_jd t1 = t0;

    printf("Test: long_jd::operator==()\n");
    if (!(t0 == t1))
    {
        printf("FAILED\n");
        return false;
    }

    printf("Test: long_jd::operator<=()\n");
    if (!(t0 <= t1))
    {
        printf("FAILED\n");
        return false;
    }

    printf("Test: long_jd::operator>=()\n");
    if (!(t0 >= t1))
    {
        printf("FAILED\n");
        return false;
    }

    t1 = ::Time::long_jd::fromUnixTime(ref.unix_time - 2);

    printf("Test: long_jd::operator!=()\n");
    if (!(t0 != t1))
    {
        printf("FAILED\n");
        return false;
    }

    printf("Test: long_jd::operator<=()\n");
    if (!(t1 <= t0))
    {
        printf("FAILED\n");
        return false;
    }

    printf("Test: long_jd::operator>=()\n");
    if (!(t0 >= t1))
    {
        printf("FAILED\n");
        return false;
    }

    printf("Test: long_jd::operator<()\n");
    if (!(t1 < t0))
    {
        printf("FAILED\n");
        return false;
    }

    printf("Test: long_jd::operator>()\n");
    if (!(t0 > t1))
    {
        printf("FAILED\n");
        return false;
    }

    ::Time::long_jd t2;
    t1 = +t0;
    t2 = -t0;

    printf("Test: long_jd::operator+()\n");
    if (!(t0 == t1))
    {
        printf("FAILED\n");
        return false;
    }

    printf("Test: long_jd::operator-()\n");
    if (!(t2 != t1))
    {
        printf("FAILED\n");
        return false;
    }
    if (!(t2 == -t1))
    {
        printf("FAILED\n");
        return false;
    }
    if (!((+::Time::long_jd(-100LL)) == (-::Time::long_jd(100LL))))
    {
        printf("FAILED\n");
        return false;
    }

    printf("Test: long_jd::operator-(const long_jd &)\n");
    t1 = ::Time::long_jd(2LL) - ::Time::long_jd(3LL);
    if (!(t1 == ::Time::long_jd(-1LL)))
    {
        printf("FAILED\n");
        return false;
    }

    printf("Test: long_jd::operator+=(const long_jd &)\n");
    t1 += ::Time::long_jd(2LL);
    if (!(t1 == ::Time::long_jd(1LL)))
    {
        printf("FAILED\n");
        return false;
    }

    printf("Test: long_jd::operator+(const long_jd &)\n");
    t2 = ::Time::long_jd(-2LL) + ::Time::long_jd(3LL);
    if (!(t2 == ::Time::long_jd(1LL)))
    {
        printf("FAILED\n");
        return false;
    }

    printf("Test: long_jd::operator-=(const long_jd &)\n");
    t2 -= ::Time::long_jd(2LL);
    if (!(t2 == ::Time::long_jd(-1LL)))
    {
        printf("FAILED\n");
        return false;
    }

    printf("Test: long_jd::operator*(long long)\n");
    t1 = t1 * 4LL;
    if (!(t1 == ::Time::long_jd(4LL)))
    {
        printf("FAILED\n");
        return false;
    }

    printf("Test: operator*(long long, const long_jd &)\n");
    t1 = 3LL * t1;
    if (!(t1 == ::Time::long_jd(12LL)))
    {
        printf("FAILED\n");
        return false;
    }

    printf("Test: long_jd::operator*(double)\n");
    t1 = t1 * 0.5;
    if (!(t1 == ::Time::long_jd(6LL)))
    {
        printf("FAILED\n");
        return false;
    }

    printf("Test: operator*(double, const long_jd &)\n");
    t1 = 0.5 * t1;
    if (!(t1 == ::Time::long_jd(3LL)))
    {
        printf("FAILED\n");
        return false;
    }

    printf("Test: long_jd::operator*=(long long)\n");
    t2 *= 3LL;
    if (!(t2 == ::Time::long_jd(-3LL)))
    {
        printf("FAILED\n");
        return false;
    }

    printf("Test: long_jd::operator/=(double)\n");
    t2 /= -0.25;
    if (!(t2 == ::Time::long_jd(12LL)))
    {
        printf("FAILED\n");
        return false;
    }

    printf("Test: long_jd::operator/=(long long)\n");
    t2 /= -2LL;
    if (!(t2 == ::Time::long_jd(-6LL)))
    {
        printf("FAILED\n");
        return false;
    }

    printf("Test: long_jd::operator/(double)\n");
    t2 = t2 / 3.0;
    if (!(t2 == ::Time::long_jd(-2LL)))
    {
        printf("FAILED\n");
        return false;
    }

    printf("Test: long_jd::operator*=(double)\n");
    t2 *= 7.0;
    if (!(t2 == ::Time::long_jd(-14LL)))
    {
        printf("FAILED\n");
        return false;
    }

    printf("Test: long_jd::operator/(long long)\n");
    t2 = t2 / -7LL;
    if (!(t2 == ::Time::long_jd(2LL)))
    {
        printf("FAILED\n");
        return false;
    }

    return true;
}

bool Tests::Time::testLongJD()
{
	printf("\nTest long_jd conversions with reference value 1.\n");
	if (!testLongJD(ref1))
		return false;

	printf("\nTest long_jd conversions with reference value 2.\n");
	if (!testLongJD(ref2))
		return false;

	printf("\nTest long_jd conversions with reference value 3.\n");
	if (!testLongJD(ref3))
		return false;

	printf("\nTest long_jd operations with reference value 1.\n");
	if (!testLongJDops(ref1))
		return false;

	printf("\nTest long_jd operations with reference value 2.\n");
	if (!testLongJDops(ref2))
		return false;

	printf("\nTest long_jd operations with reference value 3.\n");
	if (!testLongJDops(ref3))
		return false;

	printf("\nOK\n");
	return true;
}
