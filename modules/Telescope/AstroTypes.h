#pragma once

#include <string>

class TelescopeViewParams
{
public:

	TelescopeViewParams() {}

	double az, el;
	double range;

	double TimeAfterEpohe;
};

class TelescopeTrack
{
public:
	TelescopeTrack() {};

	double stepSec_;
	std::vector<TelescopeViewParams> track;
};

class TelescopObject
{
public:
	TelescopObject() {}

	int id;						// Observer number
	double lon, lat, height;	// Telescop  position
								//double x,y,z;
	double rate;				// Frame rate
	double FOVwa;				// FOV 05 width on ASC
	double FOVwb;				// FOV 05 width on Dec
	double sig;					// FOR TRACK SIMULATION ONLY :sigma of Ra,Dec measurements errors [ang.sec]
	double max_dur;				// FOR TRACK SIMULATION ONLY :maximum track duration [sec]
	std::string name;
};

class AstroConstant
{
public:
	AstroConstant() {}

	const double RG = 180.0 / M_PI;	// Radian-degree convertion
	const double GR = 1.0 / RG;

	const double DtoSec = 3600.0;
};

class TelescopeDataTime
{
public:
	int Y, M, D;
	int hh, mm, ss;

	// if set
	double ms = 0;
	double JD = 0;

	double getYYYYMMDD()
	{
		double YYYYMMDD = 0;
		YYYYMMDD = 10000.0 * Y + 100.0 * M + D;
		return YYYYMMDD;
	}

	double getHHMMSS()
	{
		double HHMMSS = 0;
		HHMMSS = 10000.0 * hh + 100.0 * mm + ss;
		return HHMMSS;
	}

	double getSecOfDay()
	{
		double sec = ss + 60.0*mm + 60.0*60.0*hh;
		return sec;
	}
};