#pragma once

#include <map>

#include "WGS84.h"
#include "InfluenceForce\InfluenceForce.h"
#include "OrbitIntegration\IPredictOrbitMod.h"
#include "OrbitIntegration\PredictOrbitMod.h"
#include "common\DataConverter.h"

extern "C" {
#include "novac/novas.h"
#include "novac/eph_manager.h"
}

class Convertor {
public:
	Convertor();
	~Convertor();

	void SetTelPosITRF(const double &x, const double &y, const double &z);
	void SetTelPosITRF(double *pos);
	void SetTelPosLatLonElev(const double &Lat, const double &Lon, const double &Elev);
	void SetTelPosLatLonElev(double *pos);
	void SetTelTP(const double &T, const double &P);
	void SetAzElevPos(const double &Jd, const double &Az, const double &Elev);
	void SetRaDecPos(const double &Jd, const double &Ra, const double &Dec);
	void SetAlphBetPos(const double &Jd, const double &Ra, const double &Dec);
	void SetConvMatr(const std::vector<double>& M);

	std::pair<double, double> GetAzElevPos();
	std::pair<double, double> GetRaDecPos();
	std::pair<double, double> GetAlphBetPos();
private:
	//CurrentPosition
	double cur_Az, cur_Elev;
	double cur_Ra, cur_Dec;
	double cur_Alph, cur_Bet;
	double cur_Jd;

	//Date
	double dAT;
	double dUT1;
	double jd_tt;
	double jd_UT1;
	double deltaT;
	double mjd0 = 2400000.50;

	//Pole
	double xp, yp;

	//TelescopePosition
	double tel_pos_ITRF[3];
	on_surface tel;

	//Integrator
	Force::InfluenceForce *IForce1;
	Orbit::PredictOrbitSat POSat1;
	DataConverter DC1;

	//Initialization
	void InitIntegrator();

	//ConvMatr
	std::vector<double> A;

	bool SetDateAndPolePos(const double &Jd);
	void Convert2AlphBet();
	void RaDec2AzElev();
	void AzElev2RaDec();
	void AzElevR2XYZ_ITRF(double* pos);
};

double Modulus(double x, double y);