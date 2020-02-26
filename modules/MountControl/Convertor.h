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

struct traject {
	double T;
	std::pair<double, double> V;
	std::pair<double, double> a;
	std::pair<double, double> startpos;
	std::pair<double, double> endpos;
};

class Convertor {
public:
	//�����������/����������
	Convertor();
	~Convertor();

	//������������� ���������� ����������
	void InitMountSpec(const double &st, const double &a, const double &V, const double &ratio, const double &steps);
	//������ ������� ��������� � ITRF
	void SetTelPosITRF(const double &x, const double &y, const double &z);
	void SetTelPosITRF(double *pos);
	//������ ������� ��������� �� lat, lon, elev
	void SetTelPosLatLonElev(const double &Lat, const double &Lon, const double &Elev);
	void SetTelPosLatLonElev(double *pos);
	//������ ����������� C � �������� mbars �� �������������� ������
	void SetTelTP(const double &T, const double &P);
	
	//������ ����������� ���������� Az, Elev
	void SetAzElevPos(const double &Jd, const double &Az, const double &Elev);
	//������ ����������� ���������� Ra, Dec
	void SetRaDecPos(const double &Jd, const double &Ra, const double &Dec);
	//������ ������� ��������� ���������� � ����� mot1, mot2
	void SetMotorsPos(const double &mot1, const double &mot2);
	//������ ������� �������������� Az, Elev -> Alph, Bet
	void SetConvMatr(const std::vector<double>& M);

	//�������������� ���������
	std::vector<double> AzElev2AlphBet(const std::vector<double>& AzElev);

	//��������� ���������� ��� ������� ����������
	traject CalcTraject(const double &inAlph, const double &inBet, const double &outAlph, const double &outBet);

	//��������� �� �������� ��������� � ��������
	traject GoToRaDec(const double &Jd, const double &outRa, const double &outDec);
	traject GoToAzElev(const double &Jd, const double &outAz, const double &outElev);

	//�������� �� �����������
	int LimitsAzElev(const double &outAz, const double &outElev);

	//�������� ������� ��������� � ������ �������� ���������
	std::pair<double, double> GetAzElevPos();
	std::pair<double, double> GetRaDecPos();
	std::pair<double, double> GetAlphBetPos();
	std::pair<double, double> GetMotorsPos();
private:
	//CurrentPosition
	double cur_Az, cur_Elev;
	double cur_Ra, cur_Dec;
	double cur_Alph, cur_Bet;
	double cur_mot1, cur_mot2;
	double cur_Jd;

	//MountSpec
	double num;
	double a_max;
	double V_max;
	double Gear_ratio;
	double circle_motor;
	double circle_mount;
	double step;

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

	//Initialization Integrator
	void InitIntegrator();

	//ConvMatr
	std::vector<double> R;

	bool SetDateAndPolePos(const double &Jd);
	void Convert2AlphBet();
	void RaDec2AzElev();
	void AzElev2RaDec();
	void AzElevR2XYZ_ITRF(double* pos);
};

double Modulus(double x, double y);
std::vector<double> AzElev2XYZ(const std::vector<double> &AzElev);
std::vector<double> XYZ2AzElev(const std::vector<double> &xyz);