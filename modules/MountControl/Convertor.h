#pragma once

#include <map>
#include <string>

#include "WGS84.h"
#include "InfluenceForce\InfluenceForce.h"
#include "OrbitIntegration\IPredictOrbitMod.h"
#include "OrbitIntegration\PredictOrbitMod.h"
#include "common\DataConverter.h"
#include "common\TLELoader.h"

extern "C" {
#include "novac/novas.h"
#include "novac/eph_manager.h"
}

//Траектория переезда
struct traject {
	double T;
	std::pair<double, double> V;
	std::pair<double, double> a;
	std::pair<double, double> startpos;
	std::pair<double, double> endpos;
};

//Структура меток для сопровождения
struct AccompPoints {
	double V0;
	double T;
	double N;
};

//Структура угловых измерений
struct Angs {
	double Jd;
	double ang1;
	double ang2;
};

//Структура положения двигателей
struct MotPos {
	double mot1;
	double mot2;
};

//Структура угловых измерений и угловых скоростей
struct AngsVels :Angs {
	double ang1_vel;
	double ang2_vel;
};

//Управляющий класс
class Convertor {
public:
	//Конструктор/Деструктор
	Convertor();
	~Convertor();

	//Инициализация параметров монтировки
	void InitMountSpec(const double &st, const double &a, const double &V, const double &ratio, const double &steps);
	//Задать позицию телескопа в ITRF
	void SetTelPosITRF(const double &x, const double &y, const double &z);
	void SetTelPosITRF(double *pos);
	//Задать позицию телескопа по lat, lon, elev
	void SetTelPosLatLonElev(const double &Lat, const double &Lon, const double &Elev);
	void SetTelPosLatLonElev(double *pos);
	//Задать температуру C и давление mbars на наблюдательном пункте
	void SetTelTP(const double &T, const double &P);
	
	//Задать направление наблюдения Az, Elev
	void SetAzElevPos(const Angs &AzElev);
	//Задать направление наблюдения Ra, Dec
	void SetRaDecPos(const Angs &RaDec);
	//Задать текущее положение двигателей в шагах mot1, mot2
	void SetMotorsPos(const MotPos &Motors);
	//Задать матрицу преобразования Az, Elev -> Alph, Bet
	void SetConvMatr(const std::vector<double>& M);

	//Преобразования координат
	Angs AzElev2AlphBet(const Angs &AzElev);
	Angs AlphBet2AzElev(const Angs &AlphBet);
	MotPos AlphBet2Motors(const Angs &AlphBet);
	MotPos AzElev2Motors(const Angs &AzElev);

	//Вычислить траекторию для шаговых двигателей
	traject CalcTraject(const Angs &inAlphBet, const Angs &outAlphBet);

	//Переехать из текущего положения в заданное
	traject GoToRaDec(const Angs &outRaDec);
	traject GoToAzElev(const Angs &outAzElev);

	//Проверка на ограничения
	int LimitsAzElev(const Angs &outAzElev);

	//Задать путь к TLE файлу
	void SetTLEparams(const std::string &path, const int &num_str);
	//Задать интервал между точками в сек
	void SetStep(const double &step);
	//Вычислить массив углов для сопровождения
	std::vector<Angs> CalculateAzElev(const int &NoradId, const double &JdStart, const double &JdEnd);
	std::vector<AngsVels> CalculateAlphBetVels(const int &NoradId, const double &JdStart, const double &JdEnd);
	//Вычисление массива точек для микроконтроллера
	std::pair<std::vector<AccompPoints>, std::vector<AccompPoints>> CalculatePoints(const std::vector<AngsVels> &AzElevArr);

	//Получить текущее положение в разных системах координат
	Angs GetAzElevPos();
	Angs GetRaDecPos();
	Angs GetAlphBetPos();
	MotPos GetMotorsPos();
private:
	//CurrentPosition
	Angs curAzElev;
	Angs curRaDec;
	Angs curAlphBet;
	MotPos curMotors;

	//TelescopePosition
	double tel_pos_ITRF[3];
	on_surface tel;

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

	//Integrator
	Force::InfluenceForce *IForce1;
	Orbit::PredictOrbitSat POSat1;
	DataConverter DC1;

	//Initialization Integrator
	void InitIntegrator();

	//TLEparams
	std::string TLEpath;
	int numb_of_line;
	TLELoader tleLoader;
	cOrbit* orb;
	double stepSec;
	double stepMin;

	//ConvMatr
	std::vector<double> R;

	//Функции внутренней конвертации
	bool SetDateAndPolePos(const double &Jd);
	void Convert2AlphBet();
	void RaDec2AzElev();
	void AzElev2RaDec();
	void AzElevR2XYZ_ITRF(double* pos);

	//Выбрать орбиту NaradId
	bool GetOrb(const int &NoradId);
};

double Modulus(double x, double y);
std::vector<double> AzElev2XYZ(const Angs &AzElev);
Angs XYZ2AzElev(const std::vector<double> &xyz);