#include "Convertor.h"
#include "WGS84.h"
#include "Mat3x3.h"

#include <iostream>
#include <math.h>

//Конструктор
Convertor::Convertor() {
	//CurrentPosition
	cur_Jd = 0;
	cur_Az = 0;
	cur_Elev = 0;
	cur_Ra = 0;
	cur_Dec = 0;
	cur_Alph = 0;
	cur_Bet = 0;
	cur_mot1 = 0;
	cur_mot2 = 0;

	//TelescopePosition
	tel_pos_ITRF[0] = 0;
	tel_pos_ITRF[1] = 0;
	tel_pos_ITRF[2] = 0;
	tel.height = 0;
	tel.latitude = 0;
	tel.longitude = 0;
	tel.pressure = 1000;
	tel.temperature = 0;

	//StandartTelescopeSpec
	num = 8;
	a_max = 1000.0;
	V_max = 10000.0;
	Gear_ratio = 720.0;
	circle_motor = 200 * num;
	circle_mount = circle_motor * Gear_ratio;
	step = 2 * pi / circle_mount;

	//StandartConv
	R.resize(9);
	R[0] = 0;
	R[1] = 0;
	R[2] = 1;
	R[3] = 0;
	R[4] = -1;
	R[5] = 0;
	R[6] = 1;
	R[7] = 0;
	R[8] = 0;

	//Standart steps for TLEpredict
	stepSec = 1.0;
	stepMin = 1.0 / 60.0;

	InitIntegrator();
}

//Деструктор
Convertor::~Convertor() {
}

//Инициализация интегратора
void Convertor::InitIntegrator() {
	IForce1 = new Force::InfluenceForce();
	IForce1->Init_CPU();
	POSat1.Init_CPU();

	//Считывание эфемерид движения планет
	short int error = 0;
	short int de_num = 0;
	double jd_beg, jd_end;
	if ((error = ephem_open("Data/eph/lnx1900.405", &jd_beg, &jd_end, &de_num)) != 0)
	{
		if (error == 1)
			printf("JPL ephemeris file not found.\n");
		else
			printf("Error reading JPL ephemeris file header.\n");
	}
	else
	{
		printf("JPL ephemeris DE%d open. Start JD = %10.2f  End JD = %10.2f\n",
			de_num, jd_beg, jd_end);
		printf("\n");
	}
}

//Инициализация характеристик монтировки
void Convertor::InitMountSpec(const double &st, const double &a, const double &V, const double &ratio, const double &steps) {
	num = st;
	a_max = a;
	V_max = V;
	Gear_ratio = ratio;
	circle_motor = steps * num;
	circle_mount = circle_motor * Gear_ratio;
	step = 2 * pi / circle_mount; 
}

//Выдать текущую позицию Az, Elev
std::pair<double, double> Convertor::GetAzElevPos() {
	return{ cur_Az, cur_Elev };
}

//Выдать текущую позицию Ra, Dec
std::pair<double, double> Convertor::GetRaDecPos() {
	return{ cur_Ra, cur_Dec };
}

//Выдать текущую позицию Alph, Bet
std::pair<double, double> Convertor::GetAlphBetPos() {
	return{ cur_Alph, cur_Bet };
}

//Выдать текущую позицию mot1, mot2
std::pair<double, double> Convertor::GetMotorsPos() {
	return{ cur_mot1, cur_mot2 };
}

//Задать позицию телескопа через x, y, z [км]
void Convertor::SetTelPosITRF(const double &x, const double &y, const double &z) {
	tel_pos_ITRF[0] = x;
	tel_pos_ITRF[1] = y;
	tel_pos_ITRF[2] = z;

	XYZ_WGS84(x, y, z, tel.height, tel.latitude, tel.longitude);
	//tel.height *= 1000;
}

//Задать позицию телескопа через массив x, y, z [км]
void Convertor::SetTelPosITRF(double *pos) {
	tel_pos_ITRF[0] = pos[0];
	tel_pos_ITRF[1] = pos[1];
	tel_pos_ITRF[2] = pos[2];

	XYZ_WGS84(pos[0], pos[1], pos[2], tel.height, tel.latitude, tel.longitude);
	//tel.height *= 1000;
}

//Задать позицию телескопа через lat, lon, elev [град, км]
void Convertor::SetTelPosLatLonElev(const double &Lat, const double &Lon, const double &Elev) {
	tel.height = Elev * 1000;
	tel.latitude = Lat;
	tel.longitude = Lon;

	WGS84_XYZ(Elev * 1000, Lat, Lon, tel_pos_ITRF[0], tel_pos_ITRF[1], tel_pos_ITRF[2]);
}

//Задать позицию телескопа через массив lat, lon, elev [град, км]
void Convertor::SetTelPosLatLonElev(double *pos) {
	tel.latitude = pos[0];
	tel.longitude = pos[1];
	tel.height = pos[2] * 1000;

	WGS84_XYZ(pos[2] * 1000, pos[0], pos[1], tel_pos_ITRF[0], tel_pos_ITRF[1], tel_pos_ITRF[2]);
}

//Задать температуру [С], давление [мбар] в точке наблюдения
void Convertor::SetTelTP(const double &T, const double &P) {
	tel.pressure = P;
	tel.temperature = T;
}

//Задать направление наблюдения через Jd, Az, Elev
void Convertor::SetAzElevPos(const double &Jd, const double &Az, const double &Elev) {
	cur_Jd = Jd;
	cur_Az = Az;
	cur_Elev = Elev;

	//Az, Elev -> Ra, Dec
	AzElev2RaDec();

	//Az, Elev -> Alph, Bet
	Convert2AlphBet();
}

//Задать направление наблюдения через Jd, Ra, Dec
void Convertor::SetRaDecPos(const double &Jd, const double &Ra, const double &Dec) {
	cur_Jd = Jd;
	cur_Ra = Ra;
	cur_Dec = Dec;

	//Ra, Dec -> Az,Elev
	RaDec2AzElev();

	//Az, Elev -> Alph, Bet
	Convert2AlphBet();
}

//Задать текущее положение двигателей mot1, mot2
void Convertor::SetMotorsPos(const double &mot1, const double &mot2) {
	cur_mot1 = mot1;
	cur_mot2 = mot2;
}

//Задать поправки ко времнеи и положение полюсов по Jd
bool Convertor::SetDateAndPolePos(const double &Jd) {
	//Get dAT
	double *tab = IForce1->TAIUTCCorrect;
	int kmax = IForce1->TAUUTC_kmax;
	bool yes = false;
	for (int i = 0; i <= kmax; i++)
	{
		int i4 = 4 * (kmax - i);
		if (tab[i4] <= Jd || i == kmax)
		{
			dAT = tab[i4 + 1];
			yes = true;
		}
	}

	if (!yes) {
		std::cerr << "Can't find dAT in TAI-UTC.DAT for current Jd = " << Jd << "!" << std::endl;
		return false;
	}

	//Calculate jd_tt
	jd_tt = Jd + (dAT + 32.184) / 86400.0;

	//Get dUT1, xp, yp
	double mjd = Jd - mjd0;
	auto it = IForce1->pole_offset.end();
	it--;
	while (abs(it->mjd_cur - mjd) > 1 && it != IForce1->pole_offset.begin()) {
		it--;
	}
	if (abs(it->mjd_cur - mjd) <= 1) {
		xp = it->px;
		yp = it->py;
		dUT1 = it->dUT1;
	}
	else {
		std::cerr << "Can't find xp and yp in finals.all for current Jd = " << Jd << "!" << std::endl;
		return false;
	}

	//Calculate jd_UT1 and delta_T
	jd_UT1 = Jd + dUT1 / 86400;
	deltaT = dAT + 32.184 - dUT1;

	return true;
}

//Задать матрицу преобразования Az,Elev -> Alph, Bet [3x3] 
void Convertor::SetConvMatr(const std::vector<double>& M) {
	for (int i = 0; i < M.size(); i++) {
		R[i] = M[i];
	}
}

//Az, Elev -> Alph, Bet
void Convertor::Convert2AlphBet() {
	std::vector<double> AzElev = { cur_Az, cur_Elev };
	
	std::vector<double> AlphBet = AzElev2AlphBet(AzElev);

	cur_Alph = AlphBet[0];
	cur_Bet = AlphBet[1];
}

//Az, Elev -> Alph, Bet
std::vector<double> Convertor::AzElev2AlphBet(const std::vector<double>& AzElev) {
	std::cout << "Az = " << AzElev[0] * 180 / M_PI << ", Elev = " << AzElev[1] * 180 / M_PI << std::endl;
	std::vector<double> xyz = AzElev2XYZ(AzElev);
	std::vector<double> xyz1 = Mat3x3XStolb3x1(R, xyz);

	std::vector<double> AlphBet(2);
	if (xyz1[2] != 1) {
		AlphBet = XYZ2AzElev(xyz1);
		//if (AzElev[0] >= 0 && AzElev[0] <= M_PI) {
		if (AlphBet[0] < 0) {
			AlphBet[0] = AlphBet[0] + M_PI;
			AlphBet[1] = AlphBet[1] - M_PI / 2;
		}
		else {
			AlphBet[1] = M_PI / 2 - AlphBet[1];
		}
	}
	else {
		AlphBet[0] = M_PI / 2;
		AlphBet[1] = 0;
	}

	//std::cout << "x = " << xyz[0] << ", y = " << xyz[1] << ", z = " << xyz[2] << std::endl;
	//std::vector<double> xyz2 = AzElev2XYZ(AlphBet);
	//std::vector<double> InvR = Inversion3x3(R);
	//std::vector<double> xyz_new = Mat3x3XStolb3x1(InvR, xyz2);
	//std::cout << "x_n = " << xyz_new[0] << ", y_n = " << xyz_new[1] << ", z_n = " << xyz_new[2] << std::endl;

	std::cout << "Alph = " << AlphBet[0] * 180 / M_PI << ", Bet = " << AlphBet[1] * 180 / M_PI << std::endl;

	return AlphBet;
}

//Ra, Dec -> Az, Elev
void Convertor::RaDec2AzElev() {
	//Set Date and Pole position
	SetDateAndPolePos(cur_Jd);

	//Ra, Dec -> Elev, Az 
	double corRa, corDec;
	double zd, az;
	double tmp_Ra, tmp_Dec;
	gcrs2equ(jd_tt, 1, 0, cur_Ra * 12 / pi, cur_Dec * 180 / pi, &tmp_Ra, &tmp_Dec);
	equ2hor(jd_UT1, deltaT, 0, xp, yp, &tel, tmp_Ra, tmp_Dec, 0, &zd, &az, &corRa, &corDec);
	cur_Elev = (90 - zd) * pi / 180;
	cur_Az = az * pi / 180;
}

//Az, Elev -> Ra, Dec
void Convertor::AzElev2RaDec() {
	//Set Date and Pole position
	SetDateAndPolePos(cur_Jd);

	//Ra, Dec -> Elev, Az
	double poster[3];
	double poscel[3];
	double Ra, Dec;
	AzElevR2XYZ_ITRF(poster);
	ter2cel(jd_UT1, 0, deltaT, 1, 0, 0, xp, yp, poster, poscel);
	vector2radec(poscel, &Ra, &Dec);
	cur_Ra = Ra / 12 * pi;
	cur_Dec = Dec / 180 * pi;
}

double Modulus(double x, double y)
{
	double modu;
	modu = x - (int)(x / y) * y;		// (int) <-> trunc() ??
	if (modu >= 0)
		return modu;
	else
		return (modu + y);
}

std::vector<double> AzElev2XYZ(const std::vector<double> &AzElev) {
	std::vector<double> xyz(3);
	xyz[0] = cos(AzElev[1]) * cos(AzElev[0]);
	xyz[1] = cos(AzElev[1]) * sin(AzElev[0]);
	xyz[2] = sin(AzElev[1]);
	return xyz;
}

std::vector<double> XYZ2AzElev(const std::vector<double> &xyz) {
	std::vector<double> AzElev(2);
	AzElev[0] = atan(xyz[1] / xyz[0]);
	AzElev[1] = asin(xyz[2]);
	return AzElev;
}

//Az, Elev -> XYZ_ITRF
void Convertor::AzElevR2XYZ_ITRF(double* pos) {
	double sindec = sin(cur_Elev) * sin(tel.latitude * pi / 180) + cos(cur_Elev) * cos(tel.latitude * pi / 180) * cos(cur_Az);
	double dec = asin(sindec);

	double lha = atan2(-sin(cur_Az) * cos(cur_Elev) / cos(dec), (sin(cur_Elev) - sin(dec) * sin(tel.latitude * pi / 180)) / (cos(dec) * cos(tel.latitude * pi / 180)));
	double ra = Modulus(tel.longitude * pi / 180 - lha, 2 * M_PI);

	double R = 1.0;
	pos[0] = R * cos(dec) * cos(ra);
	pos[1] = R * cos(dec) * sin(ra);
	pos[2] = R * sin(dec);
}

//Вычисление траектории из заданного положения inAlph, inBet в положение outAlph, outBet
traject Convertor::CalcTraject(const double &inAlph, const double &inBet, const double &outAlph, const double &outBet) {
	traject tr;
	tr.startpos.first = cur_mot1;
	tr.startpos.second = cur_mot2;
	
	//Общее кол-во шагов по двум осям
	double stepsA = (outAlph - inAlph) / step;
	double stepsB = (outBet - inBet) / step;
	tr.endpos.first = cur_mot1 + stepsA;
	tr.endpos.second = cur_mot2 + stepsB;

	double signA = sign(stepsA);
	double signB = sign(stepsB);
	stepsA = abs(stepsA);
	stepsB = abs(stepsB);

	tr.a.first = signA * a_max;
	tr.a.second = signB * a_max;

	//Время ускорения + торможения
	double Ta = 2 * (V_max / a_max);
	//Расстояние пройденное за ускорение + торможение
	double Sa = 2 * (V_max * V_max / (2 * a_max));
	
	double T_A;
	//Если выходит на максимальную скорость
	if (stepsA > Sa) {
		tr.V.first = signA * V_max;
		T_A = Ta + (stepsA - Sa) / V_max;
	}
	//Если не может достигнуть максималдьной скорости
	else {
		double V = sqrt(2 * (stepsA / 2) * a_max);
		tr.V.first = signA * V;
		T_A = 2 * V / a_max;
	}

	double T_B;
	//Если выходит на максимальную скорость
	if (stepsB > Sa) {
		tr.V.second = signB * V_max;
		T_B = Ta + (stepsB - Sa) / V_max;
	}
	//Если не может достигнуть максимальной скорости
	else {
		double V = sqrt(2 * (stepsB / 2) * a_max);
		tr.V.second = signB * V;
		T_B = 2 * V / a_max;
	}

	//Выбрать максимальное время
	if (T_A >= T_B) {
		tr.T = T_A;
	}
	else {
		tr.T = T_B;
	}
	
	return tr;
}

//Переехать из текущего положения cur_Jd, cur_Ra, cur_Dec в заданное Jd, outRa, outDec
traject Convertor::GoToRaDec(const double &Jd, const double &outRa, const double &outDec) {
	double inAlph = cur_Alph;
	double inBet = cur_Bet;
	SetRaDecPos(Jd, outRa, outDec);
	traject tr = CalcTraject(inAlph, inBet, cur_Alph, cur_Bet);
	SetMotorsPos(tr.endpos.first, tr.endpos.second);

	return tr;
}

//Переехать из текущего положения cur_Jd, cur_Az, cur_Elev в заданное Jd, outAz, outElev
traject Convertor::GoToAzElev(const double &Jd, const double &outAz, const double &outElev) {
	double inAlph = cur_Alph;
	double inBet = cur_Bet;
	SetAzElevPos(Jd, outAz, outElev);
	traject tr = CalcTraject(inAlph, inBet, cur_Alph, cur_Bet);
	SetMotorsPos(tr.endpos.first, tr.endpos.second);

	return tr;
}

//Проверка ограничений
int Convertor::LimitsAzElev(const double &outAz, const double &outElev) {
	//Проверка начального положения
	if ((cur_Az < 0 && cur_Az > 2 * pi) || (cur_Elev < 0 && cur_Elev > pi / 2)) {
		return -1;
	}
	//Проверка конечного положения
	if ((outAz < 0 && outAz > 2 * pi) || (outElev < 0 && outElev > pi / 2)) {
		return -2;
	}
	//Проверка пересечения линии севера
	std::vector<double> AlphBet = AzElev2AlphBet({ outAz, outElev });
	if ((cur_Alph < pi / 2) && (AlphBet[0] > pi / 2) || (cur_Alph > pi / 2) && (AlphBet[0] < pi / 2)) {
		return 1;
	}
	return 0;
}

//Задать TLE параметры
void Convertor::SetTLEparams(const std::string &path, const int &num_str) {
	TLEpath = path;
	numb_of_line = num_str;
	tleLoader.LoadData(TLEpath.c_str(), numb_of_line);
}

//Задать интервал между точками в сек
void Convertor::SetStep(const double &step) {
	stepSec = step;
	stepMin = stepSec / 60.0;
}

//Выбрать орбиту NaradId
bool Convertor::GetOrb(const int &NoradId) {
	tleLoader.GetFromID(NoradId);
	if (orb != NULL) {
		return true;
	}
	else {
		return false;
	}
}

//Вычислить массив углов для сопровождения
std::vector<AccompAngs> Convertor::CalculateAngs(const int &NoradId, const double &JdStart, const double &JdEnd) {
	std::vector<AccompAngs> angs;

	cSite siteView(tel.latitude, tel.longitude, tel.height / 1000);

	if (GetOrb(NoradId)) {
		double JdEpoch = orb->Epoch().Date();
		for (double Jd = JdStart; Jd < JdEnd; Jd + stepSec / 86400) {
			//Количество времени от записи tle в мин.
			double mpe = (Jd - JdEpoch) * 24 * 60;

			// Get the position of the satellite at time "mpe"
			cEciTime eci = orb->GetPosition(mpe);

			// Now get the "look angle" from the site to the satellite. 
			// Note that the ECI object "eciSDP4" contains a time associated
			// with the coordinates it contains; this is the time at which
			// the look angle is valid.
			cTopo topoLook = siteView.GetLookAngle(eci);

			angs.push_back({ Jd, topoLook.AzimuthRad(), topoLook.ElevationRad() });
		}
	}

	return angs;
}

//Вычисление массива точек для микроконтроллера
std::vector<AccompPoints> Convertor::CalculatePoints(const std::vector<AccompAngs> &ags) {
	std::vector<AccompPoints> points;

	return points;
}