#include "Convertor.h"
#include "WGS84.h"
#include "Mat3x3.h"

#include <iostream>
#include <math.h>

//Конструктор
Convertor::Convertor() {
	//CurrentPosition
	curRaDec = { 0, 0, 0 };
	curAzElev = { 0, 0, 0 };
	curAlphBet = { 0, 0, 0 };
	curMotors = { 0, 0 };

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
Angs Convertor::GetAzElevPos() {
	return curAzElev;
}

//Выдать текущую позицию Ra, Dec
Angs Convertor::GetRaDecPos() {
	return curRaDec;
}

//Выдать текущую позицию Alph, Bet
Angs Convertor::GetAlphBetPos() {
	return curAlphBet;
}

//Выдать текущую позицию mot1, mot2
MotPos Convertor::GetMotorsPos() {
	return curMotors;
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
void Convertor::SetAzElevPos(const Angs &AzElev) {
	curAzElev = AzElev;

	//Az, Elev -> Ra, Dec
	AzElev2RaDec();

	//Az, Elev -> Alph, Bet
	Convert2AlphBet();
}

//Задать направление наблюдения через Jd, Ra, Dec
void Convertor::SetRaDecPos(const Angs &RaDec) {
	curRaDec = RaDec;

	//Ra, Dec -> Az,Elev
	RaDec2AzElev();

	//Az, Elev -> Alph, Bet
	Convert2AlphBet();
}

//Задать текущее положение двигателей mot1, mot2
void Convertor::SetMotorsPos(const MotPos &Motors) {
	curMotors = Motors;
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
	Angs AzElev = curAzElev;
	
	curAlphBet = AzElev2AlphBet(AzElev);
}

//Az, Elev -> Alph, Bet
Angs Convertor::AzElev2AlphBet(const Angs &AzElev) {
	std::cout << "Az = " << AzElev.ang1 * 180 / M_PI << ", Elev = " << AzElev.ang2 * 180 / M_PI << std::endl;
	std::vector<double> xyz = AzElev2XYZ(AzElev);
	std::vector<double> xyz1 = Mat3x3XStolb3x1(R, xyz);

	Angs AlphBet;
	if (xyz1[2] != 1) {
		AlphBet = XYZ2AzElev(xyz1);
		//if (AzElev[0] >= 0 && AzElev[0] <= M_PI) {
		if (AlphBet.ang1 < 0) {
			AlphBet.ang1 = AlphBet.ang1 + M_PI;
			AlphBet.ang2 = AlphBet.ang2 - M_PI / 2;
		}
		else {
			AlphBet.ang2 = M_PI / 2 - AlphBet.ang2;
		}
	}
	else {
		AlphBet.ang1 = M_PI / 2;
		AlphBet.ang2 = 0;
	}
	AlphBet.Jd = AzElev.Jd;

	//std::cout << "x = " << xyz[0] << ", y = " << xyz[1] << ", z = " << xyz[2] << std::endl;
	//std::vector<double> xyz2 = AzElev2XYZ(AlphBet);
	//std::vector<double> InvR = Inversion3x3(R);
	//std::vector<double> xyz_new = Mat3x3XStolb3x1(InvR, xyz2);
	//std::cout << "x_n = " << xyz_new[0] << ", y_n = " << xyz_new[1] << ", z_n = " << xyz_new[2] << std::endl;

	std::cout << "Alph = " << AlphBet.ang1 * 180 / M_PI << ", Bet = " << AlphBet.ang2 * 180 / M_PI << std::endl;

	return AlphBet;
}

//Alph, Bet -> Az, Elev
Angs Convertor::AlphBet2AzElev(const Angs &AlphBet) {
	Angs tmpAlphBet = AlphBet;
	if (tmpAlphBet.ang1 > M_PI / 2) {
		tmpAlphBet.ang1 = tmpAlphBet.ang1 - M_PI;
		tmpAlphBet.ang2 = tmpAlphBet.ang2 + M_PI / 2;
	}
	else {
		tmpAlphBet.ang2 = M_PI / 2 - tmpAlphBet.ang2;
	}
	std::vector<double> xyz = AzElev2XYZ(tmpAlphBet);
	std::vector<double> InvR = Inversion3x3(R);
	std::vector<double> xyzAzElev = Mat3x3XStolb3x1(InvR, xyz);
	Angs AzElev = XYZ2AzElev(xyzAzElev);
	AzElev.Jd = AlphBet.Jd;
	
	if (AlphBet.ang2 < -M_PI / 2 || AlphBet.ang2 > M_PI / 2) {
		AzElev.ang1 += M_PI;
	}
	else if (AlphBet.ang2 < M_PI / 2 && AlphBet.ang1 < M_PI / 2) {
		AzElev.ang1 += 2 * M_PI;
	}

	return AzElev;
}

//Alph, Bet -> Motors
MotPos Convertor::AlphBet2Motors(const Angs &AlphBet) {
	MotPos Motors;
	Motors.mot1 = AlphBet.ang1 / step;
	Motors.mot2 = AlphBet.ang2 / step;
	return Motors;
}

//Az, Elev -> Motors
MotPos Convertor::AzElev2Motors(const Angs &AzElev) {
	return AlphBet2Motors(AzElev2AlphBet(AzElev));
}

//Ra, Dec -> Az, Elev
void Convertor::RaDec2AzElev() {
	//Set Date and Pole position
	SetDateAndPolePos(curRaDec.Jd);

	//Ra, Dec -> Elev, Az 
	double corRa, corDec;
	double zd, az;
	double tmp_Ra, tmp_Dec;
	gcrs2equ(jd_tt, 1, 0, curRaDec.ang1 * 12 / pi, curRaDec.ang2 * 180 / pi, &tmp_Ra, &tmp_Dec);
	equ2hor(jd_UT1, deltaT, 0, xp, yp, &tel, tmp_Ra, tmp_Dec, 0, &zd, &az, &corRa, &corDec);
	curAzElev.ang2 = (90 - zd) * pi / 180;
	curAzElev.ang1 = az * pi / 180;
	curAzElev.Jd = curRaDec.Jd;
}

//Az, Elev -> Ra, Dec
void Convertor::AzElev2RaDec() {
	//Set Date and Pole position
	SetDateAndPolePos(curAzElev.Jd);

	//Ra, Dec -> Elev, Az
	double poster[3];
	double poscel[3];
	double Ra, Dec;
	AzElevR2XYZ_ITRF(poster);
	ter2cel(jd_UT1, 0, deltaT, 1, 0, 0, xp, yp, poster, poscel);
	vector2radec(poscel, &Ra, &Dec);
	curRaDec.ang1 = Ra / 12 * pi;
	curRaDec.ang2 = Dec / 180 * pi;
	curRaDec.Jd = curAzElev.Jd;
}

double Modulus_(double x, double y)
{
	double modu;
	modu = x - (int)(x / y) * y;		// (int) <-> trunc() ??
	if (modu >= 0)
		return modu;
	else
		return (modu + y);
}

std::vector<double> AzElev2XYZ(const Angs &AzElev) {
	std::vector<double> xyz(3);
	xyz[0] = cos(AzElev.ang2) * cos(AzElev.ang1);
	xyz[1] = cos(AzElev.ang2) * sin(AzElev.ang1);
	xyz[2] = sin(AzElev.ang2);
	return xyz;
}

Angs XYZ2AzElev(const std::vector<double> &xyz) {
	Angs AzElev;
	AzElev.ang1 = atan(xyz[1] / xyz[0]);
	AzElev.ang2 = asin(xyz[2]);
	return AzElev;
}

//Az, Elev -> XYZ_ITRF
void Convertor::AzElevR2XYZ_ITRF(double* pos) {
	double sindec = sin(curAzElev.ang2) * sin(tel.latitude * pi / 180) + cos(curAzElev.ang2) * cos(tel.latitude * pi / 180) * cos(curAzElev.ang1);
	double dec = asin(sindec);

	double lha = atan2(-sin(curAzElev.ang1) * cos(curAzElev.ang2) / cos(dec), (sin(curAzElev.ang2) - sin(dec) * sin(tel.latitude * pi / 180)) / (cos(dec) * cos(tel.latitude * pi / 180)));
	double ra = Modulus_(tel.longitude * pi / 180 - lha, 2 * M_PI);

	double R = 1.0;
	pos[0] = R * cos(dec) * cos(ra);
	pos[1] = R * cos(dec) * sin(ra);
	pos[2] = R * sin(dec);
}

//Вычисление траектории из заданного положения inAlph, inBet в положение outAlph, outBet
traject Convertor::CalcTraject(const Angs &inAlphBet, const Angs &outAlphBet) {
	traject tr;
	tr.startpos.first = curMotors.mot1;
	tr.startpos.second = curMotors.mot2;
	
	//Общее кол-во шагов по двум осям
	double stepsA = (outAlphBet.ang1 - inAlphBet.ang1) / step;
	double stepsB = (outAlphBet.ang2 - inAlphBet.ang2) / step;
	tr.endpos.first = curMotors.mot1 + stepsA;
	tr.endpos.second = curMotors.mot2 + stepsB;

	double signA = sign_(stepsA);
	double signB = sign_(stepsB);
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
traject Convertor::GoToRaDec(const Angs &outRaDec) {
	Angs inAlphBet = curAlphBet;
	SetRaDecPos(outRaDec);
	traject tr = CalcTraject(inAlphBet, curAlphBet);
	SetMotorsPos({ tr.endpos.first, tr.endpos.second });

	return tr;
}

//Переехать из текущего положения cur_Jd, cur_Az, cur_Elev в заданное Jd, outAz, outElev
traject Convertor::GoToAzElev(const Angs &outAzElev) {
	Angs inAlphBet = curAlphBet;
	SetAzElevPos(outAzElev);
	traject tr = CalcTraject(inAlphBet, curAlphBet);
	SetMotorsPos({ tr.endpos.first, tr.endpos.second });

	return tr;
}

//Проверка ограничений
int Convertor::LimitsAzElev(const Angs &outAzElev) {
	//Проверка начального положения
	if ((curAzElev.ang1 < 0 && curAzElev.ang1 > 2 * pi) || (curAzElev.ang2 < 0 && curAzElev.ang2 > pi / 2)) {
		return -1;
	}
	//Проверка конечного положения
	if ((outAzElev.ang1 < 0 && outAzElev.ang1 > 2 * pi) || (outAzElev.ang2 < 0 && outAzElev.ang2 > pi / 2)) {
		return -2;
	}
	//Проверка пересечения линии севера
	Angs AlphBet = AzElev2AlphBet(outAzElev);
	if ((curAlphBet.ang1 < pi / 2) && (AlphBet.ang1 > pi / 2) || (curAlphBet.ang1 > pi / 2) && (AlphBet.ang1 < pi / 2)) {
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
std::vector<Angs> Convertor::CalculateAzElev(const int &NoradId, const double &JdStart, const double &JdEnd) {
	std::vector<Angs> angs;

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

//Вычисление массива углов и угловых скоростей AlphBet
std::vector<AngsVels> Convertor::CalculateAlphBetVels(const int &NoradId, const double &JdStart, const double &JdEnd) {
	std::vector<AngsVels> AlphBetArr;
	cSite siteView(tel.latitude, tel.longitude, tel.height / 1000);

	if (GetOrb(NoradId)) {
		double JdEpoch = orb->Epoch().Date();
		for (double Jd = JdStart; Jd < JdEnd; Jd + stepSec / 86400) {
			AngsVels AlphBet;
			double mpe = (Jd - JdEpoch) * 24 * 60;

			// Get the position of the satellite at time "mpe"
			cEciTime eci = orb->GetPosition(mpe);

			// Now get the "look angle" from the site to the satellite. 
			// Note that the ECI object "eciSDP4" contains a time associated
			// with the coordinates it contains; this is the time at which
			// the look angle is valid.
			cTopo topoLook = siteView.GetLookAngle(eci);

			Angs AzElev = { Jd, topoLook.AzimuthRad(), topoLook.ElevationRad() };
			Angs tmpAlphBet = AzElev2AlphBet(AzElev);
			AlphBet.ang1 = tmpAlphBet.ang1;
			AlphBet.ang2 = tmpAlphBet.ang2;
			AlphBet.Jd = tmpAlphBet.Jd;

			double h = 0.1;
			eci = orb->GetPosition(mpe - h / 60.0);
			topoLook = siteView.GetLookAngle(eci);
			AzElev = { Jd, topoLook.AzimuthRad(), topoLook.ElevationRad() };
			tmpAlphBet = AzElev2AlphBet(AzElev);

			eci = orb->GetPosition(mpe + h / 60.0);
			topoLook = siteView.GetLookAngle(eci);
			AzElev = { Jd, topoLook.AzimuthRad(), topoLook.ElevationRad() };
			Angs tmpAlphBet1 = AzElev2AlphBet(AzElev);

			AlphBet.ang1_vel = (tmpAlphBet1.ang1 - tmpAlphBet.ang1) / (2 * h);
			AlphBet.ang2_vel = (tmpAlphBet1.ang2 - tmpAlphBet.ang2) / (2 * h);
			AlphBetArr.push_back(AlphBet);
		}
	}

	return AlphBetArr;
}

//Вычисление массива точек для микроконтроллера
std::pair<std::vector<AccompPoints>, std::vector<AccompPoints>> Convertor::CalculatePoints(const std::vector<AngsVels> &AlphBetArr) {
	std::vector<AccompPoints> points1(AlphBetArr.size() + 2);
	std::vector<AccompPoints> points2(AlphBetArr.size() + 2);

	MotPos motors_tmp = AlphBet2Motors(AlphBetArr[0]);
	points1[0].V0 = 0;
	points1[0].T = 0;

	points2[0].V0 = 0;
	points2[0].T = 0;

	for (int i = 0; i < AlphBetArr.size(); i++) {
		MotPos motors = AlphBet2Motors(AlphBetArr[i]);
		points1[i + 1].N = motors.mot1;
		points2[i + 1].N = motors.mot2;
		points1[i + 1].V0 = AlphBetArr[i].ang1_vel / step;
		points2[i + 1].V0 = AlphBetArr[i].ang2_vel / step;
		if (i == 0) {
			points1[1].T = points1[i].V0 / a_max;
			points2[1].T = points2[i].V0 / a_max;
		}
		else {
			points1[i + 1].T = points1[i].T + 1;
			points2[i + 1].T = points2[i].T + 1;
		}
	}

	points1[0].N = points1[1].N - points1[1].V0 * points1[1].V0 / (2 * a_max);
	points2[0].N = points2[1].N - points2[1].V0 * points2[1].V0 / (2 * a_max);

	points1.back().V0 = 0;
	points1.back().T = points1[AlphBetArr.size()].T + abs(points1[AlphBetArr.size()].V0) / a_max;
	points1.back().N = points1[AlphBetArr.size()].N + (abs(points1[AlphBetArr.size()].V0) / a_max) * points1[AlphBetArr.size()].V0 / 2;

	points2.back().V0 = 0;
	points2.back().T = points2[AlphBetArr.size()].T + abs(points2[AlphBetArr.size()].V0) / a_max;
	points1.back().N = points2[AlphBetArr.size()].N + (abs(points2[AlphBetArr.size()].V0) / a_max) * points2[AlphBetArr.size()].V0 / 2;

	return{ points1, points2 };
}