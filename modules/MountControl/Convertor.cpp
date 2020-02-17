#include "Convertor.h"
#include "WGS84.h"
#include "Mat3x3.h"

#include <iostream>

Convertor::Convertor() {
	//CurrentPosition
	cur_Jd = 0;
	cur_Az = 0;
	cur_Elev = 0;
	cur_Ra = 0;
	cur_Dec = 0;
	cur_Alph = 0;
	cur_Bet = 0;

	//TelescopePosition
	tel_pos_ITRF[0] = 0;
	tel_pos_ITRF[1] = 0;
	tel_pos_ITRF[2] = 0;
	tel.height = 0;
	tel.latitude = 0;
	tel.longitude = 0;
	tel.pressure = 1000;
	tel.temperature = 0;

	//StandartConv
	A.resize(6);
	A[0] = 1.0;
	A[1] = 0.0;
	A[2] = 0.0;
	A[3] = 0.0;
	A[4] = 1.0;
	A[5] = 0.0;

	InitIntegrator();
}

Convertor::~Convertor() {
}

void Convertor::InitIntegrator() {
	IForce1 = new Force::InfluenceForce();
	IForce1->Init_CPU();
	POSat1.Init_CPU();

	//—читывание эфемерид движени€ планет
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

std::pair<double, double> Convertor::GetAzElevPos() {
	return{ cur_Az, cur_Elev };
}

std::pair<double, double> Convertor::GetRaDecPos() {
	return{ cur_Ra, cur_Dec };
}

std::pair<double, double> Convertor::GetAlphBetPos() {
	return{ cur_Alph, cur_Bet };
}

void Convertor::SetTelPosITRF(const double &x, const double &y, const double &z) {
	tel_pos_ITRF[0] = x;
	tel_pos_ITRF[1] = y;
	tel_pos_ITRF[2] = z;

	XYZ_WGS84(x, y, z, tel.height, tel.latitude, tel.longitude);
	//tel.height *= 1000;
}

void Convertor::SetTelPosITRF(double *pos) {
	tel_pos_ITRF[0] = pos[0];
	tel_pos_ITRF[1] = pos[1];
	tel_pos_ITRF[2] = pos[2];

	XYZ_WGS84(pos[0], pos[1], pos[2], tel.height, tel.latitude, tel.longitude);
	//tel.height *= 1000;
}

void Convertor::SetTelPosLatLonElev(const double &Lat, const double &Lon, const double &Elev) {
	tel.height = Elev * 1000;
	tel.latitude = Lat;
	tel.longitude = Lon;

	WGS84_XYZ(Elev * 1000, Lat, Lon, tel_pos_ITRF[0], tel_pos_ITRF[1], tel_pos_ITRF[2]);
}

void Convertor::SetTelPosLatLonElev(double *pos) {
	tel.latitude = pos[0];
	tel.longitude = pos[1];
	tel.height = pos[2] * 1000;

	WGS84_XYZ(pos[2] * 1000, pos[0], pos[1], tel_pos_ITRF[0], tel_pos_ITRF[1], tel_pos_ITRF[2]);
}

void Convertor::SetTelTP(const double &T, const double &P) {
	tel.pressure = P;
	tel.temperature = T;
}

void Convertor::SetAzElevPos(const double &Jd, const double &Az, const double &Elev) {
	cur_Jd = Jd;
	cur_Az = Az;
	cur_Elev = Elev;

	//Az, Elev -> Ra, Dec
	AzElev2RaDec();

	//Az, Elev -> Alph, Bet
	Convert2AlphBet();
}

void Convertor::SetRaDecPos(const double &Jd, const double &Ra, const double &Dec) {
	cur_Jd = Jd;
	cur_Ra = Ra;
	cur_Dec = Dec;

	//Ra, Dec -> Az,Elev
	RaDec2AzElev();

	//Az, Elev -> Alph, Bet
	Convert2AlphBet();
}

void Convertor::SetAlphBetPos(const double &Jd, const double &Alph, const double &Bet) {
	cur_Jd = Jd;
	cur_Alph = Alph;
	cur_Bet = Bet;

	//Alph, Bet -> Az, Elev
	std::vector<double> InvA(4);
	std::vector<double> A_(4);
	A_[0] = A[0];
	A_[1] = A[1];
	A_[2] = A[3];
	A_[3] = A[4];
	InvA = Inversion2x2(A_);
	cur_Az = InvA[0] * (Alph - A[2]) + InvA[1] * (Alph - A[2]);
	cur_Elev = InvA[2] * (Bet - A[5]) + InvA[3] * (Bet - A[5]);

	//Az, Elev -> Ra, Dec
	AzElev2RaDec();
}

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

void Convertor::SetConvMatr(const std::vector<double>& M) {
	for (int i = 0; i < M.size(); i++) {
		A[i] = M[i];
	}
}

void Convertor::Convert2AlphBet() {
	cur_Alph = A[0] * cur_Az + A[1] * cur_Elev + A[2];
	cur_Bet = A[3] * cur_Az + A[4] * cur_Elev + A[5];
}

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