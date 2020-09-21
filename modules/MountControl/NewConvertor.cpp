#include <iostream>

#include "NewConvertor.h"

NewConvertor::NewConvertor() {
}

NewConvertor::~NewConvertor() {
}

//Задать поправки ко времнеи и положение полюсов по Jd
bool NewConvertor::SetDateAndPolePos(double Jd) {
	//Get dAT
	double* tab = IForce1->TAIUTCCorrect;
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

Angs NewConvertor::AzElev2RaDec(double Jd, const Angs& AzElev, const on_surface& pos) {
	//Set Date and Pole position
	SetDateAndPolePos(Jd);

	//Ra, Dec -> Elev, Az
	double poster[3];
	double poscel[3];
	double Ra, Dec;
	AzElevR2ITRF(AzElev.ang1, AzElev.ang2, 1.0, pos, poster);
	ter2cel(jd_UT1, 0, deltaT, 1, 0, 0, xp, yp, poster, poscel);
	vector2radec(poscel, &Ra, &Dec);
	double ang1 = Ra / 12 * M_PI;
	double ang2 = Dec / 180 * M_PI;
	
	return { ang1, ang2 };
}

Angs NewConvertor::RaDec2AzElev(double Jd, const Angs& RaDec, const on_surface& pos) {
	//Set Date and Pole position
	SetDateAndPolePos(Jd);

	//Ra, Dec -> Elev, Az 
	double corRa, corDec;
	double zd, az;
	double tmp_Ra, tmp_Dec;
	on_surface tmp_pos = pos;
	gcrs2equ(jd_tt, 1, 0, RaDec.ang1 * 12 / M_PI, RaDec.ang2 * 180 / M_PI, &tmp_Ra, &tmp_Dec);
	equ2hor(jd_UT1, deltaT, 0, xp, yp, &tmp_pos, tmp_Ra, tmp_Dec, 0, &zd, &az, &corRa, &corDec);
	double ang1 = az * M_PI / 180;
	double ang2 = (90 - zd) * M_PI / 180;
	
	return { ang1, ang2 };
}

void NewConvertor::AzElevR2ITRF(double Az, double Elev, double R, const on_surface& pos, double* ITRF) const {
	std::vector<double> xyz = AzElevR2ITRF(Az, Elev, R, pos);
	
	ITRF[0] = xyz[0];
	ITRF[1] = xyz[1];
	ITRF[2] = xyz[2];
}

std::vector<double> NewConvertor::AzElevR2ITRF(double Az, double Elev, double R, const on_surface& pos) const {
	double sindec = sin(Elev) * sin(pos.latitude * M_PI / 180) + cos(Elev) * cos(pos.latitude * M_PI / 180) * cos(Az);
	double dec = asin(sindec);

	double lha = atan2(-sin(Az) * cos(Elev) / cos(dec), (sin(Elev) - sin(dec) * sin(pos.latitude * M_PI / 180)) / (cos(dec) * cos(pos.latitude * M_PI / 180)));
	double ra = Modulus(pos.longitude * M_PI / 180 - lha, 2 * M_PI);

	return SphereCoord2XYZ({ R, ra, dec });
}

SphereCoord NewConvertor::XYZ2SphereCoord(const std::vector<double>& XYZ) const {
	double R = sqrt(XYZ[0] * XYZ[0] + XYZ[1] * XYZ[1] + XYZ[2] * XYZ[2]);
	double ang1 = atan(XYZ[1] / XYZ[0]);
	double ang2 = asin(XYZ[2]);
	return { R, ang1, ang2 };
}

std::vector<double> NewConvertor::SphereCoord2XYZ(const SphereCoord& sphere) const {
	double x = sphere.R * cos(sphere.ang2) * cos(sphere.ang1);
	double y = sphere.R * cos(sphere.ang2) * sin(sphere.ang1);
	double z = sphere.R * sin(sphere.ang2);
	return { x, y, z };
}

void NewConvertor::WGS84_XYZ(double Hw, double Fwg, double Lwg, double& X, double& Y, double& Z) const {
	double n;
	long double cF, sF, cL, sL;

	SinCos(Fwg * GR, sF, cF);
	SinCos(Lwg * GR, sL, cL);
	n = R0Earth / sqrt(1 - E2Earth * Sqr(sF));
	Z = (n * (1 - E2Earth) + Hw) * sF;
	X = (n + Hw) * cF * cL;
	Y = (n + Hw) * cF * sL;
}

void NewConvertor::XYZ_WGS84(double X, double Y, double Z, double& Hw, double& Fwg, double& Lwg) const {
	double R1, R, sh, dol;
	R1 = sqrt(X * X + Y * Y);
	R = sqrt(R1 * R1 + Z * Z);
	sh = asin(Z / R);
	dol = acos(X / R1) * sign(Y);
	Lwg = dol * RG;
	RF_WGS84(R, sh, Hw, Fwg); // перевод из геоцентр. СК в WGS-84
}

void NewConvertor::RF_WGS84(double R, double F, double& Hz, double& Fzg) const {
	long double sF, cF, n, Z;
	double rW, FW;

	if (R <= 0) {
		Hz = 0;
		Fzg = 0;
		return;
	}
	Hz = R - RMEarth;
	Fzg = F;

	do {
		SinCos(Fzg, sF, cF);
		n = R0Earth / sqrt(1 - E2Earth * Sqr(sF));
		Z = (n * (1 - E2Earth) + Hz) * sF;
		rW = sqrt(Sqr((n + Hz) * cF) + Sqr(Z));
		FW = asin(Z / rW);
		Hz = Hz - (rW - R);
		Fzg = Fzg - (FW - F);
	} while ((fabs(rW - R) >= 1e-8));
	Fzg *= RG;
}

double Modulus(double x, double y) {
	double modu;
	modu = x - (int)(x / y) * y;		// (int) <-> trunc() ??
	if (modu >= 0)
		return modu;
	else
		return (modu + y);
}

void SinCos(double a, long double& sina, long double& cosa) {
	sina = sin(a);
	cosa = cos(a);
}

double sign(double Val) {
	if (Val == 0.)  return 0;
	if (Val > 0.)  return 1;
	else return -1;
}

double Sqr(double a) {
	return a * a;
}