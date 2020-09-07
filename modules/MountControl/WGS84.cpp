#include "WGS84.h"

#include <math.h>

void WGS84_XYZ(double Hw, double Fwg, double Lwg, double & X, double & Y, double & Z)
{
	double n;
	long double cF, sF, cL, sL;

	SinCos(Fwg*GR, sF, cF);
	SinCos(Lwg*GR, sL, cL);
	n = R0Earth / sqrt(1 - E2Earth*Sqr(sF));
	Z = (n*(1 - E2Earth) + Hw)*sF;
	X = (n + Hw)*cF*cL;
	Y = (n + Hw)*cF*sL;
}

void XYZ_WGS84(double X, double Y, double Z, double & Hw, double & Fwg, double & Lwg)
{
	double R1, R, sh, dol;
	R1 = sqrt(X*X + Y*Y);
	R = sqrt(R1*R1 + Z*Z);
	sh = asin(Z / R);
	dol = acos(X / R1)*sign(Y);
	Lwg = dol * RG;
	RF_WGS84(R, sh, Hw, Fwg); // перевод из геоцентр. СК в WGS-84
}

void RF_WGS84(double R, double F, double & Hz, double & Fzg)
{
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
		n = R0Earth / sqrt(1 - E2Earth*Sqr(sF));
		Z = (n*(1 - E2Earth) + Hz)*sF;
		rW = sqrt(Sqr((n + Hz)*cF) + Sqr(Z));
		FW = asin(Z / rW);
		Hz = Hz - (rW - R);
		Fzg = Fzg - (FW - F);
	} while ((fabs(rW - R) >= 1e-8));
	Fzg *= RG;
}

void SinCos(double a, long double & sina, long double & cosa)
{
	sina = sin(a);
	cosa = cos(a);
}

double sign(double Val)
{
	if (Val == 0.)  return 0;
	if (Val > 0.)  return 1;
	else return -1;
}

double Sqr(double a)
{
	return a*a;
}