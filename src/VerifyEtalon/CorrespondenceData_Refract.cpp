//==============================================================================//
// вычисление рефракции через температуру и давление
//==============================================================================//
#include "CorrespondenceData.h"

//==============================================================================//
//
//==============================================================================//
double CorrespondenceData::refractivity( double pMillibar, double tKelvin )
{
	double result =77.6*pMillibar/tKelvin;
	return result;
}
//==============================================================================//
//
//==============================================================================//
void CorrespondenceData::calculateF( double *F1p, double *F2p, double *F3p, double *F4p, double F0, double F1, double FF1, double FF2 )
{
	double F1_ = FF1;
	double F2_ = FF2/F1_-F1_;
	double F3_ = F2_/(F0*F0*F1_*(1+F1_/F2_)-(1+F1*F1_));
	double F4_ = F0*F1_*F3_/F2_;

	*F1p = F1_;
	*F2p = F2_;
	*F3p = F3_;
	*F4p = F4_;
}
//==============================================================================//
//
//==============================================================================//
void CorrespondenceData::calculateG( double *G1p, double *G2p, double *G3p, double *G4p, double G0, double G1, double GG1, double GG2 )
{
	double G1_ =GG1;
	double G2_ =GG2/G1_-G1_;
	double G3_ =G2_/(G1_*G0-1);
	double G4_ =G3_*G3_*G1_*G1/G2_;

	*G1p = G1_;
	*G2p = G2_;
	*G3p = G3_;
	*G4p = G4_;
}
//==============================================================================//
// INPUT VALUES:
//
// Rkm - raduis from COE to observer
// Dkm - distance to satellite
// THETAradians - observed elevation
// pMillibar - pressure at surface
// tKelvin - temperature at surface
//
// OUTPUT:
// returns refraction angle in radians
//==============================================================================//
double CorrespondenceData::dE_Theta( double pMillibar, double tKelvin, double Rkm, double Dkm, double THETAradians )
{
	double HH,P,Q,i0,i1,II1,II2,k0;
	double N0,R0,THETA0,RR; 
	//N0 - refractivity at surface, 
	//R0 - radius to observer, 
	//RR - staright line distance, 
	//THETA0 - angle of arrival

	double F0,F1,F1_,F2_,F3_,F4_,FF1,FF2;
	double I,I1_,I2_,I3_,I4_;
	double SN,CN;
	double LL;
	double DEE;

	R0 = Rkm;
	RR = Dkm;
	THETA0 = THETAradians;
	N0 = refractivity(pMillibar, tKelvin);


	HH = 1/log(N0/(N0-7.32*exp(0.005577*N0)));  //km
	P = sqrt(2*HH/R0);
	Q = 1e-6*N0*R0/HH;
	i0 = sqrt(PI)*pow(1-0.9206*Q,-0.4468);
	i1 = 2/(1-Q);
	II1 = 0.5*(1-0.5*Q);
	II2 = 0.75*(1-0.75*Q+(1/6)*Q*Q);
	k0 = sqrt(2*PI)*pow(1-0.9408*Q,-0.4759);

	F0 = i0;
	F1 = i1;
	FF1 = II1;
	FF2 = II2;

	calculateF( &F1_, &F2_, &F3_, &F4_, F0, F1, FF1, FF2 );

	I1_ = P*P*F1_;
	I2_ = P*P*F2_;
	I3_ = P*P*F3_;
	I4_ = P*F4_;

	SN = sin(THETA0);
	CN = cos(THETA0);

	I = 1/(SN+I1_/(SN+I2_/(SN+I3_/(SN+I4_))));

	LL = 1-I*SN+0.5e-6*N0*I*I;

	DEE = 1e-6*N0*CN*(I-R0*LL/RR); //Radians

	return DEE;
}
//==============================================================================//
// INPUT VALUES:
//
// Rkm - raduis from COE to observer
// Dkm - distance to satellite
// Eradians - true elevation
// pMillibar - pressure at surface
// tKelvin - temperature at surface
//
// OUTPUT:
// returns refraction angle in radians
//
//==============================================================================//
double CorrespondenceData::dE_Elev( double pMillibar, double tKelvin, double Rkm, double Dkm, double Eradians )
{
	double HH,P,Q;
	double N0,R0,E,RR;
	//N0 - refractivity at surface, 
	//R0 - radius to observer, 
	//RR - staright line distance, 
	//E - elevation
	double F0,F1,F1_,F2_,F3_,F4_,FF1,FF2;
	double G0,G1,G1_,G2_,G3_,G4_,GG1,GG2;
	double U,U1_,U2_,U3_,U4_;
	double U0,U1,U2,UU1,UU2;
	double UP, UP1_, UP2_, UP3_, UP4_;
	double SN,CN;
	double IPU,DEN,IPPU,PAREN;
	double DEE;

	R0 = Rkm;
	RR = Dkm;
	N0 = refractivity(pMillibar, tKelvin);
	E = Eradians;

	HH = 1/log(N0/(N0-7.32*exp(0.005577*N0)));  //km
	P = sqrt(2*HH/R0);
	Q = 1e-6*N0*R0/HH;
	IPU = -2*pow(1+1.482*Q,-0.3826);
	DEN = 1-0.5*Q*IPU;
	IPPU = 2*sqrt(PI)*pow(1+1.71*Q,0.1);
	U0 = sqrt(PI)*pow(1+1.4844*Q,-0.39144);
	U1 = -IPU/DEN;
	U2 = IPPU/pow(DEN,3);
	UU1 = 0.5*(1+0.5*Q);
	UU2 = 0.75*(1+7/12*Q+1/6*Q*Q);

	F0 = U0;
	F1 = U1;
	FF1 = UU1;
	FF2 = UU2;

	calculateF( &F1_, &F2_, &F3_, &F4_, F0, F1, FF1, FF2 );

	U1_ = P*P*F1_;
	U2_ = P*P*F2_;
	U3_ = P*P*F3_;
	U4_ = P*F4_;

	G0 = U1;
	G1 = U2;
	GG1 = 2*UU1;
	GG2 = 5*UU2;
	calculateG( &G1_, &G2_, &G3_, &G4_, G0, G1, GG1, GG2 );

	UP1_ = P*P*G1_;
	UP2_ = P*P*G2_;
	UP3_ = P*P*G3_;
	UP4_ = P*G4_;

	SN = sin(E);
	CN = cos(E);

	U = 1/(SN+U1_/(SN+U2_/(SN+U3_/(SN+U4_))));
	UP = -(SN*SN+UP1_)/(1+UP2_/(UP3_+SN*(UP4_+SN)));
	PAREN = 1-U*(SN-0.5e-6*N0*U);
	DEE = 1e-6*N0*CN*(U-R0/RR*PAREN*(1+1e-6*N0*UP));

	return DEE;
}
//==============================================================================//
// поправка на рефракцию
//==============================================================================//
void CorrespondenceData::refractcorrect( double optRa, double optDec, double *Telicrf, double *outRa, double *outDec, double SatDist, double rDist )
{
	double pi = 3.1415926535;
	double az = GetZenitAngle( optRa, optDec, Telicrf );

	//////////////////////////////////////////////////////////////////////
	// вычисление рефракции
	// Rkm - raduis from COE to observer
	// Dkm - distance to satellite
	// THETAradians - observed elevation
	// pMillibar - pressure at surface
	// tKelvin - temperature at surface
	// OUTPUT:
	// returns refraction angle in radians
	//function dE_Theta(pMillibar, tKelvin,Rkm,Dkm,THETAradians : Double) : Double;
	//////////////////////////////////////////////////////////////////////

	// K
	double pMillibar = 1.01325*1000.0*594.3/760.0;
	//double tKelvin = 273 - 3.0;	// T kelvin
	double tKelvin = 273.0 + 10.0;	// T kelvin
	double Rkm = 6371.0 + 2.0850013;

	// M
	//double pMillibar = 1.01325*1000.0*750.3/760.0;
	//double tKelvin = 273 + 10.0;	// T kelvin
	//double Rkm = 6371.0;

	//
	//double pMillibar = 1.01325*1000.0*755 / 660.0;
	//double tKelvin = 273 + 20.0;	// T kelvin
	//double Rkm = 6371.0 + 1.576;


	double Dkm = SatDist;
	double THETAradians = pi/2.0-az;

	double resRef = dE_Theta( pMillibar, tKelvin, Rkm, Dkm, THETAradians );
	double resRefStar = dE_Theta( pMillibar, tKelvin, Rkm, 10000000000000.0, THETAradians );
	double resDelta = resRefStar-resRef;
	//printf("p = %f\t ps = %f\tDist = %f\t zenit %f\t E = %f\n", resRef/pi*180.0*3600.0, resRefStar/pi*180.0*3600.0, SatDist, az/pi*180.0, resDelta/pi*180.0*3600.0 );

	double b = resDelta;
	//printf("B = %f\n", b/pi*180.0*3600.0 );

	// вектор телескопа
	S3DCoordinate Rtel;
	Rtel.x = Telicrf[0];
	Rtel.y = Telicrf[1];
	Rtel.z = Telicrf[2];

	// вектор наблюдения
	S3DCoordinate Rn;
	Rn.x = cos(optDec)*cos(optRa);
	Rn.y = cos(optDec)*sin(optRa);
	Rn.z = sin(optDec);

	// вектор вокруг которго надо поворачивать
	S3DCoordinate Rrot;
	Rrot = Rn^Rtel;
	Rrot = Rrot/Rrot.norm();

	double X = Rrot.x;
	double Y = Rrot.y;
	double Z = Rrot.z;

	S3DCoordinate Rnc;
	Rnc.x = (cos(b) +(1.0 - cos(b))*X*X)*Rn.x		+ ((1.0 - cos(b))*X*Y - Z*sin(b))*Rn.y	+ ((1.0 - cos(b))*X*Z + Y*sin(b))*Rn.z;
	Rnc.y = ((1.0 - cos(b))*Y*X + Z*sin(b))*Rn.x	+ (cos(b) +(1.0 - cos(b))*Y*Y)*Rn.y		+ ((1.0 - cos(b))*Y*Z - X*sin(b))*Rn.z;
	Rnc.z = ((1.0 - cos(b))*Z*X - Y*sin(b))*Rn.x	+ ((1.0 - cos(b))*Z*Y + X*sin(b))*Rn.y	+ (cos(b) +(1.0 - cos(b))*Z*Z)*Rn.z;

	double cosa = (Rtel*Rnc )/( Rtel.norm()*Rnc.norm() );
	double zenitc = acos( cosa );
	double deltaz = (zenitc-az)/pi*180.0*3600.0;
	//printf( "correct zenit Angle = %f, In Grad = %f deltaz = %f\n", zenitc, zenitc/pi*180.0, deltaz );

	double crr = atan2( Rnc.y, Rnc.x );
	if( crr < 0 )
		crr = 2.0*pi + crr;

	double crd = atan2( Rnc.z, sqrt( Rnc.x*Rnc.x + Rnc.y*Rnc.y ) );
	//printf("%f %f\n", crr/pi*180.0, crd/pi*180.0 );

	*outRa = crr;
	*outDec = crd;
}
//==============================================================================//
//
//==============================================================================//
void CorrespondenceData::TestCorrectRefract()
{
	double pi = 3.1415926535;

	printf("Test Correct Refraction, Dkm = 1000, Normal P and T\n");
	for( int it = 0; it < 19; it++ )
	{
		double az = ((double)it)*5.0/180.0*pi;

		//////////////////////////////////////////////////////////////////////
		// вычисление рефракции
		// Rkm - raduis from COE to observer
		// Dkm - distance to satellite
		// THETAradians - observed elevation
		// pMillibar - pressure at surface
		// tKelvin - temperature at surface
		// OUTPUT:
		// returns refraction angle in radians
		//function dE_Theta(pMillibar, tKelvin,Rkm,Dkm,THETAradians : Double) : Double;
		//////////////////////////////////////////////////////////////////////
		double pMillibar = 1.01325*1000.0;
		double tKelvin = 283;	// T kelvin
		double Rkm = 6371.0;
		double Dkm = 1000;
		double THETAradians = pi/2.0-az;

		double resRef = dE_Theta( pMillibar, tKelvin, Rkm, Dkm, THETAradians );
		double resRefStar = dE_Theta( pMillibar, tKelvin, Rkm, 10000000000000.0, THETAradians );
		double resDelta = resRefStar-resRef;

		printf( "Angle = %f\t C = %f\n", az/pi*180.0, resDelta/pi*180.0*3600.0 );
	}
	return;
}
//==============================================================================//
// формула Карского
//==============================================================================//
// http://rp5.ru/Погода_в_Шаджатмазе
//double Tatm = -5.0;	// в градусах Цельсия
//double Patm = 594.3;	// мм рт. ст
//double nh = 1.0 + 0.0002926*( Patm/760.0 )*( 273.0 /(Tatm + 273.0));
//printf("nh = %f\n", nh );
// moscow
/*on_surface location;
location.height = 190.0; //190m
location.latitude = 55.80364; //deg
location.longitude = 37.54794; //deg
double R = refract (&location, 1, az/pi*180.0 ); //deg

// constants
// учет рефракции
double nh0 = 1.00029255;
double p1 = 206264.8;
double z1 = az; //rad

double H = rDist - 6371.0; //km
double r = 6371.0;	//km
double Rr = R/180.0*pi; // rad
double Rs = R*3600.0; //sec

// показатель преломления
//double nh = 1.0 + 0.0002926*( 594.0/760.0 );

// формула для вычисления поправки
// стр 172
double PR = (nh0 - 1.0)*p1*sin(z1) - Rs*cos(z1) + 0.5*sin(z1)*sin(z1)*Rs*Rs/p1;
// sigma
double sigma = (r+H)/r;
sigma = sigma*sigma - 1;
// корень в знаменателе
double ch1 = cos( z1 + Rr );
ch1 = ch1*ch1 + sigma;
ch1 = sqrt( ch1 );
PR = PR/( ch1 - cos(z1+Rr) ); //sec
// корректируем измерения
// угол поворота, против часовой стрелки, вверх
double b = PR/3600.0/180.0*pi;

// test
//printf("az = %f\t p = %f\t PR = %f\n", az/pi*180.0, R*3600.0, PR );
//!!! подгон
b = b*0.82;
//printf("B = %f\n", b/pi*180.0*3600.0 );
*/
//==============================================================================//
