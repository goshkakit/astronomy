//==============================================================================//
// Andrianov N.G.
// opbit predict 
// module find Influence Force
// find force from Sun
//==============================================================================//
#include <math.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

//==============================================================================//
// расчет теневой функции, или иначе, процент закрываемого диска Солнца
//==============================================================================//
__device__ double kernalshadow_geo( double *r_ka, double *r_sun )
{
	double sh;
	double Re = 6.37813660;
	double Rs = 695.5;
	double pi = 3.141592653589793;

	double r_s[3];
	double r_e[3];
	double phi_s, phi_e, phi_se, A, x, y, r_s_abs, r_e_abs;

	sh = 0.0;
	// вектор до солнца
	r_s[0] = r_sun[0] - r_ka[0];
	r_s[1] = r_sun[1] - r_ka[1];
	r_s[2] = r_sun[2] - r_ka[2];
	// вектор до земли
	r_e[0] = -r_ka[0];
	r_e[1] = -r_ka[1];
	r_e[2] = -r_ka[2];

	r_s_abs = sqrt( r_s[0]*r_s[0] + r_s[1]*r_s[1] + r_s[2]*r_s[2] );
	r_e_abs = sqrt( r_e[0]*r_e[0] + r_e[1]*r_e[1] + r_e[2]*r_e[2] );

	double r_se = r_s[0]*r_e[0] + r_s[1]*r_e[1] + r_s[2]*r_e[2];
	// угол между вектором на солнце и вектором в центр земли
	phi_se =acos(r_se/r_s_abs/r_e_abs);

	phi_s = asin(Rs/r_s_abs);
	phi_e = asin(Re/r_e_abs);

	if (phi_e > phi_se)
	{
		sh = 1.0;
	}
	else if( phi_e+phi_s > phi_se ) 
	{
		x = (phi_se*phi_se+phi_s*phi_s-phi_e*phi_e)/2.0/phi_se;
		y = sqrt(phi_s*phi_s-x*x);
		A = phi_s*phi_s*acos(x/phi_s)+phi_e*phi_e*acos((phi_se-x)/phi_e)-phi_se*y;
		sh = A/pi/phi_s/phi_s;
	}
	return sh;
}
//==============================================================================//
// Процедура расчета ускорения от светового давления по самой простой
// модели, в которой площадь сечения аппарата считается постоянной,
// результирующее ускорение направлено от Солнца к аппарату, коэффи-
// циент sp_q выражается в долях массы Солнца.
//==============================================================================//
__device__ void kernalsp_cannonball( double *x, double *f_sp, double *pln_coords11, double mu_plan11, double sp_q )
{
	double epsilon = 1.0E-6;
	double r_sun_sc[3];
	double r_abs;

	// Подновление координат планет, если надо
	//if( abs( t - t_planets_update ) > epsilon)
	//	planets_update_geo(t);

	r_sun_sc[0] = x[0] - pln_coords11[0];
	r_sun_sc[1] = x[1] - pln_coords11[1];
	r_sun_sc[2] = x[2] - pln_coords11[2];

	r_abs = sqrt(  r_sun_sc[0]*r_sun_sc[0] +  r_sun_sc[1]*r_sun_sc[1] +  r_sun_sc[2]*r_sun_sc[2] );

	// сила давления
	// f_sp = r_sun_sc*sp_q*mu_plan(11)/r_abs/r_abs/r_abs
	double koeff = sp_q*mu_plan11/(r_abs*r_abs*r_abs);
	f_sp[0] = r_sun_sc[0]*koeff;
	f_sp[1] = r_sun_sc[1]*koeff;
	f_sp[2] = r_sun_sc[2]*koeff;
}
//==============================================================================//
// давление от солнца
// требует предварительного обновления положения плант
//==============================================================================//
__device__ void kernalsp_cannonballForce( double *x, double sp_q, double *f_sp, double *pln_coords11 )
{
	//double pln_coords11[3];
	//pln_coords( 10, pln_coords11 );

	double mu_plan11;
	mu_plan11 = 132712439.935;//mu_plan(10);

	// add shadow
	double sh = kernalshadow_geo( x, pln_coords11 );
	sp_q  = sp_q *(1.0-sh);

	kernalsp_cannonball( x, f_sp, pln_coords11, mu_plan11, sp_q );
}
//==============================================================================//
