//==============================================================================//
// Andrianov N.G.
// opbit predict 
// module find Influence Force
// Nutation Earth
// GPU version
//==============================================================================//

#include "cuda_runtime.h"
#include "device_launch_parameters.h"


#include "cuInfluiencePrecessionEarth.cu"
#include "cuInfluencePoleEarth.cu"
#include "cuInfluenceNutationEarth.cu"
//==============================================================================//
// Процедура перевода вектора состояния из EME2000 в гринвичскую вра щающуюся СК.
// процедура, переводящая целиком вектор состояния в систему, связанную с землей
// если ты ее переделал на C, то все ок
//==============================================================================//
__device__ void kernalstate_to_itrf( double t, double *x_2000, double *x_g, double *A_rot )
{
	double omega_earth = 0.072921150850;
	// обновление матрицы вращения
	//double t_epsilon = 1.0E-6;
	//if (abs(t-t_iers_update).GT.t_epsilon) 
	//	iers_update_matrix(t)

	//matVecMul( A_rot, x_2000, x_g );
	kernalmatVecMul_V6( A_rot, x_2000, x_g );

	x_g[3] = x_g[3] + omega_earth*x_g[1];
	x_g[4] = x_g[4] - omega_earth*x_g[0];
};
//==============================================================================//
// Процедура получения матрицы соотвутствующей суточному вращению Земли.
//==============================================================================//
__device__ void kernalER_mat ( double t, double d_UT1, double *A_rotat, double ajd0, double delt0, double *cuARG, double *cuAMPL )
{
	double E0 = 2451545.0;
	double THJ = 36525.0;
	double pi2 = 6.2831853071795860;
	double hyt[2];
	double ID, JD, t100, s, a, q, z, dUT1, D;

	dUT1=d_UT1/86400.0;
	ID = ajd0 + (delt0+t)/86.40 - E0 - dUT1;

	double tmpv = (delt0+t)/86.40;

	D = kernalDMOD(ajd0, 1.0) + kernalDMOD(tmpv, 1.0)-dUT1 + 0.50;
	t100 = ID/THJ;

	s = 24110.548410+ID*236.5553679080+D*86400.0 + t100*t100*(0.0931040-t100*6.2E-6);

	JD = ajd0+(delt0+t)/86.40 + dUT1;

	a = kernalE2000(JD, JD);
	kernalN2000( 106, JD, hyt, cuARG, cuAMPL );

	q = s/86400.0+hyt[0]*cos(a)/pi2;

	a = int(q);
	z = (q-a)*pi2;

	for( int it = 0; it < 9; it++ )
		A_rotat[it] = 0.0;

	A_rotat[0] = cos(z);
	A_rotat[1] = sin(z);
	A_rotat[3] = -sin(z);
	A_rotat[4] = cos(z);
	A_rotat[8] = 1.0;
}
//==============================================================================//
// Процедура получения матрицы перехода из ITRF в ICRF. Иными словами, 
// матрица перехода от гринвичской вращающейся системы координат, 
// зафиксированной на момент времени t, в систему EME2000.
//==============================================================================//
__device__ void kernaliers_mat( double t, double *A, double ajd0, double delt0, double *cuARG, double *cuAMPL, double* cufinals_tab, int cufinals_n )
{
	double jd2k = 2451545.0;
	double xyt[3];
	double A_pole[9];
	double A_rotat[9];
	double A_prc[9];
	double A_nut[9];
	double jd;

	jd = ajd0+(delt0+t)/86.40;

	// прецессия
	kernalPM2000(jd2k, jd, A_prc);
	// нутация
	kernalNM2000(jd, A_nut, cuARG, cuAMPL );
	// движение полюсов
	kernalget_xyt( t, xyt, ajd0, delt0, cufinals_tab, cufinals_n );
	kernalPM_mat(xyt[0], xyt[1], A_pole);
	// суточное вращение
	kernalER_mat(t, xyt[2], A_rotat, ajd0, delt0, cuARG, cuAMPL );

	// результирующая матрица вращения
	double T1[9], T2[9];
	kernalmatMul( A_pole, A_rotat, T1 );
	kernalmatMul( T1 ,A_nut, T2 );
	kernalmatMul( T2 ,A_prc, A );
}

//==============================================================================//
// Процедура обновления текущего значения матрицы перехода из Гринвич-
// ской вращающейся СК в EME2000. При этом матрица поворота сохраняется
// в модульной переменной A_gr.
//==============================================================================//
//__device__ void kernaliers_update_matrix( double t, double *A_rot, double ajd0, double delt0 )
//{
//	kernaliers_mat(t, A_rot, ajd0, delt0 );
//	double t_iers_update = t;
//}
//==============================================================================//