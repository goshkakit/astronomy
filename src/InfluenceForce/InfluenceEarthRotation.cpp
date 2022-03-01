//==============================================================================//
// Andrianov N.G.
// opbit predict 
// module find Influence Force
// Nutation Earth
//==============================================================================//
#include <math.h>
#include <stdio.h>

#include "InfluenceForce.h"

namespace Force
{
	// константа
	double omega_earth = 0.072921150850;
	//==============================================================================//
	// Процедура перевода вектора состояния из EME2000 в гринвичскую вра щающуюся СК.
	//subroutine state_to_itrf(t, x_2000, x_g)
	//	implicit none
	//	real*8, intent(in) :: t, x_2000(6)
	//	real*8, intent(out) :: x_g(6)
	//	real*8, parameter :: t_epsilon = 1.d-6
	//	if (abs(t-t_iers_update).GT.t_epsilon) then
	//		call iers_update_matrix(t)
	//		end if
	//		x_g(1:3) = matmul(A_rot, x_2000(1:3))
	//		x_g(4:6) = matmul(A_rot, x_2000(4:6))
	//		x_g(4) = x_g(4)+omega_earth*x_g(2)
	//		x_g(5) = x_g(5)-omega_earth*x_g(1)
	//		end subroutine state_to_itrf
	//		процедура, переводящая целиком вектор состояния в систему, связанную с землей
	//		если ты ее переделал на C, то все ок
	//==============================================================================//
	void InfluenceForce::state_to_itrf( double t, double *x_2000, double *x_g, double *A_rot )
	{
		// обновление матрицы вращения
		//double t_epsilon = 1.0E-6;
		//if (abs(t-t_iers_update).GT.t_epsilon) 
		//	iers_update_matrix(t)

		//matVecMul( A_rot, x_2000, x_g );
		matVecMul_V6( A_rot, x_2000, x_g );

		x_g[3] = x_g[3] + omega_earth*x_g[1];
		x_g[4] = x_g[4] - omega_earth*x_g[0];
	};
	//==============================================================================//
	// Процедура обновления текущего значения матрицы перехода из Гринвич-
	// ской вращающейся СК в EME2000. При этом матрица поворота сохраняется
	// в модульной переменной A_gr.
	//subroutine iers_update_matrix(t)
	//  implicit none
	//  real*8, intent(in) :: t
	//  call iers_mat(t, A_rot)
	//  t_iers_update = t
	//end subroutine iers_update_matrix
	//==============================================================================//
	void InfluenceForce::iers_update_matrix( double t, double *A_rot, double ajd0, double delt0 )
	{
		iers_mat(t, A_rot, ajd0, delt0 );
		double t_iers_update = t;
	}
	//==============================================================================//
	// Процедура получения матрицы соотвутствующей суточному вращению Земли.
	//==============================================================================//
	void InfluenceForce::ER_mat ( double t, double d_UT1, double *A_rotat, double ajd0, double delt0 )
	{
		double E0 = 2451545.0;
		double THJ = 36525.0;
		double pi2 = 6.2831853071795860;
		double hyt[2];
		double ID, JD, t100, s, a, q, z, dUT1, D;

		dUT1=d_UT1/86400.0;
		ID = ajd0 + (delt0+t)/86.40 - E0 - dUT1;

		double tmpv = (delt0+t)/86.40;

		D = DMOD(ajd0, 1.0) + DMOD(tmpv, 1.0)-dUT1 + 0.50;
		t100 = ID/THJ;

		s = 24110.548410+ID*236.5553679080+D*86400.0 + t100*t100*(0.0931040-t100*6.2E-6);

		JD = ajd0+(delt0+t)/86.40 + dUT1;

		a = E2000(JD, JD);
		N2000( 106, JD, hyt );

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
	//
	// ICRF - это средний экватор и среднее равноденствие эпохи j2000
	// NORAD работает и дает элементы орбиты в системе True equatorm, mean equinox
	// т.е. истиный экватор и среднее равноденствие
	// причем эта система уникальна для каждого момента времени 
	//
	// J2000
	// One commonly used ECI frame is defined with the Earth's Mean Equator and
	// Equinox at 12:00 Terrestrial Time on 1 January 2000. 
	// It can be referred to as J2000 or EME2000. The x-axis is aligned with the mean equinox.
	// The z-axis is aligned with the Earth's spin axis or celestial North Pole. 
	// The y-axis is rotated by 90° East about the celestial equator.
	//==============================================================================//
	void InfluenceForce::iers_mat( double t, double *A, double ajd0, double delt0 )
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
		PM2000(jd2k, jd, A_prc);
		// нутация
		NM2000(jd, A_nut);
		// движение полюсов
		get_xyt( t, xyt, ajd0, delt0 );
		PM_mat(xyt[0], xyt[1], A_pole);
		// суточное вращение
		ER_mat(t, xyt[2], A_rotat, ajd0, delt0 );

		// результирующая матрица вращения
		double T1[9], T2[9];
		matMul( A_pole, A_rotat, T1 );
		matMul( T1 ,A_nut, T2 );
		matMul( T2 ,A_prc, A );
	}
	//==============================================================================//
	// Процедура расчитывающая матрицу перехода из J2000 в TEME (true equa-
	// tor mean equinox)
	// http://en.wikipedia.org/wiki/Earth-centered_inertial
	// Earth-centered inertial (ECI)
	// The ECI frame used for the NORAD two-line elements is sometimes called true equator,
	// mean equinox (TEME) although it does not use the conventional mean equinox.
	//==============================================================================//
	void InfluenceForce::GetTemeMatrix( double t, double *A, double ajd0, double delt0 )
	{
		double A_nut[9];
		double A_prc[9];
		double A_r[9];
		double hyt[2];

		double de, jd;
		jd = (t+delt0)/86.40 + ajd0;

		PM2000( 2451545.0, jd, A_prc );

		NM2000( jd, A_nut );

		N2000(106, jd, hyt);

		de = hyt[0];

		for( int it = 0; it < 9; it++ )
			A_r[it] = 0.0;

		A_r[0] = cos(de);
		A_r[1] = sin(de);
		A_r[3] = -sin(de);
		A_r[4] = cos(de);
		A_r[8] = 1.0;

		double T1[9];
		matMul( A_nut, A_prc, T1 );
		matMul( A_r, T1, A );
	}
	//==============================================================================//
};