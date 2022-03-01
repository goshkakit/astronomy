//==============================================================================//
// Andrianov N.G.
// opbit predict 
// module find Influence Force
// find force from Earth harmonics
// gpu version
//==============================================================================//
#include <math.h>

#include "cuda.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

//==============================================================================//
// Процедура расчета возмущений от нецентральности гравитационного поля Земли.
// Вектор x передается в связанной с Землей СК.
// Re - earht equatorial radius
//==============================================================================//
__device__ void kernalGetF_Harm_egm96( double *x_b, int n_harm, double *f_harm, double *CEGM96 )
{
	// константы
	double muE = 398.60044150;
	double Re  = 6.3781360;

	double eps = 2.0E-12;
	int n, m, n1, n21, pos_mm, pos_nm1, pos_nm;
	double x[3];
	double ro, rv, w1, w2, w3, w4, dw2, dw3, dw;
	double cosL, sinL, cosmL, sinmL, cosm1L, sinm1L;
	double R, R_n, R_m, tmp, P_nm, P_n1m, P_mm, PR;

	ro = x_b[0]*x_b[0]+x_b[1]*x_b[1];
	rv = sqrt(x_b[2]*x_b[2]+ro);

	x[0] = x_b[0]/rv;
	x[1] = x_b[1]/rv;
	x[2] = x_b[2]/rv;

	f_harm[0] = 0.0;
	f_harm[0] = 0.0;
	f_harm[0] = 0.0;

	ro = sqrt(ro);
	if(ro > eps)
	{
		cosL = x_b[0]/ro;
		sinL = x_b[1]/ro;
	}
	else
	{
		cosL = 1.0;
		sinL = 0.0;
	}

	ro = ro/rv;
	R  = Re/rv;
	w1 = 0.0;
	w2 = 0.0;
	w3 = 0.0;
	w4 = 0.0;
	// Значения cos(ml), sin(ml), cos((m-1)l), sin((m-1)l) для начального значения m = 1
	cosm1L = 1.0;
	sinm1L = 0.0;
	cosmL  = cosL;
	sinmL  = sinL;
	// Позиции в массиве коэффициентов egm96(pos_nm) = C_nm, egm96(pos_nm+1) = S_nm
	// за исключением S_n0 = 0, которые не входят в массив.

	//############ индексы исправить для egm96,
	//############ сейчас фортрановская нумерация

	pos_mm  = 2;
	pos_nm1 = 1;
	P_mm = 1.0;
	R_m  = R;
	//do m = 1,n_harm 
	for( m = 1; m <=n_harm; m++ )
	{
		// Внешний цикл, проход по m от 1 до числа рассматриваемых гармоник
		// Начальные значения присоединенных функций Лежандра
		P_nm   = P_mm;
		P_n1m  = 0.0;
		pos_nm = pos_mm;
		R_n = R_m;
		n1  = m+1;
		dw2 = 0.0;
		dw3 = 0.0;
		// do n = m,n_harm
		for( n = m; n <=n_harm; n++ )
		{
			//! Внутренний цикл, проходи по n от m до числа рассматриваемых гармоник
			PR = P_nm*R_n;
			//dw  = PR*( cuegm96[pos_nm]*cosmL + cuegm96[pos_nm+1]*sinmL); 
			dw  = PR*( CEGM96[pos_nm]*cosmL + CEGM96[pos_nm+1]*sinmL);
			w1  = w1 + ((double)(n1))*dw;
			dw2 = dw2 + dw;
			//dw3 = dw3 + PR*(-cuegm96[pos_nm]*sinmL + cuegm96[pos_nm+1]*cosmL);
			//w4  = w4 + PR*( cuegm96[pos_nm1]*cosm1L + cuegm96[pos_nm1+1]*sinm1L);
			dw3 = dw3 + PR*(-CEGM96[pos_nm]*sinmL + CEGM96[pos_nm+1]*cosmL);
			w4  = w4 + PR*( CEGM96[pos_nm1]*cosm1L + CEGM96[pos_nm1+1]*sinm1L);
			R_n = R_n*R;
			tmp = P_n1m;
			P_n1m = P_nm;
			n21 = n+n1;
			P_nm = ( ((double)(n21))*x[2]*P_n1m - ((double)(n+m))*tmp)/((double)(n1-m));
			n1 = n1+1;
			pos_nm1 = pos_nm1 + n21;
			pos_nm = pos_nm + n21;
		}
		w2 = ((double)(m))*dw2+w2;
		w3 = ((double)(m))*dw3+w3;
		R_m = R_m*R;
		cosm1L = cosmL;
		sinm1L = sinmL;
		cosmL = cosm1L*cosL - sinm1L*sinL;
		sinmL = cosm1L*sinL + sinm1L*cosL;
		n21 = m+m+1;
		P_mm = ((double)(n21))*ro*P_mm;
		pos_nm1 = n21 + pos_mm;
		pos_mm = pos_nm1 + 2;
	}
	//! Суммирование членов с m=0
	P_n1m = 1.0;
	P_nm = x[2];
	pos_mm = 1;
	R_n = R;
	w1 = w1*ro;
	//do n = 1,n_harm
	for( n = 1; n <=n_harm; n++ )
	{
		//w1 = w1 + ((double)(n+1))*cuegm96[pos_mm]*R_n*P_nm;
		w1 = w1 + ((double)(n+1))*CEGM96[pos_mm]*R_n*P_nm;
		if (n != n_harm )
		{
			n1 = n+1;
			R_n = R_n*R;
			tmp = P_n1m;
			P_n1m = P_nm;
			n21 = n+n1;
			P_nm = ( ((double)(n21))*x[2]*P_n1m - ((double)(n))*tmp )/n1;
			pos_mm = pos_mm +n21;
		}
	}
	//! Вычисление действующего ускорения
	tmp = -muE/rv/rv;
	w2 = w4*ro - w2*x[2];
	f_harm[2] = (w1*x[2] - w2*ro)*tmp;
	w2 = w2*x[2];
	f_harm[0] = (w2*cosL + w3*sinL + w1*x[0])*tmp;
	f_harm[1] = (w2*sinL - w3*cosL + w1*x[1])*tmp;
}

//==============================================================================//
// получение ускорения вызванного гармониками
//==============================================================================//
__device__ void kernalGetHarmForce( double *x, double *Fharm, double *cuegm96 )
{
	//int n = 75;
	kernalGetF_Harm_egm96( x, 75, Fharm, cuegm96 );
};
//==============================================================================//