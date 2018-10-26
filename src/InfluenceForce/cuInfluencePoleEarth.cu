//==============================================================================//
// Andrianov N.G.
// opbit predict 
// module find Influence Force
// Motion Pole Earth
// GPU version
//==============================================================================//
#include <stdio.h>
#include <math.h>

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

//==============================================================================//
// ѕроцедура получени€ значени€ координат полюса «емли и сдвига по вре-
// мени в момент t.
//==============================================================================//
__device__ void kernalget_xyt ( double t, double *xyt, double ajd0, double delt0, double* cufinals_tab, int cufinals_n )
{
	double jd0000 = 2400000.0;
	double xyt4[4];
	double jd, jd00, res;
	int k, k5;

	xyt[0] = 0.0;
	xyt[1] = 0.0;
	xyt[2] = 0.0;

	jd   = ajd0+(delt0+t)/86.40;
	jd00 = cufinals_tab[1] + jd0000;

	if ((jd < 23.0E+5) || (jd > 25.0E+5)) 
	{
		printf( "Module IERS -> get_xyt : unexpected Julian date. Return\n" );
		return;
	}
	k = int(jd-jd00)+1;
	if (k <= 0) 
	{
		printf( "jd = %f jd00 = %f\n", jd, jd00 );
		printf( "Module IERS -> get_xtr : Time is less than first record of finals.dat table. Return\n" );
	}
	if (k > cufinals_n-2) 
	{
		k = cufinals_n-2;
	}
	else if (k < 1) 
	{
		k = 1;
	}
	res = jd - ( jd00 + (double)(k-1) );
	//res = jd - ( jd00 + (double)(k) );
	k5 = (k-1)*5;
	//k5 = (k)*5;

	//xyt4 = finals_tab(k5+2:k5+5)*(1.d0-res)+finals_tab(k5+7:k5+10)*res
	xyt4[0] = cufinals_tab[k5+2]*(1.0-res) + cufinals_tab[k5+7]*res;
	xyt4[1] = cufinals_tab[k5+3]*(1.0-res) + cufinals_tab[k5+8]*res;
	xyt4[2] = cufinals_tab[k5+4]*(1.0-res) + cufinals_tab[k5+9]*res;
	xyt4[3] = cufinals_tab[k5+5]*(1.0-res) + cufinals_tab[k5+10]*res;

	//xyt = xyt4([1, 2, 4])
	xyt[0] = xyt4[0];
	xyt[1] = xyt4[1];
	xyt[2] = xyt4[3];
}
//==============================================================================//
// ѕроцедура расчета матрицы поворота соответствующей смещению полюсов
//==============================================================================//
__device__ void kernalPM_mat ( double x, double y, double *A_pole )
{
	double a, qx, qy, gamma, x2y2;
	qx   = x*x;
	qy   = y*y;
	x2y2 = qx+qy;
	if( (x2y2) < 1.0E-18)
	{
		// —мещение полюсов нулевое или почти нулевое:
		// искома€ матрица - единична€
		for( int it = 0; it < 9; it++ )
			A_pole[it] = 0.0;

		A_pole[0] = 1.0;
		A_pole[4] = 1.0;
		A_pole[8] = 1.0;
	}
	else
	{
		// «адание матрицы, соответствующей смещению полюсов
		a     = sqrt(1.0-x2y2);
		gamma = -x*y*(a-1.0)/(x2y2);
		A_pole[0] = (qy+qx*a)/(x2y2);
		A_pole[4] = (qx+qy*a)/(x2y2);
		A_pole[1] = gamma;
		A_pole[3] = gamma;
		A_pole[8] = a;

		A_pole[2] = x;
		A_pole[6] = -x;
		A_pole[5] = -y;
		A_pole[7] = y;
	}
}
//==============================================================================//