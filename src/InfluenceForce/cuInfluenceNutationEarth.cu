//==============================================================================//
// Andrianov N.G.
// opbit predict 
// module find Influence Force
// Nutation Earth GPU
//==============================================================================//
#include <math.h>
#include <stdio.h>

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "cuAdditionalFunction.cu"
//==============================================================================//
// расчет среднего наклона эклиптики
//==============================================================================//
__device__ double kernalE2000( double E1, double E2)
{
	double E0, THJ, S, T1, T2, DT, T12, T13, DT2, DT3;
	E0 = 2451545.0;
	THJ = 36525.0;
	S = 206264.806;                                                       
	T1=(E1-E0)/THJ;                                     
	T2=(E2-E0)/THJ;                                     
	DT=T2-T1;                                           
	T12=T1*T1;
	T13=T12*T1;                                          
	DT2=DT*DT;                                          
	DT3=DT*DT2;
	//NAKLON  EkLIPTIKI k EKBATORY  
	double rE2000 =  ( 84381.4480 - 46.81500*T1 - 0.000590*T12 + 0.0018130*T13)
		+(-46.81500 - 0.001170*T1 + 0.0054390*T12)*DT
		+(-0.000590 + 0.0054390*T1)*DT2+0.0018130*DT3;
	rE2000=rE2000/S; 
	return rE2000;                                            
} 
//==============================================================================//
//
//==============================================================================//
__device__ void kernalFA2000( double AED, double *FA )
{
	double E0, THJ, PII, PI, T, DT2, DT3, A, B, C, P, R;
	E0 = 2451545.0;
	THJ = 36525.0;

	PII=6.283185307179586;
	PI=3.141592653589793;
	T=(AED-E0)/THJ;
	DT2=T*T;
	DT3=T*DT2;

	A=PI*2650.0+3.470890870;
	FA[0]=2.3555483930+A*T+0.1517952E-3*DT2+0.3103E-6*DT3;

	B=PI*198.0+6.26661061;
	FA[1]=6.24003594+B*T-0.27974E-5*DT2-0.582E-7*DT3;

	C=PI*2684.0+1.431476080;
	FA[2]=1.627901930+C*T-0.642717E-4*DT2+0.533E-7*DT3;

	P=PI*2472.0+5.36010650;
	FA[3]=5.198469510+P*T-0.34085E-4*DT2+0.921E-7*DT3;

	R=PI*10.0+2.34111940;
	FA[4]=2.182438620-R*T+0.361429E-4*DT2+0.388E-7*DT3;

	for( int it = 0; it < 5; it++ )
		FA[it]= kernalDMOD( FA[it], PII ) + PII*kernalDDIM(- kernalDSIGN( 1.0, FA[it] ), 0.0 );

}
//==============================================================================//
// вычисление поправки нутации
//==============================================================================//
__device__ void kernalN2000( int N, double AJD, double *HYT, double *cuARG, double *cuAMPL )
{
	double FA[5];
	double B[2];

	double EO, THJ, S, T, A;

	EO = 2451545.0;
	THJ = 36525.0;
	S = 206264.806;

	T = (AJD - EO)/THJ;

	kernalFA2000( AJD, FA );

	HYT[0] = 0.0;
	HYT[1] = 0.0;
	int L = 0;
	for( int i=0; i < N; i++ )
	{
		A = ( cuARG[L]*FA[0] + cuARG[L+1]*FA[1] + cuARG[L+2]*FA[2] + cuARG[L+3]*FA[3] + cuARG[L+4]*FA[4] );

		HYT[0] = ( cuAMPL[L+1] + cuAMPL[L+2]*T)*sin(A) + HYT[0];
		HYT[1] = ( cuAMPL[L+3] + cuAMPL[L+4]*T)*cos(A) + HYT[1];
		B[0] = HYT[0]/S;
		B[1] = HYT[1]/S;
		// PRINT 2,(B(K),K=1,2)
		L=L+5;
	}

	// HYTAWIJ B DOLGOTE B PAD.
	HYT[0]=HYT[0]/S;
	// HYTAWIJ B HAKLOHE B PAD.
	HYT[1]=HYT[1]/S;
	// 1    FORMAT (1X,'BFA200',5f18.11)
	// 2    FORmat (1X,2f20.11)
	// 3    format (1x,e10.3)
}
//==============================================================================//
// общая матрица нутации
//==============================================================================//
__device__ void kernalNM2000( double E, double *HUT, double *cuARG, double *cuAMPL )
{
	double A, Epsi, COS1, COS2, SIN1, SIN2, SIN3;
	int N;
	double HYT[2];

	//INTERFACE 
	//      FUNCTION E2000(E1,E2)
	//      REAL(KIND=8) :: E1
	//      REAL(KIND=8) :: E2
	//      REAL(KIND=8) :: E2000
	//      END FUNCTION E2000
	//   END INTERFACE 

	N = 106;      
	A = kernalE2000( E, E );
	kernalN2000( N, E, HYT, cuARG, cuAMPL );

	Epsi = A + HYT[1];
	COS1 = cos( A );
	COS2 = cos( Epsi );
	SIN1 = sin( HYT[0] );
	SIN2 = sin( A );
	SIN3 = sin( Epsi );

	//MATPITSA HYTATSII
	HUT[0] = cos( HYT[0] );
	HUT[1] = -SIN1*COS1;
	HUT[2] = -SIN1*SIN2;
	HUT[3] = SIN1*COS2;
	HUT[4] = HUT[0]*COS2*COS1+SIN3*SIN2;
	HUT[5] = HUT[0]*COS2*SIN2-SIN3*COS1;
	HUT[6] = SIN1*SIN3;
	HUT[7] = HUT[0]*SIN3*COS1-COS2*SIN2;
	HUT[8] = HUT[0]*SIN3*SIN2+COS2*COS1;
}
//==============================================================================//
