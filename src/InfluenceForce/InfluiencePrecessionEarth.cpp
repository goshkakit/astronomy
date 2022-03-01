//==============================================================================//
// Andrianov N.G.
// opbit predict 
// module find Influence Force
// Precession Earth
//==============================================================================//
#include <math.h>
#include <stdio.h>

#include "InfluenceForce.h"

namespace Force
{
	//==============================================================================//
	// расчет матрицы прицессии  
	//==============================================================================//
	void InfluenceForce::PM2000( double E1, double E2, double *P )
	{
		double T1, T2, DT, T12, DT2, DT3;
		double A, B, X, Y, Z, SX;
		double CX, SY, CY, SZ, CZ;

		double E0 = 2451545.0;
		double THJ = 36525.0;           
		double S = 206264.806;

		T1=(E1-E0)/THJ;                                     
		T2=(E2-E0)/THJ;                                     
		DT=T2-T1;                                           
		T12=T1*T1;                                          
		DT2=DT*DT;                                          
		DT3=DT*DT2;                                         
		A = (2306.2181+1.39656*T1-0.000139*T12)*DT;        
		B = (2004.3109-0.85330*T1-0.000217*T12)*DT; 

		X = A+(0.30188-0.000345*T1)*DT2+0.017998*DT3;       
		Y = A+(1.09468+0.000066*T1)*DT2+0.018203*DT3;       
		Z = B+(-0.42665-0.000217*T1)*DT2-0.041833*DT3; 

		X = X/S;                                              
		Y = Y/S;                                              
		Z = Z/S;

		SX = sin(X);                                        
		CX = cos(X);                                         
		SY = sin(Y);                                         
		CY = cos(Y);                                         
		SZ = sin(Z);                                         
		CZ = cos(Z);   

		P[0] = CY*CZ*CX-SY*SX;                               
		P[3] = SY*CZ*CX+CY*SX;                               
		P[6] = SZ*CX;                                        
		P[1] =-CY*CZ*SX-SY*CX;                               
		P[4] =-SY*CZ*SX+CY*CX;                               
		P[7] =-SZ*SX;                                        
		P[2] =-CY*SZ;                                        
		P[5] =-SY*SZ;
		P[8] = CZ;                                                                                               
	} 
	//==============================================================================//
};