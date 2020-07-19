//==============================================================================//
// Andrianov N.G.
// opbit predict 
// module find Influence Force
// Additional Function
//==============================================================================//
#include <math.h>
#include <stdio.h>

#include "InfluenceForce.h"
#include "common\DataConverter.h"

namespace Force
{
	//==============================================================================//
	// дополнительные модули и функции умножение векторов и матриц
	// mat*mat
	//==============================================================================//
	void InfluenceForce::matMul( double *in1, double *in2, double *out )
	{
		for( int j = 0; j < 3; j++ )
		{
			for( int i = 0; i < 3; i++ )
			{
				out[j*3+i] = 0;
				for( int k = 0; k < 3; k++ )
				{
					out[j*3+i] += in1[j*3+k]*in2[k*3 + i];
				}
			}
		}
	}
	//==============================================================================//
	// mat*vec
	//==============================================================================//
	void InfluenceForce::matVecMul( double *inMat, double *inVec, double *outVec )
	{
		for( int j = 0; j < 3; j++ )
		{
			outVec[j] = 0;
			for( int k = 0; k < 3; k++ )
			{
				outVec[j] += inMat[j*3+k]*inVec[k];
			}
		}
	}
	//==============================================================================//
	// mat*pos and velocity
	//==============================================================================//
	void InfluenceForce::matVecMul_V6( double *inMat, double *inVec, double *outVec )
	{
		for( int j = 0; j < 3; j++ )
		{
			outVec[j] = 0;
			for( int k = 0; k < 3; k++ )
			{
				outVec[j] += inMat[j*3+k]*inVec[k];
			}
		}

		for( int j = 0; j < 3; j++ )
		{
			outVec[j+3] = 0;
			for( int k = 0; k < 3; k++ )
			{
				outVec[j+3] += inMat[j*3+k]*inVec[k+3];
			}
		}
	}
	//==============================================================================//
	// транспонирование матриц
	//==============================================================================//
	void InfluenceForce::transpose( double *A, double *B )
	{
		for( int j = 0; j < 3; j++ )
		{
			for( int i = 0; i < 3; i++ )
			{
				B[i*3+j] = A[j*3+i];
			}
		}
	}
	//==============================================================================//
	// дополнительные функции преобразования чиисле и округления
	//==============================================================================//
	double InfluenceForce::DMOD( double X, double Y )
	{
		int s = (int)(X/Y);
		double res = X - ((double)s)*Y;

		return res;
	}
	double InfluenceForce::DDIM( double X, double Y  )
	{
		double res = X-Y;
		if( res < 0 )
			res = 0;

		return res;
	}
	double InfluenceForce::DSIGN( double X, double Y )
	{
		double sig = 1;
		if( Y < 0 )
			sig = -1;
		if( Y == 0 )
			sig = 0;

		X = sig*X;

		return X;
	}
	double InfluenceForce::Dmax( double X, double Y )
	{
		if( X > Y )
			return X;
		else
			return Y;
	}
	double InfluenceForce::Dmin( double X, double Y )
	{
		if( X < Y )
			return X;
		else
			return Y;
	}
	//==============================================================================//


	//==============================================================================//
	// перевод спрогнозированного положения в координаты RA DEC
	//==============================================================================//
	void InfluenceForce::ConvertXYZtoRADEC(double *posICRF, double *TelICRF, double *Ra, double *Dec)
	{
		// входные данные задаются в ICRF
		// вектор направления в системе ICRF
		double x = posICRF[0] - TelICRF[0];
		double y = posICRF[1] - TelICRF[1];
		double z = posICRF[2] - TelICRF[2];

		double r = atan2(y, x);
		double d = atan2(z, sqrt(x*x + y*y));

		double pi = 3.1415926535;

		if (r < 0)
			r = 2.0*pi + r;

		*Ra = r;
		*Dec = d;
	}

	// jd full date
	void InfluenceForce::ITRFToICRF(double jd, double *posITRF, double *posICRF)
	{
		DataConverter Dconv;
		// MDB
		jd = jd + 0.125;
		double dataMDB = Dconv.JDtoYYYYMMDD(jd);
		double timeMDB = Dconv.SECtoHHMMSS(dataMDB, jd);

		// установка времени
		double int1, ajd1, delt1;
		set_time(dataMDB, timeMDB, &ajd1, &delt1, &int1);

		// матрица перевода в земную систему
		double Arot[9];
		iers_update_matrix(int1, Arot, ajd1, delt1);
		// матрица перехода из земной в нормальную систему
		double invArot[9];
		transpose(Arot, invArot);

		matVecMul(invArot, posITRF, posICRF);
	}

	void InfluenceForce::ICRFToITRF(double jd, double *posICRF, double *posITRF)
	{
		DataConverter Dconv;
		// MDB
		jd = jd + 0.125;
		double dataMDB = Dconv.JDtoYYYYMMDD(jd);
		double timeMDB = Dconv.SECtoHHMMSS(dataMDB, jd);

		// установка времени
		double int1, ajd1, delt1;
		set_time(dataMDB, timeMDB, &ajd1, &delt1, &int1);

		// матрица перевода в земную систему
		double Arot[9];
		iers_update_matrix(int1, Arot, ajd1, delt1);

		matVecMul(Arot, posICRF, posITRF);
	}
};