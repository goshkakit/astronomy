//==============================================================================//
// Andrianov N.G.
// opbit predict 
// module find Influence Force
// Additional Function
//==============================================================================//
#include <math.h>
#include <stdio.h>

#include "InfluenceForce.h"

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
};