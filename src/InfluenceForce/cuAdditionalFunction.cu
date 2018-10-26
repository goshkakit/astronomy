//==============================================================================//
// дополнительные модули и функции умножение векторов и матриц
// mat*mat
//==============================================================================//
__device__ void kernalmatMul( double *in1, double *in2, double *out )
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
__device__ void kernalmatVecMul( double *inMat, double *inVec, double *outVec )
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
__device__ void kernalmatVecMul_V6( double *inMat, double *inVec, double *outVec )
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
__device__ void kernaltranspose( double *A, double *B )
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
__device__ double kernalDMOD( double X, double Y )
{
	int s = (int)(X/Y);
	double res = X - ((double)s)*Y;

	return res;
}

__device__ double kernalDDIM( double X, double Y  )
{
	double res = X-Y;
	if( res < 0 )
		res = 0;

	return res;
}

__device__ double kernalDSIGN( double X, double Y )
{
	double sig = 1;
	if( Y < 0 )
		sig = -1;
	if( Y == 0 )
		sig = 0;

	X = sig*X;

	return X;
}
__device__ double kernalDmax( double X, double Y )
{
	if( X > Y )
		return X;
	else
		return Y;
}
__device__ double kernalDmin( double X, double Y )
{
	if( X < Y )
		return X;
	else
		return Y;
}
//==============================================================================//