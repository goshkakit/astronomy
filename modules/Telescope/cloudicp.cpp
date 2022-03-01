//=======================================================================//
// Andrianov N.G.
// icp for two cloud
//=======================================================================//
#include <conio.h>
#include <stdio.h>
#include "cloudicp.h"

#include "svd/ap.h"
#include "svd/svd.h"

//=======================================================================//
// Construct
//=======================================================================//
template <typename T>
MakeIcpForCloud<T>::MakeIcpForCloud()
{
}

template <typename T>
MakeIcpForCloud<T>::~MakeIcpForCloud()
{
}
//=======================================================================//
// load cloud from .3d file
//=======================================================================//
template <typename T>
int MakeIcpForCloud<T>::loadCloud( char *fname, T* arrayPoint )
{
	FILE *file;
	file = fopen( fname, "r" );

	std::vector< S3DCoordinate > tmpPoints;
	std::vector< SRGB > tmpColors;

	//.3d file
	int R, G, B;
	float x, y, z;

	if (file)
	{

		int result = fscanf(file, "%f %f %f %d %d %d\n", &x, &y, &z, &R, &G, &B);
		while(result != EOF)
		{
			tmpPoints.push_back(S3DCoordinate(x, y, z));
			tmpColors.push_back(SRGB(R, G, B));			
			result = fscanf(file, "%f %f %f %d %d %d\n", &x, &y, &z, &R, &G, &B);
		}
		fclose(file);
	}
	else
	{
		printf( " Couldn't open file %s for reading\n", fname);
		return -1;
	}

	int s = tmpPoints.size();
	for( int it = 0; it < s; it++ )
	{
		arrayPoint[3*it+0] = (T)tmpPoints[it].x;
		arrayPoint[3*it+1] = (T)tmpPoints[it].y;
		arrayPoint[3*it+2] = (T)tmpPoints[it].z;
	}
	printf( "Load Size cloud = %d\n", s );
	
	return 0;
};
//=======================================================================//
// поиск ближайших точек на CPU перебором
//=======================================================================//
template <typename T>
T MakeIcpForCloud<T>::FindClosestPoint( unsigned int numSrc, T* PointsSrc, unsigned int numDst, T* PointsDst, T* distansces, unsigned int *CPointIndex )
{
	double TotalDist = 0;
	for( unsigned int it = 0; it < numSrc; it++ )
	{
		// get point
		S3DCoordinate Pm;
		Pm.x = PointsSrc[3*it+0];
		Pm.y = PointsSrc[3*it+1];
		Pm.z = PointsSrc[3*it+2];

		// init result
		unsigned int index = 0;
		double dist = 100000000.0; 

		for( unsigned int k = 0; k < numDst; k++ )
		{
			S3DCoordinate Pc;
			Pc.x = PointsDst[3*k+0];
			Pc.y = PointsDst[3*k+1];
			Pc.z = PointsDst[3*k+2];

			double d = (Pm-Pc).norm();
			if( d < dist )
			{
				dist = d;
				index = k;
			}
		}
		CPointIndex[it] = index;
		distansces[it] = dist;

		TotalDist += dist;
	}

	TotalDist = TotalDist/numSrc;
	printf("Dist = %f\n", TotalDist );
	return TotalDist;
};
//=======================================================================//
// вычисление расстояния, сдвига и поворота
//=======================================================================//
template <typename T>
T MakeIcpForCloud<T>::getShiftRotationFromICP( unsigned int numSrc, T* PointsSrc, unsigned int numDst, T* PointsDst, unsigned int *CPointIndex, S3DCoordinate &totalTn, S3DMatrix &totalRn )
{
	using namespace ap;
	S3DCoordinate fcenter(0,0,0);
	S3DCoordinate ccenter(0,0,0);

	double TotalDist = 0;
	S3DCoordinate Pm;
	S3DCoordinate Pc;
	int closestPointsN = 0;

	for( unsigned int i = 0; i < numSrc; i++)
	{
		unsigned int index = CPointIndex[i];
		if( index >= 0 && index < numDst )
		{		
			Pm.x = PointsSrc[3*i+0];
			Pm.y = PointsSrc[3*i+1];
			Pm.z = PointsSrc[3*i+2];
			
			Pc.x = PointsDst[3*index+0];
			Pc.y = PointsDst[3*index+1];
			Pc.z = PointsDst[3*index+2];
			
			double pdist = (Pm-Pc).norm();

			// глобальное ограничение
			if( pdist > 200 )
				continue;


			closestPointsN++;
			
			TotalDist += pdist;
			fcenter = fcenter + Pm;
			ccenter = ccenter + Pc;
		}
	}

	// ошибка числа точек
	if( closestPointsN == 0 )
	{
		printf("Error\n");
		return -1;
	}
	fcenter = fcenter / closestPointsN;
	ccenter = ccenter / closestPointsN;
	TotalDist = TotalDist/closestPointsN;

	ap::real_2d_array Cov3;
	ap::real_1d_array S;
	ap::real_2d_array U, Vt;

	Cov3.setbounds(1,3,1,3);
	for (int i=1; i<=3; i++)
	{
		for (int j=1; j<=3; j++)
		{
			Cov3(i,j)=0;
		}
	}

	double Yx, Yy, Yz, Xx, Xy, Xz;
	for( unsigned int i = 0; i < numSrc; i++)
	{
		Yx =  PointsSrc[3*i+0] - fcenter.x; 
		Yy =  PointsSrc[3*i+1] - fcenter.y; 
		Yz =  PointsSrc[3*i+2] - fcenter.z; 
	
		int index = CPointIndex[i];
		Xx = PointsDst[3*index+0] - ccenter.x; 
		Xy = PointsDst[3*index+1] - ccenter.y; 
		Xz = PointsDst[3*index+2] - ccenter.z; 

		// ограничение на расстояние до ближайших
		S3DCoordinate Pm;
		Pm.x = PointsSrc[3*i+0];
		Pm.y = PointsSrc[3*i+1];
		Pm.z = PointsSrc[3*i+2];
		
		S3DCoordinate Pc;
		Pc.x = PointsDst[3*index+0];
		Pc.y = PointsDst[3*index+1];
		Pc.z = PointsDst[3*index+2];

		double distXY = (Pm-Pc).norm();
		if( distXY > 40 )
			continue;

		Cov3(1,1)+=Xx*Yx;
		Cov3(2,2)+=Xy*Yy;
		Cov3(3,3)+=Xz*Yz;

		Cov3(1,2)+=Xx*Yy;
		Cov3(2,1)+=Xy*Yx;

		Cov3(1,3)+=Xx*Yz;
		Cov3(3,1)+=Xz*Yx;

		Cov3(2,3)+=Xy*Yz;
		Cov3(3,2)+=Xz*Yy;
	}
	//находим собственные вектора		
	svddecomposition(Cov3, 3, 3, 2, 2, 2, S, U, Vt);
	// матрица поворота
	S3DMatrix RCur;
	for (int i=1; i<=3; i++){
		for (int j=1; j<=3; j++){
			for (int k=1; k<=3; k++){
				RCur.a[i-1][k-1]+=U(i,j)*Vt(j,k);
			}
		}
	}
	// вектор переноса
	S3DCoordinate TCur = ccenter-RCur*fcenter;

	// дополнительный поворот и перенос
	totalTn = TCur+(RCur*totalTn);
	totalRn = RCur*totalRn;

	return TotalDist;
}
//=======================================================================//
// поворот маски в облако точек
// без копирования цветов
//=======================================================================//
template <typename T>
int MakeIcpForCloud<T>::shiftMaskToCloud( unsigned int numSrc, T* PointsSrc, T* PointsDst,  S3DCoordinate &totalT, S3DMatrix &totalR )
{
	for( unsigned int it = 0; it < numSrc; it++ )
	{
		S3DCoordinate pt;
		pt.x = PointsSrc[3*it+0];
		pt.y = PointsSrc[3*it+1];
		pt.z = PointsSrc[3*it+2];

		pt = totalR*pt + totalT;

		PointsDst[3*it+0] = pt.x;
		PointsDst[3*it+1] = pt.y;
		PointsDst[3*it+2] = pt.z;
	}
	return 0;
}
//=======================================================================//
// save cloud
//=======================================================================//
template <typename T>
int MakeIcpForCloud<T>::SaveCloud( char *fname, T* arrayCloud, unsigned int size, bool addtoFile, unsigned char R, unsigned char G, unsigned char B )
{
	FILE *file;

	if( addtoFile == false )
		file = fopen( fname, "w" );
	else
		file = fopen( fname, "at" );

	if (file) 
	{
		for( unsigned int it = 0; it < size; it++ )
		{
			S3DCoordinate pt;
			pt.x = arrayCloud[3*it+0];
			pt.y = arrayCloud[3*it+1];
			pt.z = arrayCloud[3*it+2];
			fwprintf(file, L"%f %f %f %d %d %d\n", pt.x, pt.y, pt.z, R, G, B );
		}
		fclose(file);
		//printf( "Save Size = %d %s\n", size, fname );
	}
	else 
	{
		printf( "ICP2Detector::SaveCloud: Couldn't open file %s for writing\n", fname);
		return -1;
	}
	return 0;
};
//=======================================================================//
// run icp
//=======================================================================//
template <typename T>
int MakeIcpForCloud<T>::icp( unsigned int numSrc, T* PointsSrc, unsigned int numDst, T* PointsDst, T* distansces )
{
	unsigned int NUMBIT = 10;

	ptr_Cpu<T> cloudX;
	cloudX.alloc( 3*numSrc );
	memcpy( cloudX.a, PointsSrc, 3*numSrc*sizeof( T ) );

	ptr_Cpu< unsigned int> CPointIndex;
	CPointIndex.alloc( numSrc );

	S3DCoordinate Tn = S3DCoordinate( 0, 0, 0 );
	S3DMatrix Mn = S3DMatrix::E();

	for( unsigned int it = 0; it < NUMBIT; it++ )
	{
		FindClosestPoint( numSrc, cloudX.a, numDst, PointsDst, distansces, CPointIndex.a );
		getShiftRotationFromICP( numSrc, cloudX.a, numDst, PointsDst, CPointIndex.a, Tn, Mn );
		shiftMaskToCloud( numSrc, PointsSrc, cloudX.a, Tn, Mn );

		//char fnameit[256];
		//sprintf( fnameit, "%d_out.3d", it );
		//SaveCloud( fnameit, cloudX.a, numSrc, false, 255, 255, 255 );
		//SaveCloud( fnameit, PointsDst, numDst, true, 0, 255, 0 );
	}
	FindClosestPoint( numSrc, cloudX.a, numDst, PointsDst, distansces, CPointIndex.a );

	//SaveCloud( "out.3d", cloudX.a, numSrc, false, 255, 255, 255 );
	//SaveCloud( "out.3d", PointsDst, numDst, true, 0, 255, 0 );

	return 0;
}
//=======================================================================//
//  example
//=======================================================================//
#define TYPE double
int main_test()
{
	printf( "Start\n" );

	MakeIcpForCloud<TYPE> MIFC;

	// example cloud
	unsigned int numSrc = 708;
	TYPE *PointsSrc = new TYPE[ 3*numSrc ];
	MIFC.loadCloud( "t1.3d", PointsSrc );

	unsigned int numDst = 14187;
	TYPE *PointsDst = new TYPE[ 3*numDst ];
	MIFC.loadCloud( "avgface.3d", PointsDst );

	// example rotation and offset
	S3DCoordinate T = S3DCoordinate( -30, 10, 0 );
	S3DMatrix M = S3DMatrix::E();
	double angle = 40/180.0*3.14;
	M.a[0][0] = cos( angle );
	M.a[1][1] = cos( angle );
	M.a[1][0] = -sin( angle );
	M.a[0][1] = sin( angle );

	// rotate Dst
	for( unsigned int it = 0; it < numDst; it++ )
	{
		S3DCoordinate pt;
		pt.x = PointsDst[3*it+0];
		pt.y = PointsDst[3*it+1];
		pt.z = PointsDst[3*it+2];

		pt = M*pt + T;

		PointsDst[3*it+0] = pt.x;
		PointsDst[3*it+1] = pt.y;
		PointsDst[3*it+2] = pt.z;
	}

	TYPE *distansces = new TYPE[ numSrc ];

	MIFC.SaveCloud( "in.3d", PointsSrc, numSrc, false, 255, 255, 255 );
	MIFC.SaveCloud( "in.3d", PointsDst, numDst, true, 0, 255, 0 );

	MIFC.icp( numSrc, PointsSrc, numDst, PointsDst, distansces );

	FILE *fre = fopen( "dist.txt", "w" );
	for( unsigned int it = 0; it < numSrc; it++ )
	{
		fprintf( fre, "%f\n", distansces[it] );
	}
	fclose( fre );

	printf( "end\n" );
	getch();
	return 0;
}
//=======================================================================//