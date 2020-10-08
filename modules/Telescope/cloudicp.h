//=======================================================================//
// Andrianov N.G.
// icp for two cloud
//=======================================================================//
#include <vector>

#include "common/mytypes.h"
using namespace VecMath;

template <typename T>
class ptr_Cpu
{
public:

	T* a;
	int size;

	ptr_Cpu()
	{
	};
	~ptr_Cpu()
	{
		delete[] a;
		printf( "delete arr\n" );
	};

	int alloc( int ln )
	{
		a = new T[ln];
		size = ln;
		printf( "init arr\n" );

		return 0;
	};
};

template <typename T>
class MakeIcpForCloud
{
private:

	T FindClosestPoint( unsigned int numSrc, T* PointsSrc, unsigned int numDst, T* PointsDst, T* distansces, unsigned int *CPointIndex );
	T getShiftRotationFromICP( unsigned int numSrc, T* PointsSrc, unsigned int numDst, T* PointsDst, unsigned int *CPointIndex, S3DCoordinate &totalTn, S3DMatrix &totalRn );
	int shiftMaskToCloud( unsigned int numSrc, T* PointsSrc, T* PointsDst, S3DCoordinate &totalT, S3DMatrix &totalR );

public:

	MakeIcpForCloud();
	~MakeIcpForCloud();

	int loadCloud( char *fname, T* arrayPoint );
	int SaveCloud( char *fname, T* arrayCloud, unsigned int size, bool addtoFile, unsigned char R, unsigned char G, unsigned char B  );
	int icp( unsigned int numSrc, T* PointsSrc, unsigned int numDst, T* PointsDst, T* distansces );
};
