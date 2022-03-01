#include "cutilNP.h"

#ifdef GPUCOMPILE

void cutilSafeCall( int err )
{
	if( err != 0 )
	 printf( "ERRR\n" );
};

#endif
