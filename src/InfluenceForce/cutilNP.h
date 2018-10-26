//--------------------------------------------------------------------------------//
// функция обработки ошибок
//--------------------------------------------------------------------------------//

#include "DefineParam.h"

#ifdef GPUCOMPILE

#include <cuda_runtime.h>
#include <stdio.h>

/*
#define cutilSafeCall( err )	{ \
									if( err != 0 ) { \
										printf("CUDA ERROR\n"); \
										char erro[512]; \
										sprintf( erro, "Error at %s : %d : %s\n", __FILE__, __LINE__ , cudaGetErrorString( err ) ); \
										printf( "%s\n", erro ); \
									} \
								}
*/

void cutilSafeCall( int err );

#endif
