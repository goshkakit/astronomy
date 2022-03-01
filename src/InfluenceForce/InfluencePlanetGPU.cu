//==============================================================================//
// Andrianov N.G.
// opbit predict 
// module find Influence Force
// find force from planet
//==============================================================================//
#include <stdio.h>
#include "InfluenceForce.h"

#ifdef GPUCOMPILE

#include <cuda.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

// name
namespace Force
{
	//==============================================================================//
	// положение планет
	// загрузка эфемерид на GPU
	//==============================================================================//
	void InfluenceForce::cuDE403LoadFileToGPU()
	{
		printf("Load DE403 to GPU\n");
		int sizeMem = (2668442 + 10000)*sizeof( double );

		printf( "GPU MemSize = %d B\n", sizeMem );
		cudaError_t error;
		error = cudaMalloc( (void**)&d_FileMemDE403, sizeMem );
		if (error != cudaSuccess) {
			printf("Failed cudaMalloc (error code %s)!\n", cudaGetErrorString(error));
			exit(EXIT_FAILURE);
		}
		error = cudaMemcpy( d_FileMemDE403, FileMemDE403, sizeMem, cudaMemcpyHostToDevice );
		if (error != cudaSuccess) {
			printf("Failed cudaMalloc (error code %s)!\n", cudaGetErrorString(error));
			exit(EXIT_FAILURE);
		}
	}
}
//==============================================================================//

#endif

