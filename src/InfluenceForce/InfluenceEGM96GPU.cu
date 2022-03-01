//==============================================================================//
// Andrianov N.G.
// opbit predict 
// module find Influence Force
// find force from Earth harmonics
// copy array to gpu
//==============================================================================//
#include <math.h>
#include <fstream>
#include <stdio.h>
#include "InfluenceForce.h"

#ifdef GPUCOMPILE

#include "cutilNP.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

namespace Force
{
	//==============================================================================//
	// Load To GPU
	//==============================================================================//
	void InfluenceForce::cuLoadegm96ToGPU()
	{
		int memoryEgm = 6561*sizeof( double );
		cutilSafeCall( cudaMalloc( (void**)&d_egm96, memoryEgm ) );
		cutilSafeCall( cudaMemcpy( d_egm96, h_egm96, memoryEgm, cudaMemcpyHostToDevice )) ;
	};
	//==============================================================================//
};

#endif