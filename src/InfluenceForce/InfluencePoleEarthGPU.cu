//==============================================================================//
// Andrianov N.G.
// opbit predict 
// module find Influence Force
// Motion Pole Earth
// set data and copy to gpu
//==============================================================================//
#include <math.h>
#include <stdio.h>
#include "InfluenceForce.h"

#ifdef GPUCOMPILE

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cutilNP.h"

namespace Force
{
	//==============================================================================//
	// копирование на GPU массивов с поправками времени
	//==============================================================================//
	void InfluenceForce::cuiers_init_gpu( double *taiutc, int sizeTai, double *table, int sizeTable, int infinals_n )
	{
		printf("Run iers_init_gpu\n");
		int memoryTAI = sizeTai*sizeof( double );
		cutilSafeCall( cudaMalloc( (void**)&d_TAIUTCDATA, memoryTAI ) );
		cutilSafeCall( cudaMemcpy( d_TAIUTCDATA, taiutc, memoryTAI, cudaMemcpyHostToDevice )) ;

		int memoryTable = sizeTable*sizeof( double );
		cutilSafeCall( cudaMalloc( (void**)&d_finals_tab, memoryTable ) );
		cutilSafeCall( cudaMemcpy( d_finals_tab, table, memoryTable, cudaMemcpyHostToDevice )) ;

		// копирование эффективного числа элементов в таблице
		cutilSafeCall( cudaMalloc( (void**)&d_finals_n, sizeof(  int ) ) );
		cutilSafeCall( cudaMemcpy( d_finals_n, &infinals_n, sizeof(  int ), cudaMemcpyHostToDevice )) ;
	}
	//==============================================================================//
};

#endif