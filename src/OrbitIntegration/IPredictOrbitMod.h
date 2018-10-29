//==============================================================================//
// Andrianov Nikolay
// 2013.09
// Predict orbit module
// For windows
// no cudart.lib;
// post build echo copy "$(CudaToolkitBinDir)\cudart*.dll" "$(OutDir)" copy "$(CudaToolkitBinDir)\cudart*.dll" "$(OutDir)"
//==============================================================================//
#include "PredictOrbitModType.h"

#ifdef WIN32

#ifndef _IPREDICTSAT_H_
#define _IPREDICTSAT_H_

#ifndef interface
#define __STRUCT__ struct
#define interface __STRUCT__
#endif

#ifndef PURE
#define PURE = 0
#endif

#ifdef GPU_EXPORT
#define GPU_API __declspec(dllexport)
#else
#ifdef GPU_STATIC
#define GPU_API
#else
#define GPU_API //__declspec(dllimport)
#endif
#endif

//==============================================================================//
// interface
//==============================================================================//
interface IPredictOrbitMod
{
	virtual ~IPredictOrbitMod(){ printf( "~interface PredictOrbitMod()\n" ); };
	virtual int _stdcall Init() PURE;
	virtual int _stdcall DeInit() PURE;
	virtual int _stdcall RunTest() PURE;
	virtual int _stdcall TestAllForce() PURE;
	virtual int _stdcall IntegrationList() PURE;
	virtual int _stdcall GetNewPosition( SatParamToPredict &sptr ) PURE;

};

//==============================================================================//
// Creator
//==============================================================================//
class CPredictOrbitMod
{
public:
	static IPredictOrbitMod* _stdcall Create();
};
//==============================================================================//

GPU_API IPredictOrbitMod* CreatePredictOrbitMod();
GPU_API void FreePredictOrbitMod( IPredictOrbitMod* is );

#endif

#endif