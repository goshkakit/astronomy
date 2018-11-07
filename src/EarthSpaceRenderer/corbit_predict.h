#pragma once

#include "iorbit.h"

#include <stdio.h>
#include "OrbitIntegration/IPredictOrbitMod.h"

class COrbit_Predict : public Scene::IOrbit
{
public:
	COrbit_Predict(IPredictOrbitMod * _Predict, double jd, double statevec[], double sun = 0.0, double atm = 0.0);

	virtual struct Scene::vec3d_t GetNewPosition(double time);

protected:
	IPredictOrbitMod * const Predict;
	const double jd_start;
	SatParamToPredict param;
};
