#pragma once

#include "iorbit.h"

#include <stdio.h>
#include "OrbitIntegration/IPredictOrbitMod.h"

class COrbit_Predict : public Scene::IOrbit
{
public:
	COrbit_Predict(IPredictOrbitMod * _Predict, double jd, double statevec[], double sun = 0.0, double atm = 0.0);

	virtual struct Scene::vec3d_t GetNewPosition(double time);

	static double estimate_a(double statevec[]);
	static double estimate_T(double a);
	static void init_param(SatParamToPredict & param, double jd, double statevec[], double sun, double atm);

	bool need_recache(double jd);
	void recache(double jd, double statevec[]);

protected:
	IPredictOrbitMod * const Predict;
	const double jd_start;

	struct parambuf_t
	{
		double jd;
		double T;
		SatParamToPredict param;

		parambuf_t() {}
		parambuf_t(double jd, double statevec[], double sun, double atm);
	} params[2];
};
