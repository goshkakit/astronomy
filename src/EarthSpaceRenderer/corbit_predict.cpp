#include "corbit_predict.h"
#include "space_defines.h"

#include "common/DataConverter.h"

#include <math.h>

COrbit_Predict::COrbit_Predict(IPredictOrbitMod * _Predict, double jd, double statevec[], double sun /*= 0.0 */, double atm /*= 0.0 */)
	: Predict(_Predict)
	, jd_start(jd)
{
	parambuf_t initial(jd, statevec, sun, atm);
	params[0] = initial;
	params[1] = initial;
}

struct Scene::vec3d_t COrbit_Predict::GetNewPosition(double time)
{
	const double MIN = 1.0 / 1440.0;
	double jd = time * MIN + jd_start;

	DataConverter dc;
	params[0].param.data_end = dc.JDtoYYYYMMDD(jd);
	params[0].param.time_end = dc.SECtoHHMMSS(params[0].param.data_start, jd);

	Predict->GetNewPosition(params[0].param);

	Scene::vec3d_t result(
		  params[0].param.outX
	);

	if (need_recache(jd))
	{
		recache(jd, params[0].param.outX);
	}

	return result;
}

double COrbit_Predict::estimate_a(double statevec[])
{
	double r2 = (statevec[0] * statevec[0]) + (statevec[1] * statevec[1]) + (statevec[2] * statevec[2]); // (10^3 km)^2
	double v2 = (statevec[3] * statevec[3]) + (statevec[4] * statevec[4]) + (statevec[5] * statevec[5]); // (km/s)^2
	double ai = (2e-3/sqrt(r2)) - (v2/EARTH_GRAVITATIONAL_PARAMETER); // 1/km
	return 1. / ai; // km
}

double COrbit_Predict::estimate_T(double a)
{
	return sqrt(a*a*a / EARTH_GRAVITATIONAL_PARAMETER) * (2.0*M_PI); // s
}

void COrbit_Predict::init_param(SatParamToPredict & param, double jd, double statevec[], double sun, double atm)
{
	DataConverter dc;
	double data = dc.JDtoYYYYMMDD(jd);
	double time = dc.SECtoHHMMSS(data, jd);

	param.inX[0] = statevec[0];
	param.inX[1] = statevec[1];
	param.inX[2] = statevec[2];
	param.inX[3] = statevec[3];
	param.inX[4] = statevec[4];
	param.inX[5] = statevec[5];
	param.sun = sun;
	param.atm = atm;
	param.data_start = data;
	param.time_start = time;
	param.outX[0] = statevec[0];
	param.outX[1] = statevec[1];
	param.outX[2] = statevec[2];
	param.outX[3] = statevec[3];
	param.outX[4] = statevec[4];
	param.outX[5] = statevec[5];
	param.data_end = data;
	param.time_end = time;
}

bool COrbit_Predict::need_recache(double jd)
{
	return ((jd - params[1].jd)*86400.0) > (params[1].T*0.125);
}

void COrbit_Predict::recache(double jd, double statevec[])
{
	parambuf_t renewed(jd, statevec, params[0].param.sun, params[0].param.atm);
	params[0] = params[1];
	params[1] = renewed;
}

COrbit_Predict::parambuf_t::parambuf_t(double jd, double statevec[], double sun, double atm)
{
	static const double T_min = estimate_T(EARTH_RADIUS_KM);

	this->jd = jd;
	this->T = fmax(estimate_T(estimate_a(statevec)), T_min);
	init_param(this->param, jd, statevec, sun, atm);
}
