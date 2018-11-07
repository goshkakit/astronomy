#include "corbit_predict.h"

#include "common/DataConverter.h"

COrbit_Predict::COrbit_Predict(IPredictOrbitMod * _Predict, double jd, double statevec[], double sun /*= 0.0 */, double atm /*= 0.0 */)
	: Predict(_Predict)
	, jd_start(jd)
{
	DataConverter dc;
	param.inX[0] = statevec[0];
	param.inX[1] = statevec[1];
	param.inX[2] = statevec[2];
	param.inX[3] = statevec[3];
	param.inX[4] = statevec[4];
	param.inX[5] = statevec[5];
	param.sun = sun;
	param.atm = atm;
	param.data_start = dc.JDtoYYYYMMDD(jd);
	param.time_start = dc.SECtoHHMMSS(param.data_start, jd);
}

struct Scene::vec3d_t COrbit_Predict::GetNewPosition(double time)
{
	const double MIN = 1.0 / 1440.0;
	double jd = time * MIN + jd_start;

	DataConverter dc;
	param.data_end = dc.JDtoYYYYMMDD(jd);
	param.time_end = dc.SECtoHHMMSS(param.data_start, jd);

	Predict->GetNewPosition(param);

	return Scene::vec3d_t(
		  param.outX[0]
		, param.outX[1]
		, param.outX[2]
	);
}
