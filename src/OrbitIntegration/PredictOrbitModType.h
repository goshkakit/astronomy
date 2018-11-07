#pragma once

struct SatParamToPredict
{
	double inX[6];
	double sun;
	double atm;

	double data_start;
	double time_start;

	double outX[6];
	double data_end;
	double time_end;
};
