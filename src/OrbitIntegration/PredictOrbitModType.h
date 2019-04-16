
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

// JD - full Julian date
struct SatParamToPredictJD
{
	double JD_start;
	double inX[6];
	double sun;
	double atm;
	double outX[6];
	double JD_end;
};