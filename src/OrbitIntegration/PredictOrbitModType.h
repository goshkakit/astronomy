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

	double atm_load_time;
};

// JD - full Julian date
struct SatParamToPredictJD
{
	double JD_start;
	double inX[6];
	double sun;	//0.5E-05;
	double atm; //0.3E-2;
	double outX[6];
	double JD_end;
	double atm_load_time;
};

// Example use
//IPredictOrbitMod *TestP = CreatePredictOrbitMod();
//TestP->Init();
//
//SatParamToPredictJD sptrjd;
//sptrjd.inX[0] = -0.45425909352353E+01;
//sptrjd.inX[1] = -0.012375643941553E+00;
//sptrjd.inX[2] = 0.56266667231952E+01;
//sptrjd.inX[3] = 0.56546982216547E+01;
//sptrjd.inX[4] = -0.02569613222602E+00;
//sptrjd.inX[5] = 0.55120179863211E+01;
//sptrjd.atm = 0.3E-2;
//sptrjd.sun = 0.5E-05;
//sptrjd.JD_start = 2456193.22378722;
//sptrjd.JD_end = 2456194.22378472;
//TestP->GetNewPositionForJD(sptrjd);

//double xfi[6];
//xfi[0] = 9.58586514817;		// 10^3 km
//xfi[1] = -0.01077540693;	// 10^3 km
//xfi[2] = -1.33992536579;	// 10^3 km
//xfi[3] = -0.80880292687;	//km/sec
//xfi[4] = 0.02038761856;		//km/sec
//xfi[5] = -5.82097466522;	//km/sec
