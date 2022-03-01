#include <stdio.h>
#include <math.h>
#include <iostream>
#include <chrono>
#include <vector>
#include "IPredictOrbitMod.h"
#include "DataConverter.h"

//=========================================================================//
// main
// set define WINDOWSCOMPILE in project
//=========================================================================//
int main()
{
	
	IPredictOrbitMod *TestP = CreatePredictOrbitMod();
	TestP->Init();
	/*
	for (int i = 20120501; i <= 20120515; i++) {
		TestP->AtmMap(7, 0.0, i, 6378136 + 120000, 6378136 + 1500000, 20000);
	}
	for (int i = 20120501; i <= 20120515; i++) {
		TestP->AtmMap(8, 0.0, i, 6378136 + 120000, 6378136 + 1500000, 20000);
	}
	for (int i = 20120501; i <= 20120515; i++) {
		TestP->AtmMap(9, 0.0, i, 6378136 + 120000, 6378136 + 1500000, 20000);
	}
	*/
	//TestP->AtmMap(6, 0.0, 20100511, 6378136 + 120000, 6378136 + 1500000, 20000);
	//TestP->AtmMap(6, 0.0, 20100512, 6378136 + 120000, 6378136 + 1500000, 20000);
	//TestP->AtmMap(6, 0.0, 20100513, 6378136 + 120000, 6378136 + 1500000, 20000);
	//TestP->AtmMap(7, 60000.0, 20100711, 6878136, 6888136, 10000);
	//TestP->AtmMap(7, 60000.0, 20100717, 6878136, 6888136, 10000);
	//TestP->AtmMap(7, 60000.0, 20100724, 6878136, 6888136, 10000);
	//TestP->AtmMap(7, 60000.0, 20100811, 6878136, 6888136, 10000);
	//TestP->AtmMap(7, 60000.0, 20100110, 6878136, 6888136, 10000);
	//TestP->AtmMap(7, 180000.0, 20100510, 6878136, 6888136, 10000);
	//TestP->AtmMap(7, 0.0, 20100511, 6878136, 6888136, 10000);
	//TestP->GravMap(5, 75, 6378100, 11378100, 100000);
	//TestP->GravMap(6, 75, 6378100, 11378100, 100000);
	//TestP->GravMap(7, 75, 6378100, 11378100, 100000);
	//TestP->GravMap(9, 75, 6378100, 11378100, 100000);
	/*
	double x[3];
	x[0] = -7.83355;
	x[1] = 0.379381;
	x[2] = -1.5914;
	TestP->countFtest(x);
	*/

	//TestP->RunTest();
	//TestP->TestAllForce();
	
	//TestP->IntegrationList();

	/*
	DataConverter DC;
	double tm = 20.0*60.0*60.0 + 22.0*60.0 + 15.21594247;
	// 73335.21594247
	double data = DC.YYYYMMDDtoJD(20120922.0) - 0.5;	// юлианская дата начала суток
	double JDmdt = data + tm / 86400.0;					// полная юлианская дата по московскому времени
	double JDutc = JDmdt - 0.125;						// полная юлианская дата по UTC

	double JDorb = JDutc + 0.125;						// полная юлианская дата по московскому времени
	double data_mdt = DC.JDtoYYYYMMDD(JDorb);
	double time_mdt = DC.SECtoHHMMSS(data_mdt, JDorb);
	printf("JDutc = %.8f\n", JDutc);
	printf("JDorb = %.8f\n", JDorb);
	printf("data_mdt = %.8f\n", data_mdt);
	printf("time_mdt = %.8f\n", time_mdt);
	*/


	/*
	//--------------------------------------// 
	SatParamToPredictJD sptrjd;
	sptrjd.inX[0] = -0.45425909352353E+01;
	sptrjd.inX[1] = -0.012375643941553E+00;
	sptrjd.inX[2] = 0.56266667231952E+01;
	sptrjd.inX[3] = 0.56546982216547E+01;
	sptrjd.inX[4] = -0.02569613222602E+00;
	sptrjd.inX[5] = 0.55120179863211E+01;
	sptrjd.atm =  0.3E-2;
	sptrjd.sun =  0.5E-05;
	sptrjd.JD_start = 2456193.22378722;
	sptrjd.JD_end = 2456194.22378472;
	clock_t tc;
	tc = clock();
	TestP->GetNewPositionForJD(sptrjd);
	tc = clock() - tc;
	std::cout << "JD Calc time is " << (double)tc / CLOCKS_PER_SEC << " sec" << std::endl;
	//--------------------------------------//
	SatParamToPredict sptr;
	sptr.inX[0] = -0.45425909352353E+01;
	sptr.inX[1] = -0.012375643941553E+00;
	sptr.inX[2] =  0.56266667231952E+01;
	sptr.inX[3] =  0.56546982216547E+01;
	sptr.inX[4] = -0.02569613222602E+00;
	sptr.inX[5] =  0.55120179863211E+01;
	sptr.atm = 0.3E-2;
	sptr.sun = 0.5E-05;
	sptr.data_start = 20120922.0;
	sptr.time_start = 202215.21594247; //(UTC+3:00)
	sptr.data_end = 20120923.0;
	sptr.time_end = 202215.0;
	tc = clock();
	TestP->GetNewPosition( sptr );
	tc = clock() - tc;
	std::cout << "Calc time is " << (double)tc / CLOCKS_PER_SEC << " sec" << std::endl;
	//--------------------------------------//

	double xfi[6];
	xfi[0] = 9.58586514817;		// 10^3 km
	xfi[1] = -0.01077540693;	// 10^3 km
	xfi[2] = -1.33992536579;	// 10^3 km
	xfi[3] = -0.80880292687;	//km/sec
	xfi[4] = 0.02038761856;		//km/sec
	xfi[5] = -5.82097466522;	//km/sec

	double dist = sqrt( (xfi[0]-sptr.outX[0])*(xfi[0]-sptr.outX[0]) + (xfi[1]-sptr.outX[1])*(xfi[1]-sptr.outX[1]) + (xfi[2]-sptr.outX[2])*(xfi[2]-sptr.outX[2]) );
	printf("GetDist = %e All good if result = 5.594739e-008\n", dist );

	dist = sqrt((xfi[0] - sptrjd.outX[0])*(xfi[0] - sptrjd.outX[0]) + (xfi[1] - sptrjd.outX[1])*(xfi[1] - sptrjd.outX[1]) + (xfi[2] - sptrjd.outX[2])*(xfi[2] - sptrjd.outX[2]));
	printf("GetDist = %e All good if result = 5.594739e-008\n", dist);
	//--------------------------------------//
	*/
	
	
	SatParamToPredict sptr;
	sptr.inX[0] = 5.40494;
	sptr.inX[1] = -3.72941;
	sptr.inX[2] = 1.23233;
	sptr.inX[3] = 2.00749;
	sptr.inX[4] = 0.45771;
	sptr.inX[5] = -7.4196;
	sptr.atm = 0.3E-2;
	sptr.sun = 0.5E-05;
	sptr.data_start = 20120501.0;
	sptr.time_start = 040000.0; //(UTC+3:00)
	sptr.data_end = 20120504.0;
	sptr.time_end = 040000.0; //202215.0
	
	auto s = std::chrono::high_resolution_clock::now();
	TestP->GetNewPosition(sptr);
	auto e = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> time = e - s;
	std::cout << "Atm load time: " << sptr.atm_load_time <<" Calc time: " << time.count() - sptr.atm_load_time << " s" << std::endl;

	double npx = 5.56581;
	double npy = -3.51253;
	double npz = -1.10962;

	if (sqrt(sptr.outX[0] * sptr.outX[0] + sptr.outX[1] * sptr.outX[1] + sptr.outX[2] * sptr.outX[2]) < 6.378136)
		std::cout << "Earth collision" << std::endl;
	else if(sqrt(sptr.outX[0] * sptr.outX[0] + sptr.outX[1] * sptr.outX[1] + sptr.outX[2] * sptr.outX[2]) > 6.378136 + 1.5)
		std::cout << "Out of atm" << std::endl;

	double dist = 1e3 * sqrt((sptr.outX[0] - npx)* (sptr.outX[0] - npx) + (sptr.outX[1] - npy) * (sptr.outX[1] - npy) + (sptr.outX[2] - npz) * (sptr.outX[2] - npz));
	std::cout << "Dist is: " << dist << " km" <<  std::endl;
	std::cout << "New position is: " << sptr.outX[0] << " " << sptr.outX[1] << " " << sptr.outX[2] << std::endl;
	
	TestP->DeInit();
	FreePredictOrbitMod( TestP );

	//int tmp;
	//scanf( "%d", &tmp );

	return 0;
}
//=========================================================================//
