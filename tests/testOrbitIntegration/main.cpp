#include <stdio.h>
#include <math.h>
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
	//TestP->RunTest();
	//TestP->TestAllForce();

	//TestP->IntegrationList();

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

	//--------------------------------------// 
	SatParamToPredictJD sptrjd;
	sptrjd.inX[0] = -0.45425909352353E+01;
	sptrjd.inX[1] = -0.012375643941553E+00;
	sptrjd.inX[2] = 0.56266667231952E+01;
	sptrjd.inX[3] = 0.56546982216547E+01;
	sptrjd.inX[4] = -0.02569613222602E+00;
	sptrjd.inX[5] = 0.55120179863211E+01;
	sptrjd.atm = 0.3E-2;
	sptrjd.sun = 0.5E-05;
	sptrjd.JD_start = 2456193.22378722;
	sptrjd.JD_end = 2456194.22378472;
	TestP->GetNewPositionForJD(sptrjd);

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
	TestP->GetNewPosition( sptr );
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
	
	TestP->DeInit();
	FreePredictOrbitMod( TestP );

	//int tmp;
	//scanf( "%d", &tmp );

	return 0;
}
//=========================================================================//
