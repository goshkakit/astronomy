#ifdef WINDOWSCOMPILE

// make dll

#else

#include <stdio.h>
#include <math.h>
#include "PredictOrbitMod.h"
//=========================================================================//
// main
//=========================================================================//
int main()
{
	printf("Run Linux Mod\n");
	PredictOrbitMod *TestP = new PredictOrbitMod();
	TestP->Init();
	TestP->RunTest();
	TestP->TestAllForce();
	TestP->IntegrationList();

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

	double xfi[6];
	xfi[0] = 9.58586514817;		// 10^3 km
	xfi[1] = -0.01077540693;	// 10^3 km
	xfi[2] = -1.33992536579;	// 10^3 km
	xfi[3] = -0.80880292687;	//km/sec
	xfi[4] = 0.02038761856;		//km/sec
	xfi[5] = -5.82097466522;	//km/sec

	double dist = sqrt( (xfi[0]-sptr.outX[0])*(xfi[0]-sptr.outX[0]) + (xfi[1]-sptr.outX[1])*(xfi[1]-sptr.outX[1]) + (xfi[2]-sptr.outX[2])*(xfi[2]-sptr.outX[2]) );
	printf("GetDist = %e\nAll good if result = 5.594739e-008\n", dist );
	//--------------------------------------//

	TestP->DeInit();
	delete TestP;

	int tmp;
	scanf( "%d", &tmp );

	return 0;
}
//=========================================================================//

#endif
