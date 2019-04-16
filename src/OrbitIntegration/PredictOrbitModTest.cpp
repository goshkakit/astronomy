#include "PredictOrbitMod.h"
#include <math.h>
#include <stdio.h>
#include <time.h>

//=========================================================================//
// run test etalon predict
//=========================================================================//
int PredictOrbitMod::RunTest()
{
	printf("___________________________________________________________________\n");
	printf("_________________________RUN TEST__________________________________\n");
	printf("___________________________________________________________________\n");

	//char *fname1 = "data/test/orbit/ic_low.txt";
	//char *fname2 = "data/test/orbit/ic_low_res.txt";
	//char *fname3 = "data/test/orbit/ic_res.txt"; // 7.178520e-008 or 5.594739e-008

	char *fname1 = "data/test/orbit/leo.txt";
	char *fname2 = "data/test/orbit/leo+24h.txt";
	char *fname3 = "data/test/orbit/leo_res.txt"; //9.940702e-009

	//char *fname1 = "data/test/orbit/geo.txt";
	//char *fname2 = "data/test/orbit/geo+240h.txt";
	//char *fname3 = "data/test/orbit/geo+240h_res.txt"; // 3.334190e-005


	//---------------------------------------------------------//
	// IN
	FILE *fre = fopen( fname1, "r" );
	const int NN = 12;
	double LOADIN[NN];
	for( int it = 0; it < NN; it++ )
	{
		double tmp = 0;
		fscanf( fre, "%lf", &tmp );
		LOADIN[it] = tmp;
	}
	fclose( fre );
	printf("Load param\n");
	for( int it = 0; it < NN; it++ )
		printf("%.14e\n", LOADIN[it] );

	// начальный вектор состояния
	double date = LOADIN[2];
	double time = LOADIN[3];
	double x[6];
	x[0] = LOADIN[4];
	x[1] = LOADIN[5];
	x[2] = LOADIN[6];
	x[3] = LOADIN[7];
	x[4] = LOADIN[8];
	x[5] = LOADIN[9];
	const double Satm = LOADIN[10];
	const double Ssun = LOADIN[11];

	double XS[6];
	for( int it = 0; it < 6; it++ )
		XS[it] = x[it];
	//---------------------------------------------------------//

	//---------------------------------------------------------//
	// RES
	fre = fopen( fname2, "r" );
	double RESULT[NN];
	for( int it = 0; it < NN; it++ )
		fscanf( fre, "%lf", &RESULT[it] );
	fclose( fre );
	printf("Res param\n");
	for( int it = 0; it < NN; it++ )
		printf("%lf\n", RESULT[it] );

	// конечный вектор
	double date_end = RESULT[2];
	double time_end = RESULT[3];
	double xfi[6];
	xfi[0] = RESULT[4];		
	xfi[1] = RESULT[5];	
	xfi[2] = RESULT[6];	
	xfi[3] = RESULT[7];	
	xfi[4] = RESULT[8];		
	xfi[5] = RESULT[9];	
	//---------------------------------------------------------//

	//---------------------------------------------------------//
	// set time
	double t, ajd0, delt0;
	double t_e;

	IForce->set_time(date, time, &ajd0, &delt0, &t );
	t_e = IForce->get_time( date_end, time_end, ajd0, delt0 );
	IForce->S_ajd0 = ajd0;
	IForce->S_delt0 = delt0;
	//---------------------------------------------------------//

	//---------------------------------------------------------//
	printf("___________________________________________________________________\n");
	// устанавливаем коэффициенты
	IForce->SetSigmaAtm( Satm );
	IForce->SetSigmaSun( Ssun );

	// вычисление прогноза по начальному приближению
	double XPOS[6];
	for( int it = 0; it < 6; it++ )
		XPOS[it] = XS[it];

	printf("Ts = %f, Te = %f\n", t, t_e );

	POSat.CalcNewPosition( t, XPOS, t_e, IForce );
	double d1 = POSat.GetDistTwoPoints( XPOS[0], XPOS[1], XPOS[2], xfi[0], xfi[1], xfi[2] );
	printf("Error Etalon XPOS = %e\n", d1 );

	fre = fopen( fname3, "w" );
	fprintf( fre, "%.14e\n", LOADIN[0] );
	fprintf( fre, "%.14e\n", LOADIN[1] );
	fprintf( fre, "%.14e\n", date_end );
	fprintf( fre, "%.14e\n", time_end );
	for( int it = 0; it < 6; it++ )
		fprintf( fre, "%.14e\n", XPOS[it] );
	fprintf( fre, "%.14e\n", Satm );
	fprintf( fre, "%.14e\n", Ssun );
	fclose( fre );
	printf("___________________________________________________________________\n");

	/*
	printf("___________________________________________________________________\n");
	printf("Test interval predoct\n");
	printf("Ts = %f, Te = %f\n", t, t_e );
	int sizePr = 4;
	double *Tarr = new double[sizePr];
	Tarr[0] = t_e-40.0;
	Tarr[1] = t_e-30.0;
	Tarr[2] = t_e-20.0;
	Tarr[3] = t_e;

	// in
	for( int it = 0; it < 6; it++ )
		XPOS[it] = XS[it];
	// out
	double *resArr = new double[6*sizePr];

	POSat.CalcNewArrayPosition( t, XPOS, IForce, Tarr, resArr, sizePr );

	// final res
	double resx = resArr[ 3*6 + 0 ];
	double resy = resArr[ 3*6 + 1 ];
	double resz = resArr[ 3*6 + 2 ];

	double d2 = POSat.GetDistTwoPoints( resx, resy, resz, xfi[0], xfi[1], xfi[2] );
	printf("Error Etalon XPOS = %e\n All good if result = 5.594739e-008\n", d2 );

	delete Tarr;
	delete resArr;
	printf("___________________________________________________________________\n\n");
	//---------------------------------------------------------//
	*/
	return 0;
};
//=========================================================================//
// run test all force
//=========================================================================//
int PredictOrbitMod::TestAllForce()
{
	printf("___________________________________________________________________\n");
	printf("___________________ TEST ALL FORCE ON CPU _________________________\n");
	printf("___________________________________________________________________\n");

	double date_end = 20120923.0;
	double time_end = 202215.0;
	double date = 20120922.0;
	double time = 202215.21594247; //(UTC+3:00)

	// reverse
	//double date = 20120923.0;
	//double time = 202215.0;
	//double date_end = 20120922.0;
	//double time_end = 202215.21594247; //(UTC+3:00)

	// число переменных
	int n = 6;
	// начальный вектор состояния
	double x[6];
	// первый тест для проверки всего 
	x[0] = -0.45425909352353E+01;
	x[1] = -0.012375643941553E+00;
	x[2] =  0.56266667231952E+01;
	x[3] =  0.56546982216547E+01;
	x[4] = -0.02569613222602E+00;
	x[5] =  0.55120179863211E+01;

	double XS[6];
	for( int it = 0; it < 6; it++ )
		XS[it] = x[it];


	double xfi[6];
	// с учетом всех влияний
	xfi[0] = 9.58586514817;		// 10^3 km
	xfi[1] = -0.01077540693;	// 10^3 km
	xfi[2] = -1.33992536579;	// 10^3 km
	xfi[3] = -0.80880292687;	//km/sec
	xfi[4] = 0.02038761856;		//km/sec
	xfi[5] = -5.82097466522;	//km/sec

	double t, ajd0, delt0;
	IForce->set_time(date, time, &ajd0, &delt0, &t );
	double t_e = IForce->get_time( date_end, time_end, ajd0, delt0 );
	IForce->S_ajd0 = ajd0;
	IForce->S_delt0 = delt0;

	IForce->TestRotation();

	printf("___________________________________________________________________\n");

	//#############################################//
	const double Satm = 0.3E-2;
	const double Ssun = 0.5E-05;

	// устанавливаем коэффициенты
	IForce->SetSigmaAtm( Satm );
	IForce->SetSigmaSun( Ssun );

	// вычисление прогноза по начальному приближению
	double XPOS[6];
	for( int it = 0; it < 6; it++ )
	{
		XPOS[it] = XS[it];
	}
	
	printf("Test interval predoct\n");
	printf("Ts = %f, Te = %f\n", t, t_e );
	int sizePr = 4;
	double *Tarr = new double[sizePr];
	Tarr[0] = t_e-40.0;
	Tarr[1] = t_e-30.0;
	Tarr[2] = t_e-20.0;
	Tarr[3] = t_e;

	double *resArr = new double[6*sizePr];
	POSat.CalcNewArrayPosition( t, XPOS, IForce, Tarr, resArr, sizePr );

	double resx = resArr[ 3*6 + 0 ];
	double resy = resArr[ 3*6 + 1 ];
	double resz = resArr[ 3*6 + 2 ];

	double d2 = POSat.GetDistTwoPoints( resx, resy, resz, xfi[0], xfi[1], xfi[2] );
	printf("Error Etalon XPOS = %e \n All good if result = 5.594739e-008\n", d2 );
	
	delete Tarr;
	delete resArr;

	printf("___________________________________________________________________\n\n");
	return 0;
};
//=========================================================================//