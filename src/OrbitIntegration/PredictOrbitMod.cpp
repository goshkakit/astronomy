#include "PredictOrbitMod.h"
#include <math.h>
#include <stdio.h>
#include <time.h>

#ifdef WIN32
//==============================================================================//
// Creator
//==============================================================================//
IPredictOrbitMod *CPredictOrbitMod::Create()
{
	return (IPredictOrbitMod *) new PredictOrbitMod();
};
//==============================================================================//
// Alternative creator
//==============================================================================//
IPredictOrbitMod* CreatePredictOrbitMod()
{
	return CPredictOrbitMod::Create();
}
void FreePredictOrbitMod( IPredictOrbitMod* is )
{
	delete is;
}
#endif

//=========================================================================//
// init all
//=========================================================================//
int PredictOrbitMod::Init()
{
	// init
	IForce = new Force::InfluenceForce;
	IForce->Init_CPU();
	POSat.Init_CPU( );

	return 0;
}
//=========================================================================//
// delete all
//=========================================================================//
int PredictOrbitMod::DeInit()
{
	// delete
	IForce->DeInit_CPU();
	delete IForce;
	POSat.DeInit_CPU();
	return 0;
}
//=========================================================================//
// одиночный прогноз
//=========================================================================//
int PredictOrbitMod::GetNewPosition( SatParamToPredict &sptr )
{
	printf("___________________________________________________________________\n");
	printf("__________________________ RUN PREDICT ____________________________\n");
	printf("___________________________________________________________________\n");

	double date_end = sptr.data_end;
	double time_end = sptr.time_end;
	double date = sptr.data_start;
	double time = sptr.time_start;

	// начальный вектор состояния
	double XS[6];
	for( int it = 0; it < 6; it++ )
		XS[it] = sptr.inX[it];

	double t, ajd0, delt0;
	IForce->set_time(date, time, &ajd0, &delt0, &t );
	double t_e = IForce->get_time( date_end, time_end, ajd0, delt0 );
	IForce->S_ajd0 = ajd0;
	IForce->S_delt0 = delt0;

	const double Satm = sptr.atm;
	const double Ssun = sptr.sun;
	IForce->SetSigmaAtm( Satm );
	IForce->SetSigmaSun( Ssun );

	// вычисление прогноза по начальному приближению
	POSat.CalcNewPosition( t, XS, t_e, IForce );

	for( int it = 0; it < 6; it++ )
	{
		sptr.outX[it] = XS[it];
	}
	printf("___________________________________________________________________\n\n");
	return 0;
}
//=========================================================================//
// заполнение списка спутников
// интегрирование списка на CPU и на GPU
//=========================================================================//
int PredictOrbitMod::IntegrationListNorad()
{
	printf("___________________________________________________________________\n");
	printf("____________________ TEST LIST ON CPU AND GPU  ____________________\n");
	printf("___________________________________________________________________\n");
	double date_end = 20120923.0;
	double time_end = 202215.0;

	//double date_end = 20120922.0;
	//double time_end = 222215.0;

	double date = 20120922.0;
	double time = 202215.21594247; //(UTC+3:00)

	const double Satm = 0.3E-2;
	const double Ssun = 0.5E-05;

	// модуль воздействий
	double t_e;
	double t, ajd0, delt0;
	IForce->set_time(date, time, &ajd0, &delt0, &t );
	t_e = IForce->get_time( date_end, time_end, ajd0, delt0 );
	IForce->S_ajd0 = ajd0;
	IForce->S_delt0 = delt0;
	printf("IForce->set_time: %f %f %f CH = %f\n", t, ajd0, delt0, t/3.6 );

	// устанавливаем коэффициенты
	IForce->SetSigmaAtm( Satm );
	IForce->SetSigmaSun( Ssun );

	// результаты		
	Orbit::SatelliteArray ListSat;
	Orbit::SatelliteArray ListSatResCpu;
	Orbit::SatelliteArray ListSatResGpu;
	ListSatResCpu.AllocMemory( 6, CU_BlockXYZ );
	ListSatResGpu.AllocMemory( 6, CU_BlockXYZ );

	// 1 list predict
	// задание положения для проверки всех воздействий
	POSat.SetListSatellite_T1( t, t_e, ListSat );
	POSat.CalcNewPositionArray( t_e, IForce, ListSat, ListSatResCpu );
	//POSat.CalcNewPositionArrayGPU( t, t_e, IForce, ListSat, ListSatResGpu ); // 1.071956e-007
	//printf("CPU\n");
	//POSat.CheckListSatellite_T1( ListSatResCpu );
	//printf("GPU\n");
	//POSat.CheckListSatellite_T1( ListSatResGpu );
	POSat.CheckListSatelliteCPU_GPU( ListSatResCpu, ListSatResGpu ); //1.366815e-008

	// 2 find minimum dist
	//POSat.SetListSatelliteFromFile( t, t_e, ListSat );
	//POSat.CalcSatCloseApproach_CPU( t_e, IForce, ListSat, ListSatResCpu );
	//POSat.CalcSatCloseApproach_GPU( t, t_e, IForce, ListSat, ListSatResGpu );
	//POSat.CheckListSatelliteCPU_GPU( ListSatResCpu, ListSatResGpu );

	ListSat.DeleteMemory();
	ListSatResGpu.DeleteMemory();
	ListSatResCpu.DeleteMemory();
	
	printf("___________________________________________________________________\n");
	return 0;
}
//=========================================================================//
// file
// orb131007.sv
// example N D T x y x vx vy vz atm, sun;
// 90184 20131008.0 030000.000000  3.91004097437e+01 -1.36725037769e+01  1.92007247555e+00  9.92869677713e-01  2.90981049268e+00  4.31101216322e-01  3.00000000000e-02  1.00000000000e-05
//=========================================================================//
int PredictOrbitMod::IntegrationList()
{
	//IntegrationListNorad();
	//return 0;

	printf( "********************************************************\n");
	printf( "************************RUN LIST************************\n");
	printf( "********************************************************\n");

	// массивы орбит
	Orbit::SatelliteArray ListSat;
	Orbit::SatelliteArray ListSatResCpu;
	Orbit::SatelliteArray ListSatResGpu;
	ListSat.AllocMemory( 6, CU_BlockXYZ );
	ListSatResCpu.AllocMemory( 6, CU_BlockXYZ );
	ListSatResGpu.AllocMemory( 6, CU_BlockXYZ );

	// for test
	Orbit::SatelliteArray ListSatVerify;
	int numbsatverify = 1;
	ListSatVerify.AllocMemory( 6, numbsatverify );
	
	// время старта
	double date = 0;
	double time = 0;
	double date_end = 0;
	double time_end = 0;
	double t_e;
	double t, ajd0, delt0;
	double Satm;
	double Ssun;

	// файл с орбитами
	//FILE *fre = fopen( "Data/orb131007.sv", "r" );
	FILE *fre = fopen( "Data/orb140224.sv", "r" );
	for( int it = 0; it < CU_BlockXYZ; it++ )
	{
		double xvec[6];
		int idSat;
		double D, T, x, y, z, vx, vy, vz, atm, sun;
		fscanf( fre, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &idSat, &D, &T, &x, &y, &z, &vx, &vy, &vz, &atm, &sun );

		if( it < 3 )
			printf( "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", idSat, D, T, x, y, z, vx, vy, vz, atm, sun );

		// первая итерация - получаем время
		if( it == 0 )
		{
			date = D;
			time = T;
			date_end = date + 7.0; // на сутки вперед
			time_end = time;
			IForce->set_time(date, time, &ajd0, &delt0, &t );
			t_e = IForce->get_time( date_end, time_end, ajd0, delt0 );
			IForce->S_ajd0 = ajd0;
			IForce->S_delt0 = delt0;

			Satm = atm;
			Ssun = sun;
			IForce->SetSigmaAtm( Satm );
			IForce->SetSigmaSun( Ssun );
		}

		// запись спутника в список
		xvec[0] = x;
		xvec[1] = y;
		xvec[2] = z;
		xvec[3] = vx;
		xvec[4] = vy;
		xvec[5] = vz;
		ListSat.SetSatellite( it, xvec, t );
		ListSat.SetSatelliteID( it, idSat );
		ListSat.SetCoeffAtmSun( it, atm, sun );

		// list for verify
		if( it < numbsatverify )
		{
			ListSatVerify.SetSatellite( it, xvec, t );
			ListSatVerify.SetSatelliteID( it, idSat );
			ListSatVerify.SetCoeffAtmSun( it, atm, sun );
		}
	}
	fclose(fre);

	printf( "********************************************************\n");
	printf( "Set Param:\n Satm = %f\n Ssun = %f\n", Satm, Ssun );
	printf( " ajd0 = %f\n delt0 = %f\n", ajd0, delt0 );
	printf( " t = %f\n te = %f\n", t, t_e );
	printf( " Date = %f\n Time = %f\n", date, time );
	printf( " Date_end = %f\n Time_end = %f\n", date_end, time_end );
	printf( " Size = %d\n", CU_BlockXYZ );
	printf( "********************************************************\n");

	// 2 find minimum dist
	//POSat.CalcSatCloseApproach_CPU( t, t_e, IForce, ListSat, ListSatResCpu );
	//POSat.CalcSatCloseApproach_CPU_list( t, t_e, IForce, ListSat, ListSatVerify );

	//POSat.CalcSatCloseApproach_GPU( t, t_e, IForce, ListSat, ListSatResGpu );
	//POSat.CalcSatCloseApproach_GPU_list( t, t_e, IForce, ListSat, ListSatVerify );
	//POSat.CheckListSatelliteCPU_GPU( ListSatResCpu, ListSatResGpu );

	// 1 list predict
	// задание положения для проверки всех воздействий
	//POSat.CalcNewPositionArray( t_e, IForce, ListSat, ListSatResCpu );
	//POSat.CalcNewPositionArrayGPU( t, t_e, IForce, ListSat, ListSatResGpu );
	//POSat.CheckListSatelliteCPU_GPU( ListSatResCpu, ListSatResGpu );

	printf( "********************************************************\n");
	// free result
	ListSat.DeleteMemory();
	ListSatResGpu.DeleteMemory();
	ListSatResCpu.DeleteMemory();

	return 0;
};
//=========================================================================//
