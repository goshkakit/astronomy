//==============================================================================//
// Andrianov N.G.
// opbit predict module
// main
//==============================================================================//

#include <time.h>
#include <stdio.h>
#include "PredictOrbitSat.h" 

// name
namespace Orbit
{
	//==============================================================================//
	// create and delete
	//==============================================================================//
	// constructor
	PredictOrbitSat::PredictOrbitSat()
	{
		printf("PredictOrbitSat()\n");
	}
	// delete
	PredictOrbitSat::~PredictOrbitSat()
	{
		printf("~PredictOrbitSat()\n");
	}
	//==============================================================================//
	// init
	//==============================================================================//
	int PredictOrbitSat::Init_CPU( )
	{
		printf("Run Init Predict\n");
		Nxyz = 6;
		Tolerance = new double[ Nxyz ];
		
		// init Constant
		InitIntegrationConstant();

#ifdef GPUCOMPILE
		// init GPU
		InitIntegrationGPU();
#endif
		return 0;
	}
	//==============================================================================//
	// free
	//==============================================================================//
	int PredictOrbitSat::DeInit_CPU()
	{
		printf("Run Deinit Predict\n");

		delete Tolerance;
		DeleteIntegrationConstant();

#ifdef GPUCOMPILE
		// delete gpu
		DeleteIntegrationGPU();
#endif
		return 0;
	}
	//==============================================================================//
	// установка шага интегрирования
	// шаг по умолчанию 5.0E-5;
	//==============================================================================//
	double PredictOrbitSat::SetStepIntegration( double h )
	{
		if ( h > 0 )
			initStepH = h;
		else
			initStepH = 5.0E-5; //	default value

		return initStepH;
	};
	//==============================================================================//
	// установка допустимых ошибок
	// значение по умолчанию 1.0E-11;
	//==============================================================================//
	void PredictOrbitSat::SetTolerance()
	{
		// default tolerance for variations (or no tolerance at all
		// ode_tol = 1.d50 
		// default tolerance for position and velocities 1.0E-11
		for( int it = 0; it < Nxyz; it++ )
			Tolerance[it] = 1.0E-11;
	};
	//==============================================================================//
	// вычисление нового положения спутника
	//==============================================================================//
	void PredictOrbitSat::CalcNewPosition( double t0, double *x0, double t, Force::InfluenceForce *inIF, double *OrbitPointsArray )
	{
		printf("Calc New Position CPU\n");
		IF = inIF;
		SetStepIntegration(0);
		SetTolerance();
		RunIntegration( t0, x0, t, OrbitPointsArray );
	};
	//==============================================================================//
	// запуск численного интегрирования
	//==============================================================================//
	void PredictOrbitSat::RunIntegration( double t0, double *x0, double t, double *OrbitPointsArray )
	{
		//int t1 = clock();
		// структура для хранения данных
		OrbitPoints OP;
		OP.AllocMemory( Nxyz );
		// заполнение структуры с параметрами
		OP.SetInitPosition( x0, t0 );
	
		// 1
		// получение начального приближения
		GetStartOrbitPoints( OP, initStepH );
		// параметры интегрирования
		CurrentIntegrateParam Iparam;
		Iparam.init( initStepH );
		
		// 2
		// интегрирование орбиты
		//OrbitIntegration( OP, t, Iparam, OrbitPointsArray );
		OrbitIntegration_step( OP, t, Iparam, OrbitPointsArray );

		// сохранение значений
		for ( int ii = 0; ii < OP.NY; ii++ ) 
			x0[ii] = OP.X0[ii];

		// очистка структуры
		OP.DeleteMemory();
		
		//int t2 = clock();
		//printf("TIME GetStartPoint and OrbitIntegration Time %f  ms\n", (double)(t2-t1)/CLOCKS_PER_SEC*1000.0 );
	};
	//==============================================================================//
	// прогноз с сохранением в нескольких точках
	// вектор по 6 координат
	//==============================================================================//
	void  PredictOrbitSat::CalcNewArrayPosition( double t0, double *x0, Force::InfluenceForce *inIF, double *Tarr, double *resArr, int size )
	{
		printf("Calc New Array Position CPU\n");
	
		// начальные параметры
		IF = inIF;
		SetStepIntegration(0);
		SetTolerance();

		int t1 = clock();

		// структура для хранения данных
		OrbitPoints OP;
		OP.AllocMemory( Nxyz );
		// заполнение структуры с параметрами
		OP.SetInitPosition( x0, t0 );
		
		//1
		// получение начального приближения
		GetStartOrbitPoints( OP, initStepH );
		// параметры интегрирования
		CurrentIntegrateParam Iparam;
		Iparam.init( initStepH );

		// 2
		// идем по всем отметкам
		for( int it = 0; it < size; it++ )
		{
			double te = Tarr[it];
			// интегрирование орбиты
			OrbitIntegration( OP, te, Iparam );

			// сохранение значений
			for ( int ii = 0; ii < OP.NY; ii++ ) 
				resArr[ 6*it + ii ] = OP.X0[ii];

		}
		// очистка структуры
		OP.DeleteMemory();
		
		int t2 = clock();
		printf("TIME GetStartPoint and OrbitIntegration Time %f  ms\n", (double)(t2-t1)/CLOCKS_PER_SEC*1000.0 );
	};
	//==============================================================================//
	// поиск расстояния между точками
	//==============================================================================//
	double PredictOrbitSat::GetDistTwoPoints( double x1, double y1, double z1, double x2, double y2, double z2 )
	{
		double dist = sqrt( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1) );
		return dist;
	}
	//==============================================================================//
	// вычисление прогноза для списка на CPU
	//==============================================================================//
	void PredictOrbitSat::CalcNewPositionArray( double t, Force::InfluenceForce *inIF, SatelliteArray &ListSat, SatelliteArray &ListSatResCpu )
	{
		printf("***********************START CPU ***********************\n");
		printf("Calc Position for Array satellite\n");
	
		IF = inIF;

		double initX[6];
		double initT[1];
		
		int t1 = clock();
		for( int it = 0; it < CU_BlockXYZ; it++ )
		{
			// чтение спутника из списка
			ListSat.GetSatellite( it, initX, initT );

			// запускаем прогноз
			SetStepIntegration(0);
			SetTolerance();
			RunIntegration( initT[0], initX, t );
		
			// записываем результат
			ListSatResCpu.SetSatellite( it, initX, t );

			// проценты выполнения
			double per = 100.0*it/((double)CU_BlockXYZ);
			printf( "persent %f %%\r", per );
		}
		int t2 = clock();
		printf("TIME for all array CPU %f  ms\n", (double)(t2-t1)/CLOCKS_PER_SEC*1000.0 );
	};
#ifdef GPUCOMPILE
	//==============================================================================//
	// вычисление прогноза для списка на GPU
	//==============================================================================//
	void PredictOrbitSat::CalcNewPositionArrayGPU( double t_s, double t, Force::InfluenceForce *inIF, SatelliteArray &ListSat, SatelliteArray &ListSatResGpu )
	{
		printf("***********************START GPU ***********************\n");
		printf("Calc Position for Array satellite\n");
		int t1 = clock();
	
		IF = inIF;
		SetStepIntegration(0);
		SetTolerance();
	
		// запускаем прогноз
		RunIntegrationGpu( t_s, t, initStepH, Tolerance[0], ListSat, ListSatResGpu );

		int t2 = clock();
		printf("TIME for all array GPU %f  ms\n", (double)(t2-t1)/CLOCKS_PER_SEC*1000.0 );
	};
	//==============================================================================//
#endif
};