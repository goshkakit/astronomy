//==============================================================================//
// Andrianov N.G.
// opbit predict module
// setup predict list satellite
//==============================================================================//

#include <time.h>
#include <stdio.h>
#include "PredictOrbitSat.h" 

// name
namespace Orbit
{
	//==============================================================================//
	// установка списка для прогнозирования
	// один спутник, во всех нитях, проверяется правильность прогноза
	//==============================================================================//
	void PredictOrbitSat::SetListSatellite_T1(  double ts, double te, SatelliteArray &ListSat  )
	{
		// начальный вектор состояния
		double x[6];
		double setx[6];

		// первый тест для проверки всего кроме атмосферы и тени от земли
		//x[0] = 0.25031033233355E+01;
		//x[1] = -0.51326544393122E+01;
		//x[2] = 0.41266667231952E+01;
		//x[3] = 0.75762130942797E+01;
		//x[4] = -0.68269673689801E+00;
		//x[5] = -0.33021161467365E+01;

		// для проверки атмосферы
		x[0] = -0.45425909352353E+01;
		x[1] = -0.012375643941553E+00;
		x[2] =  0.56266667231952E+01;
		x[3] =  0.56546982216547E+01;
		x[4] = -0.02569613222602E+00;
		x[5] =  0.55120179863211E+01;

		ListSat.AllocMemory( 6, CU_BlockXYZ );
		for( int it = 0; it < CU_BlockXYZ; it++ )
		{
			setx[0] = x[0];
			setx[1] = x[1];
			setx[2] = x[2];
			setx[3] = x[3];
			setx[4] = x[4];
			setx[5] = x[5];

			ListSat.SetSatellite( it, setx, ts );
		}
	};
	//==============================================================================//
	// установка списка для прогнозирования
	// загрузка из файла
	//==============================================================================//
	void PredictOrbitSat::SetListSatelliteFromFile(  double ts, double te, SatelliteArray &ListSat  )
	{
		printf("Load From File\n CU_BlockXYZ = %d\n", CU_BlockXYZ);

		FILE *fre;
		fre = fopen( "Data/pos.txt", "r" );

		float setx[6];
		double x[6];
		
		ListSat.AllocMemory( 6, CU_BlockXYZ );
		for( int it = 0; it < CU_BlockXYZ; it++ )
		{
			fscanf( fre, "%f %f %f %f %f %f\n", &setx[0], &setx[1], &setx[2], &setx[3], &setx[4], &setx[5] );

			if( it < 3 )
				printf( "%f %f %f %f %f %f\n", setx[0], setx[1], setx[2], setx[3], setx[4], setx[5] );

			x[0] = setx[0]/1000.0;
			x[1] = setx[1]/1000.0;
			x[2] = setx[2]/1000.0;
			x[3] = setx[3];
			x[4] = setx[4];
			x[5] = setx[5];

			ListSat.SetSatellite( it, x, ts );
		}
		fclose( fre );
		printf("Load OK\n");
	}
	//==============================================================================//
	// проверка ошибок
	// проверка результатов и сравнение их с результатами на FORTRAN
	//==============================================================================//
	void PredictOrbitSat::CheckListSatellite_T1( SatelliteArray &ListSatResCpu )
	{
		printf("***************** ORIG CHECK *******************\n");
		printf("Check List Satellite\n");
		// вектор состояния
		double x[6];
		double t[1];

		// результат вычислений
		//double xf1 =  8.90261372055;	// 10^3 km
		//double xf2 =  -0.95550183205;	// 10^3 km
		//double xf3 =  -4.09982860233;	// 10^3 km

		// для проверки атмосферы
		double xf1 = 9.58586514817;		// 10^3 km
		double xf2 = -0.01077540693;	// 10^3 km
		double xf3 = -1.33992536579;	// 10^3 km

		
		// проверка всех спутников
		printf("Error\n");
		int k = 0;
		for( int it = 0; it < CU_BlockXYZ; it++ )
		{
			ListSatResCpu.GetSatellite( it, x, t );

			double d =  GetDistTwoPoints( x[0], x[1], x[2], xf1, xf2, xf3 );
			
			printf("Find %d error = %e\n", it, d );
			k++;

			//if( k > 30 )
			//	break;
		}
	};
	//==============================================================================//
	// проверка ошибок
	//==============================================================================//
	void PredictOrbitSat::CheckListSatelliteCPU_GPU( SatelliteArray &ListSatResCpu, SatelliteArray &ListSatResGpu )
	{
		printf("**************** START CHECK CPU_GPU *******************\n");
		printf( "Check List Satellite 1 - 2 size = %d\n", CU_BlockXYZ );
		// вектор состояния
		double x1[6];
		double x2[6];
		double t[1];

		// проверка всех спутников
		for( int it = 0; it < CU_BlockXYZ; it++ )
		{
			ListSatResCpu.GetSatellite( it, x1, t );
			ListSatResGpu.GetSatellite( it, x2, t );
			
			if( it%5 == 0 && it < 300 )
			{
				double d =  GetDistTwoPoints( x1[0], x1[1], x1[2], x2[0], x2[1], x2[2] );
				printf("Find %d error = %e\n", it, d );
			}
		}
	};
	//==============================================================================//

};