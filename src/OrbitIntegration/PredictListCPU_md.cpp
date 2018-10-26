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
	// ���������� �������� ��� ������ �� GPU
	//==============================================================================//
#ifdef GPUCOMPILE
	void PredictOrbitSat::CalcSatCloseApproach_GPU_list( double t_s, double t, Force::InfluenceForce *inIF, SatelliteArray &ListSat, SatelliteArray &ListSatVerify )
	{
		FILE *Fres = fopen( "DistInterpolate_gpu_std.txt", "w" );
		int t1 = clock();

		// ������ ���������
		std::vector< PointWithMinDist > result;

		// ��������� ����������
		IF = inIF;
		SetStepIntegration(0);
		SetTolerance();

		RunIntegrationGpu_Approach( t_s, t, initStepH, Tolerance[0], ListSat, ListSatVerify, result );

		printf("GPU OK List result = %d\n", result.size() );

		double initT[1];
		double Xverify[6];
		ListSatVerify.GetSatellite( 0, Xverify, initT );
		double atm = ListSatVerify.GetCoeffAtm( 0 );
		double sun = ListSatVerify.GetCoeffSun( 0 );
		CalcSatCloseApproach_CPU_verify( t_s, t, Xverify, atm, sun, inIF, ListSat, result );

		// ���������� ������
		fprintf( Fres, "N = %d\t NoradID = %d\n", 0, ListSatVerify.GetSatelliteID( 0 ) );

		for( int it = 0; it < result.size(); it++ )
		{
			int iarr = result[it].Nlist;
			int idn = result[it].norad;
			double d = result[it].d;
			double dv = result[it].d_verify;
			double t = result[it].t;
			double dd = d- dv;
			fprintf( Fres, "%d\t %d\t %f\t %f\t %f\t %.12f\n", iarr, idn, t, d, dv, dd );  
		}
				
		fclose( Fres );
		int t2 = clock();
		printf("TIME for all array GPU %f  ms\n", (double)(t2-t1)/CLOCKS_PER_SEC*1000.0 );
	};
#endif
	//==============================================================================//
	// �������� ������ ��������� �� ���������
	// t_s - ����� ������
	// t - ����� �����
	// inIF - �����������
	// ListSat - ������ ��������� ��� �������� - �������
	// ListSatVerify - ���������� ��������, �� ����� ��������� �� ���� ���������
	//==============================================================================//
	void PredictOrbitSat::CalcSatCloseApproach_CPU_list( double t_s, double t, Force::InfluenceForce *inIF, SatelliteArray &ListSat, SatelliteArray &ListSatVerify )
	{
		double StepSec = STEPDTIME;
		
		FILE *Fres = fopen( "DistInterpolate_cpu_std.txt", "w" );

		for( int gi = 0; gi < ListSatVerify.Block; gi++ )
		{
			// ������ ���������
			std::vector< PointWithMinDist > result;

			// ����������� ������� �� ������ ����������
			double initT[1];
			double Xverify[6];
			double atm, sun;
			ListSatVerify.GetSatellite( gi, Xverify, initT );
			atm = ListSatVerify.GetCoeffAtm( gi );
			sun = ListSatVerify.GetCoeffSun( gi );
			// set coeff for current verify satellite
			fprintf( Fres, "N = %d\t NoradID = %d atm = %f\t sun = %f\n Time start = %f\t Time end = %f\n", gi, ListSatVerify.GetSatelliteID( gi ), atm, sun, t_s, t );

			// ��������� ��������
			CalcSatCloseApproach_CPU_onetoall( t_s, t, Xverify, atm, sun, inIF, ListSat, result, StepSec );

			// ���������� ������
			for( int it = 0; it < result.size(); it++ )
			{
				int iarr = result[it].Nlist;
				int idn = result[it].norad;
				double d = result[it].d;
				double dv = result[it].d_verify;
				double t = result[it].t;
				double dd = d- dv;
				fprintf( Fres, "%d\t %d\t %f\t %f\t %f\t %.12f\n", iarr, idn, t, d, dv, dd );  
			}
		}
		fclose( Fres );
	};
	//==============================================================================//
	// ���������� �������� ��� ������ �� CPU
	// ����������� ���� ������� Xverify �� ������� ���������
	// ��������� ������ ����������� ���������� � ���������� ��������� ��������� � ��� �����
	// t_s - ����� ������, ����� � ����
	// t - ����� ��������� ������
	// inIF - �����������
	// ������ ��������� ��� ��������
	//==============================================================================//
	void PredictOrbitSat::CalcSatCloseApproach_CPU_onetoall( double t_s, double t, double *Xverify, double atmv, double sunv, Force::InfluenceForce *inIF, SatelliteArray &ListSat, std::vector< PointWithMinDist > &result, double StepSec )
	{
		int t1 = clock();

		// ��������� �� ����� � �������������
		IF = inIF;
		double initX[6];
		double initT[1];

		// points in array
		double dt = abs( t - t_s );
		dt = dt*1000.0;
		printf( "Numb sec = %f\n", dt );
		const int MaxPointN = dt/StepSec + 10; // for memory

		//alloc memory
		double *OrbitPointsArray_S0 = new double[4*MaxPointN+1];	// Nall x y z t, ....
		double *OrbitPointsArray_S1 = new double[4*MaxPointN+1];

		// ��������� ������� ������ ��������, � ������� ��� ����������
		SetStepIntegration(0);
		SetTolerance();
		double Xverify_copy[6];
		for( int i = 0; i < 6; i++ )
			Xverify_copy[i] = Xverify[i];
		
		inIF->SetSigmaAtm( atmv );
		inIF->SetSigmaSun( sunv );

		RunIntegration( t_s, Xverify_copy, t, OrbitPointsArray_S0 );

		// ��������� ������� ���� ��������� �� ������� � ���������� ����� ���������
		for( int it = 0; it < CU_BlockXYZ; it++ )
		{
			printf("Nit = %d\n", it );

			// ������ �������� �� ������
			ListSat.GetSatellite( it, initX, initT );
			double atm = ListSat.GetCoeffAtm( it );
			double sun = ListSat.GetCoeffSun( it );
			inIF->SetSigmaAtm( atm );
			inIF->SetSigmaSun( sun );

			// ��������� �������
			SetStepIntegration(0);
			SetTolerance();
			RunIntegration( initT[0], initX, t, OrbitPointsArray_S1 );

			FindSmallDistant_cpu( OrbitPointsArray_S0, OrbitPointsArray_S1, it, MaxPointN, result );
			
		}

		CalcSatCloseApproach_CPU_verify( t_s, t, Xverify, atmv, sunv, inIF, ListSat, result );

		// free memory
		delete OrbitPointsArray_S0;
		delete OrbitPointsArray_S1;

		int t2 = clock();
		printf("TIME for all array CPU %f  ms\n", (double)(t2-t1)/CLOCKS_PER_SEC*1000.0 );
	};
	//==============================================================================//
	// ����� ������ ��������, ������ �� ��������� ��������
	// ������� � ����� �������� �� CPU � ������� ��������� ����������� ���������
	// ������� ������ ����
	//==============================================================================//
	void PredictOrbitSat::CalcSatCloseApproach_CPU_verify( double t_s, double t, double *Xverify, double atmv, double sunv, Force::InfluenceForce *inIF, SatelliteArray &ListSat, std::vector< PointWithMinDist > &result )
	{
		// �������� ������������ � ������ ������ � ��������� ����������
		for( int it = 0; it < result.size(); it++ )
		{
			int iarr = result[it].Nlist;				// ������ � �������
			int idn = ListSat.GetSatelliteID( iarr );	// ����������� �����

			result[it].norad = idn;

			// ������ ��������
			double Time = result[it].t;

			//--------------------------------------//
			// �������, � ������� ���� ���������
			double Xverify_copy[6];
			for( int i = 0; i < 6; i++ )
				Xverify_copy[i] = Xverify[i];

			// ��������� �������
			SetStepIntegration(0);
			SetTolerance();
			inIF->SetSigmaAtm( atmv );
			inIF->SetSigmaSun( sunv );

			RunIntegration( t_s, Xverify_copy, Time );

			//--------------------------------------//
			// ������� �� ������
			double initX1[6];
			double initT1[1];
			ListSat.GetSatellite( iarr, initX1, initT1 );
			double atm = ListSat.GetCoeffAtm( iarr );
			double sun = ListSat.GetCoeffSun( iarr );
			// ��������� �������
			SetStepIntegration(0);
			SetTolerance();
			inIF->SetSigmaAtm( atm );
			inIF->SetSigmaSun( sun );
			RunIntegration( initT1[0], initX1, Time );


			//--------------------------------------//
			// ���������� ����� ����������
			double d =  GetDistTwoPoints( Xverify_copy[0], Xverify_copy[1], Xverify_copy[2], initX1[0], initX1[1], initX1[2] );
			d = d*1000.0; // to km
			result[it].d_verify = d;
			printf("Find min dist = %e\n", d );

			double dD = result[it].d - d;
			printf("Find min dist error = %e\n", dD );
		}
	};
	//==============================================================================//
	// ����� �������� � ���� ��������
	//==============================================================================//
	void PredictOrbitSat::FindSmallDistant_cpu( double *OrbitPointsArray_S0, double *OrbitPointsArray_S1, int iln, int MaxPointN, std::vector< PointWithMinDist > &result )
	{
		// ���������� �������
		double minDist = 10000000000;
		double minDistTime = 0;

		// ��� ������ ��������
		double LastDist;
		int Direction = 0;
		PointDist pt1; // �������
		PointDist pt2; // �������
		PointDist pt3; // ��� ���� �����
		pt1.init();
		pt2.init();
		pt3.init();

		for( int di = 0; di < MaxPointN; di++ )
		{
			int Np0 = (int)OrbitPointsArray_S0[0];				// ����� �����
			int Np1 = (int)OrbitPointsArray_S1[0];				// ����� �����
			//if( di <= Np0 && di <= Np1 )
			if( di < Np0 && di < Np1 )
			{
				double t0 = OrbitPointsArray_S0[ di*4+1 ];
				double x0 = OrbitPointsArray_S0[ di*4+1+1 ];
				double y0 = OrbitPointsArray_S0[ di*4+1+2 ];
				double z0 = OrbitPointsArray_S0[ di*4+1+3 ];

				double t1 = OrbitPointsArray_S1[ di*4+1 ];
				double x1 = OrbitPointsArray_S1[ di*4+1+1 ];
				double y1 = OrbitPointsArray_S1[ di*4+1+2 ];
				double z1 = OrbitPointsArray_S1[ di*4+1+3 ];

				if( t0 != t1 )
				{
					printf("Error T %f %f\n", t0, t1 );
				}

				// �������� ���������� ����� �������
				double d =  GetDistTwoPoints( x0, y0, z0, x1, y1, z1 );
				d = d*1000.0; // ��

				if( d < minDist )
				{
					minDist = d;
					minDistTime = t0;
				}

				//==========================================//
				// ���� ����������� ��������
				bool flagReverse = false;
				if( di > 0 )
				{
					// ����������� �� ���������� �������������
					int nowDirection;
					if( d < LastDist )
						nowDirection = 1;
					else
						nowDirection = -1;

					// ����������� ����������
					if( nowDirection < 0 && Direction > 0 )
					{
						flagReverse = true;
					}
					// ��������� �����������
					Direction = nowDirection;
				}
				// ��������� ���������� �������
				LastDist = d;
				//==========================================//

				// ����� �����
				pt3 = pt2;
				pt2 = pt1;
				// ����� ��������� ����� ����������� �����������
				pt1.t = t0;
				pt1.d = d;

				// ���� ���� ��������� �����������
				// ������������� ����������
				double InterpolateDist = 0;
				double InterpolateTime = 0;
				if( flagReverse == true && di > 1 )
				{
					double x1 = pt3.t;
					double x2 = pt2.t;
					double x3 = pt1.t;

					double y1 = pt3.d;
					double y2 = pt2.d;
					double y3 = pt1.d;

					double a = y3 - ( x3*(y2-y1) + x2*y1 - x1*y2 )/(x2 - x1);
					a = a/( x3*( x3 - x1 - x2 ) + x1*x2 );

					double b = (y2 - y1)/(x2 - x1) - a*( x1 + x2 );

					double c = (x2*y1 - x1*y2)/(x2 - x1) + a*x1*x2;

					double discr = b*b - 4.0*a*c;

					double yv = 0;
					double xv = 0;
					if( a!= 0 )
					{
						yv = -discr/4.0/a;
						xv = -b/2.0/a;
					}

					InterpolateDist = yv;
					InterpolateTime = xv;

					// add small distant
					if( InterpolateDist < PROGMINDIST )
					{
						PointWithMinDist ptd;
						ptd.Nlist = iln;
						ptd.d = InterpolateDist;
						ptd.t = InterpolateTime;
						result.push_back( ptd );
					}
				}
				//==========================================//
			}
		}
		// ok
	}
	//==============================================================================//
}