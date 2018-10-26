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
	// ���������� �������� ��� ������ �� CPU
	//==============================================================================//
	void PredictOrbitSat::CalcSatCloseApproach_CPU( double t_s, double t, Force::InfluenceForce *inIF, SatelliteArray &ListSat, SatelliteArray &ListSatResCpu )
	{
		printf("********************************************************\n");
		printf("************ CalcSatCloseApproach CPU ******************\n");
		int t1 = clock();

		FILE *fre2 = fopen( "LogMinDist.txt", "w" );
		fprintf( fre2, "START CPU\n" );
		fclose( fre2 );

		// ��������� �� ����� � �������������
		IF = inIF;
		double initX[6];
		double initT[1];

		// ����� ������ � ��������
		double dt = abs( t - t_s );
		dt = dt*1000.0;
		printf( "Numb sec = %f\n", dt );
		// points in array
		const int MaxPointN = dt/2+2;

		//alloc memory
		double *OrbitPointsArray_S0 = new double[4*MaxPointN+1];	// Nall x y z t, ....
		double *OrbitPointsArray_S1 = new double[4*MaxPointN+1];
		double *ResultFindMinDist = new double[3*MaxPointN+1];		// Nall, ni Ti Di, ....
		ResultFindMinDist[0] = 0;

		// for all orbit
		for( int it = 0; it < CU_BlockXYZ; it++ )
		{
			// ������ �������� �� ������
			ListSat.GetSatellite( it, initX, initT );

			// ��������� �������
			SetStepIntegration(0);
			SetTolerance();
			if( it == 0 )
				RunIntegration( initT[0], initX, t, OrbitPointsArray_S0 );
			else
				RunIntegration( initT[0], initX, t, OrbitPointsArray_S1 );

			// ���������� ���������
			ListSatResCpu.SetSatellite( it, initX, t );

			if(it > 0 )
			{
				printf("Nit = %d\n", it );
				//FindSmallDistant( OrbitPointsArray_S0, OrbitPointsArray_S1, it, ResultFindMinDist );

				
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
					if( di <= Np0 && di <= Np1 )
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
							printf("Error T %f %f\n", t0, t1 );

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
							if( InterpolateDist < 1000 )
							{
								int npos = (int)ResultFindMinDist[0];
								int indexpos = 3*npos + 1;
								ResultFindMinDist[ indexpos ] = it;
								ResultFindMinDist[ indexpos + 1] = InterpolateTime;
								ResultFindMinDist[ indexpos + 2 ] = InterpolateDist;
								ResultFindMinDist[0] += 1;
							}
						}
						//==========================================//
					}
				}
				FILE *fre2 = fopen( "LogMinDist.txt", "at" );
				fprintf( fre2, "%f\t%f\n", minDistTime, minDist );
				fclose( fre2 );
			}
		}

		// save min dist
		FILE *Fres = fopen( "DistInterpolate_cpu.txt", "w" );
		int nf = (int)ResultFindMinDist[0];
		for( int it = 0; it < nf; it++ )
		{
			int index = 3*it + 1;
			int iarr = ResultFindMinDist[index];
			int idn = ListSat.GetSatelliteID( iarr );
			fprintf( Fres, "%d\t %d\t %f\t %f\n", iarr, idn, ResultFindMinDist[index+1], ResultFindMinDist[index+2] );  
		}
		fclose( Fres );

		printf("********************************************************\n");
		printf("********************** Verify result *******************\n");
		printf("********************************************************\n");
		nf = (int)ResultFindMinDist[0];
		for( int it = 0; it < nf; it++ )
		{
			int index = 3*it + 1;
			int iarr = ResultFindMinDist[index];
			int idn = ListSat.GetSatelliteID( iarr );
			double Time = ResultFindMinDist[index+1];

			printf( "%d\t %d\t %f\t %f\n", iarr, idn, ResultFindMinDist[index+1], ResultFindMinDist[index+2] );  

			// ������ �������� �� ������
			double initX1[6];
			double initT1[1];
			ListSat.GetSatellite( iarr, initX1, initT1 );
			// ��������� �������
			SetStepIntegration(0);
			SetTolerance();
			RunIntegration( initT1[0], initX1, Time );

			// ������ �������� �� ������
			double initX0[6];
			double initT0[1];
			ListSat.GetSatellite( 0, initX0, initT0 );
			// ��������� �������
			SetStepIntegration(0);
			SetTolerance();
			RunIntegration( initT0[0], initX0, Time );

			double d =  GetDistTwoPoints( initX0[0], initX0[1], initX0[2], initX1[0], initX1[1], initX1[2] );
			d = d*1000.0; // to km
			printf("Find min dist = %e\n", d );

			double dD = ResultFindMinDist[index+2] - d;
			printf("Find min dist error = %e\n", dD );
		}
		printf("********************************************************\n");

		// free memory
		delete OrbitPointsArray_S0;
		delete OrbitPointsArray_S1;
		delete ResultFindMinDist;

		int t2 = clock();
		printf("TIME for all array CPU %f  ms\n", (double)(t2-t1)/CLOCKS_PER_SEC*1000.0 );
		printf("********************************************************\n");
	};
	//==============================================================================//
	// ���������� �������� ��� ������ �� GPU
	//==============================================================================//
#ifdef GPUCOMPILE
	void PredictOrbitSat::CalcSatCloseApproach_GPU( double t_s, double t, Force::InfluenceForce *inIF, SatelliteArray &ListSat, SatelliteArray &ListSatResGpu )
	{
		printf("***********************START GPU ***********************\n");

		FILE *fre2 = fopen( "LogMinDist.txt", "at" );
		fprintf( fre2, "START GPU\n" );
		fclose( fre2 );

		printf("CalcNewPositionArray\n");
		int t1 = clock();

		IF = inIF;
		SetStepIntegration(0);
		SetTolerance();

		//RunIntegrationGpu_Approach( t_s, t, initStepH, Tolerance[0], ListSat, ListSatResGpu );
		
		int t2 = clock();
		printf("TIME for all array GPU %f  ms\n", (double)(t2-t1)/CLOCKS_PER_SEC*1000.0 );
	};
#endif
	//==============================================================================//
	// ������� ������������ �� ��������
	//==============================================================================//
	double PredictOrbitSat::GetInterpolatePoint( double *x, double *y, int np, double xf )
	{
		double res = 0.0;
		// ��������� ������� ��������
		for( int it = 0; it < np; it++ )
		{
			double tmp = 1.0;
			// ��������� �������
			for( int j = 0; j < np; j++ )
			{
				// ���������� ���������
				if( j == it )
					j++; 

				tmp = tmp*(xf-x[j])/(x[it]-x[j]);
			}

			// �������� �� �������� � ���� �����
			tmp = tmp*y[it];

			// ���������
			res = res + tmp;
		}
		return res;
	};
	//==============================================================================//
	// �������� ������� ���������
	//==============================================================================//
	int PredictOrbitSat::FindSmallDistant( double *OrbitPointsArray_S0, double *OrbitPointsArray_S1, int ns, double *ResultFindMinDist )
	{
		// ��������� ������
		// ����� ����� ��� ������������
		const int Nlograng = 5;
		const int Nlograng_half = 2;
		// ������� ��� ������, � ������� �����
		int Ns = 3;
		// ��������� �� �����
		int Ne = 3;
		// ��� �� ������� � ��� ��������
		double dT = 1.0/1000.0;	
		double minDist = 10000000000;
		double minDistTime = 0;


		double Ts = OrbitPointsArray_S0[ Ns*4+1 ];			// ����� - ����� ������, ��������� �� ������
		int Np0 = (int)OrbitPointsArray_S0[0];				// ����� �����
		int Np1 = (int)OrbitPointsArray_S1[0];				// ����� �����
		double Te0 = OrbitPointsArray_S0[ (Np0-1-Ne)*4+1 ];	// ����� ���������, ��������� �� �����
		double Te1 = OrbitPointsArray_S1[ (Np1-1-Ne)*4+1 ];	// ����� ���������, ��������� �� �����
		// �������� ����������� �����, ����� �� ������� �� ������� ������� �����
		double Te;
		if( Te0 < Te1 )
			Te = Te0;
		else 
			Te = Te1;

		// ���������
		printf("%d\t %d\t %d\t %f\t %f\t %f\n", ns, Np0, Np1, Ts, Te, dT );

		// ������� ������� - ������� � ������� ������
		int crN0 = Ns;  // ��� ����� �� �����
		int crN1 = 1;	// ��� ����� ����� ��������
		
		// ��� ������ ��������
		double LastDist;
		int Direction = 0;
		int iterFind = 0;
		PointDist pt1; // �������
		PointDist pt2; // �������
		PointDist pt3; // ��� ���� �����
		pt1.init();
		pt2.init();
		pt3.init();

		// �������
		double *PT0 = new double[Nlograng];
		double *PT1 = new double[Nlograng];
		
		double *SX0 = new double[Nlograng];
		double *SY0 = new double[Nlograng];
		double *SZ0 = new double[Nlograng];

		double *SX1 = new double[Nlograng];
		double *SY1 = new double[Nlograng];
		double *SZ1 = new double[Nlograng];

		// ���� �� ������� � �����
		for( double iT = Ts; iT < Te; iT += dT )
		{
			// iT - ������� �����
			// �� ������ ��������� � ������ ������
			// �������� �� ����� �� �� �������
			if( crN0 >= Np0-1-Ne || crN1 >= Np1-1-Ne )
				break;

			// ��������� ����� ����� ������ �������
			while( 1 ){
				// ���������� ��������� ����� � �������
				if( iT > OrbitPointsArray_S0[ (crN0-1)*4+1 ] && iT < OrbitPointsArray_S0[ (crN0+1)*4+1 ] )	break;
				// �������� ������� ������ �������, ��������� �������
				if( iT <= OrbitPointsArray_S0[ (crN0-1)*4+1 ] )	crN0--;
				// �������� ������� ������ ������� ��������, ����������� �������
				if( iT >= OrbitPointsArray_S0[ (crN0+1)*4+1 ] )	crN0++;
			}
			while( 1 )	{
				// ���������� ��������� ����� � �������
				if( iT > OrbitPointsArray_S1[ (crN1-1)*4+1 ] && iT < OrbitPointsArray_S1[ (crN1+1)*4+1 ] )	break;
				// �������� ������� ������ �������, ��������� �������
				if( iT <= OrbitPointsArray_S1[ (crN1-1)*4+1 ] )	crN1--;
				// �������� ������� ������ ������� ��������, ����������� �������
				if( iT >= OrbitPointsArray_S1[ (crN1+1)*4+1 ] )	crN1++;
			}

			// �������� �� ����� �� �� �������
			if( crN0 >= Np0-1-Ne || crN1 >= Np1-1-Ne )
				break;

			if( crN0 <= Nlograng_half || crN1 <= Nlograng_half )
				continue;

			// ������ ��� ������������
			// crN0, crN1 - ������� ����������� �����
			// iT - �����
			// �������� ����� ��� ������������, ����� � ��� �� �� ��������
			int k = 0;
			for( int it = crN0 - Nlograng_half; it <= crN0 + Nlograng_half; it++)
			{
				PT0[k] = OrbitPointsArray_S0[ it*4+1 ];
				SX0[k] = OrbitPointsArray_S0[ it*4+1+1 ];
				SY0[k] = OrbitPointsArray_S0[ it*4+1+2 ];
				SZ0[k] = OrbitPointsArray_S0[ it*4+1+3 ];
				k++;
			}

			k = 0;
			for( int it = crN1 - Nlograng_half; it <= crN1 + Nlograng_half; it++)
			{
				PT1[k] = OrbitPointsArray_S1[ it*4+1 ];
				SX1[k] = OrbitPointsArray_S1[ it*4+1+1 ];
				SY1[k] = OrbitPointsArray_S1[ it*4+1+2 ];
				SZ1[k] = OrbitPointsArray_S1[ it*4+1+3 ];
				k++;
			}
			// ����� �������� ����� � �������������
			double x0 = GetInterpolatePoint( PT0, SX0, Nlograng, iT );
			double y0 = GetInterpolatePoint( PT0, SY0, Nlograng, iT );
			double z0 = GetInterpolatePoint( PT0, SZ0, Nlograng, iT );

			double x1 = GetInterpolatePoint( PT1, SX1, Nlograng, iT );
			double y1 = GetInterpolatePoint( PT1, SY1, Nlograng, iT );
			double z1 = GetInterpolatePoint( PT1, SZ1, Nlograng, iT );

			// �������� ���������� ����� �������
			double d =  GetDistTwoPoints( x0, y0, z0, x1, y1, z1 );
			d = d*1000.0; // ��

			if( d < minDist )
			{
				minDist = d;
				minDistTime = iT;

			}

			//==========================================//
			// ���� ����������� ��������
			bool flagReverse = false;
			if( iterFind > 0 )
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

			// ����� �����
			pt3 = pt2;
			pt2 = pt1;
			// ����� ��������� ����� ����������� �����������
			pt1.t = iT;
			pt1.d = d;

			// ���� ���� ��������� �����������
			// ������������� ����������
			double InterpolateDist = 100000000000;
			double InterpolateTime = 0;
			if( flagReverse == true && iterFind > 2 )
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
			}
			//==========================================//
			// ��������� ����
			if( d >= 100 && d < 300 )
				dT = 1.0/1000.0;
			else if( d < 100 )
				dT = 0.1/1000.0;
			else
				dT = 5.0/1000.0;

			if( flagReverse == true && InterpolateDist < 1000 )
			{
				int npos = (int)ResultFindMinDist[0];
				int indexpos = 3*npos + 1;
				ResultFindMinDist[ indexpos ] = ns;
				ResultFindMinDist[ indexpos + 1] = InterpolateTime;
				ResultFindMinDist[ indexpos + 2 ] = InterpolateDist;
				ResultFindMinDist[0] += 1;
			}
			iterFind ++;
		}

		// free all
		delete PT0;
		delete PT1;

		delete SX0;
		delete SY0;
		delete SZ0;

		delete SX1;
		delete SY1;
		delete SZ1;

		// write min dist
		FILE *fre2 = fopen( "LogMinDist.txt", "at" );
		fprintf( fre2, "%f\t%f\n", minDistTime, minDist );
		fclose( fre2 );

		return Np1;
	};
	//==============================================================================//


	//##############################################################################//
	//				���������� ���� - ���������� ��� 
	//##############################################################################//
	//==============================================================================//
	// �������� ������� ���������
	//
	// �������� � ��� � ��� ������
	// double Tstart
	// double Tend
	// double Tstep 
	//==============================================================================//
	int PredictOrbitSat::ReCalculateStep( double *OrbitPointsArray_S0, double *OrbitPointsArray_new, double Tstart, double Tend, double Tstep )
	{
		// ����� ����� ��� ������������
		const int Nlograng = 5;
		const int Nlograng_half = 2;

		int Np0 = (int)OrbitPointsArray_S0[0];			// ����� �����
		double Tas = OrbitPointsArray_S0[ 1 ];			// ����� ������ ������
		double Tae = OrbitPointsArray_S0[ (Np0-1)*4+1 ];// ����� ��������� ������

		// �������� �� ����� �������
		if( Np0 <= Nlograng ){
			printf( "small size array = %d\n", Np0 );
			return -1;
		}
		// ��������� ����� ������, ������ � ������
		if( Tstart < Tas || Tstart > Tae )
		{
			printf( "Tstart not correct\n" );
			return -1;
		}
		// ����� �����
		int lsize = (Tend-Tstart)/Tstep;
		printf("lsize = %d\n", lsize );

		// �������
		double *PT0 = new double[Nlograng];
		double *SX0 = new double[Nlograng];
		double *SY0 = new double[Nlograng];
		double *SZ0 = new double[Nlograng];
		
		// ���� �� ������� � �����
		int index = 0;
		int crN0 = 0; 
		for( double iT = Tstart; iT < Tend; iT += Tstep )
		{
			// ��������� ����� ����� ������ �������
			bool fail = false;
			while( 1 ){
				// ���������� ��������� ����� � �������
				if( iT >= OrbitPointsArray_S0[ crN0*4+1 ] && iT <= OrbitPointsArray_S0[ (crN0+1)*4+1 ] ) break;
				// �������� ������� ������ �������, ��������� �������
				if( iT < OrbitPointsArray_S0[ crN0*4+1 ] )	crN0--;
				// �������� ������� ������ ������� ��������, ����������� �������
				if( iT > OrbitPointsArray_S0[ (crN0+1)*4+1 ] )	crN0++;

				// ����� �� �������
				if( crN0 < 0 || crN0+1 >= Np0 ) 
				{
					fail = true;
					break;
				}
			}

			// ������ ������ 
			if( fail )
			{
				printf("find point fail\n");
				break;
			}

			// ������ ��� ������������
			// crN0 - ������� ����������� �����
			// iT - �����
			// �������� ����� ��� ������������, ����� � ��� �� �� ��������, ��� ����� ��������
			int posStart = crN0 - Nlograng_half;
			int posEnd = crN0 + Nlograng_half;
		
			if( posStart < 0 ){
				posStart = 0;
				posEnd = posStart + (Nlograng-1);
			}
		
			if( posEnd > Np0-1 ){
				posEnd = Np0-1;
				posStart = posEnd - (Nlograng-1);
			}
			
			// �������� ������ �����
			int k = 0;
			for( int it = posStart; it <= posEnd; it++)
			{
				PT0[k] = OrbitPointsArray_S0[ it*4+1 ];
				SX0[k] = OrbitPointsArray_S0[ it*4+1+1 ];
				SY0[k] = OrbitPointsArray_S0[ it*4+1+2 ];
				SZ0[k] = OrbitPointsArray_S0[ it*4+1+3 ];
				k++;
			}

			// ����� �������� ����� � �������������
			double x0 = GetInterpolatePoint( PT0, SX0, Nlograng, iT );
			double y0 = GetInterpolatePoint( PT0, SY0, Nlograng, iT );
			double z0 = GetInterpolatePoint( PT0, SZ0, Nlograng, iT );

			// ��������� �����
			OrbitPointsArray_new[ index*4+1   ] = iT;
			OrbitPointsArray_new[ index*4+1+1 ] = x0;
			OrbitPointsArray_new[ index*4+1+2 ] = y0;
			OrbitPointsArray_new[ index*4+1+3 ] = z0;
			index++;
		}
		// ����� �����
		OrbitPointsArray_new[0] = index;

		printf("size index = %d\n", index );

		// free all
		delete PT0;
		delete SX0;
		delete SY0;
		delete SZ0;
		
		return index;
	};
	//==============================================================================//
	// ����� �� ������ � ���������� �����
	//==============================================================================//
	int PredictOrbitSat::FindSmallDistantStep( double *OrbitPointsArray_S0, double *OrbitPointsArray_S1, int ns, double *ResultFindMinDist )
	{
		double Ts = OrbitPointsArray_S0[ 1 ];		// ����� - ����� ������, ��������� �� ������
		int Np0 = (int)OrbitPointsArray_S0[0];		// ����� �����
		int Np1 = (int)OrbitPointsArray_S1[0];		// ����� �����
		
		// �������� ����������� ����� ������
		int numbPt = Np0;
		if( numbPt > Np1 )
			numbPt = Np1;

		// min dist
		double minDist = 10000000000.0;
		double minDistTime = 0;
		int minit = 0;

		// ��� ������ ��������
		double LastDist = 10000000000.0;
		int Direction = 0;
		int iterFind = 0;
		PointDist pt1; // �������
		PointDist pt2; // �������
		PointDist pt3; // ��� ���� �����
		pt1.init();
		pt2.init();
		pt3.init();

		// ���� �� ������� � �����
		for( int it = 0; it < numbPt; it++ )
		{
			// ����� �������� ����� � �������������
			double iT = OrbitPointsArray_S0[ it*4+1+0 ];
			double x0 = OrbitPointsArray_S0[ it*4+1+1 ];
			double y0 = OrbitPointsArray_S0[ it*4+1+2 ];
			double z0 = OrbitPointsArray_S0[ it*4+1+3 ];

			double x1 = OrbitPointsArray_S1[ it*4+1+1 ];
			double y1 = OrbitPointsArray_S1[ it*4+1+2 ];
			double z1 = OrbitPointsArray_S1[ it*4+1+3 ];

			// �������� ���������� ����� �������
			double d =  GetDistTwoPoints( x0, y0, z0, x1, y1, z1 );
			d = d*1000.0; // ��

			if( d < minDist )
			{
				minDist = d;
				minDistTime = iT;
				minit = it;
			}

			//==========================================//
			// ���� ����������� ��������
			bool flagReverse = false;
			if( iterFind > 0 )
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

			// ����� �����
			pt3 = pt2;
			pt2 = pt1;
			// ����� ��������� ����� ����������� �����������
			pt1.t = iT;
			pt1.d = d;

			// ���� ���� ��������� �����������
			// ������������� ����������
			double InterpolateDist = 100000000000;
			double InterpolateTime = 0;
			if( flagReverse == true && iterFind > 2 )
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
			}

			if( flagReverse == true && InterpolateDist < 100 )
			{
				int npos = (int)ResultFindMinDist[0];
				int indexpos = 3*npos + 1;
				ResultFindMinDist[ indexpos ] = ns;
				ResultFindMinDist[ indexpos + 1] = InterpolateTime;
				ResultFindMinDist[ indexpos + 2 ] = InterpolateDist;
				ResultFindMinDist[0] += 1;
			}
			iterFind ++;
		}

		// write min dist
		FILE *fre2 = fopen( "LogMinDist.txt", "at" );
		fprintf( fre2, "%f\t%f\t%d\n", minDistTime, minDist, minit );
		fclose( fre2 );

		return 0;
	}
	//==============================================================================//
}

/*
const int MaxMinDist = 100000;

		double *Orb0 = new double[4*100000];
		double *Orb1 = new double[4*100000];

		double *OrbitPointsArray_S0 = new double[4*MaxPointN+1]; // 320KB �� 10��� �����
		double *OrbitPointsArray_S1 = new double[4*MaxPointN+1];
		// Nall, ni Ti Di, ....
		double *ResultFindMinDist = new double[3*MaxMinDist+1]; 
		ResultFindMinDist[0] = 0;

		for( int it = 0; it < CU_BlockXYZ; it++ )
		{
			// ������ �������� �� ������
			ListSat.GetSatellite( it, initX, initT );

			// ��������� �������
			SetStepIntegration(0);
			SetTolerance();
			if( it == 0 )
				RunIntegration( initT[0], initX, t, OrbitPointsArray_S0 );
			else
				RunIntegration( initT[0], initX, t, OrbitPointsArray_S1 );

			if(it > 0 )
			{
				FindSmallDistant( OrbitPointsArray_S0, OrbitPointsArray_S1, it, ResultFindMinDist );

				// new
				//double t1;
				//double t2;
				//double dt = 86.400/4.0;
				//
				//int NN = 0;
				//while(1)
				//{
				//	t1 = initT[0] + NN*dt;
				//	t2 = initT[0] + (NN+1)*dt;

				//	if( t1 >= t )
				//		break;

				//	ReCalculateStep( OrbitPointsArray_S0, Orb0, t1, t2, 0.001 );
				//	ReCalculateStep( OrbitPointsArray_S1, Orb1, t1, t2, 0.001 );
				//	FindSmallDistantStep( Orb0, Orb1, it, ResultFindMinDist );
				//	NN++;
				//}
			}
			// ���������� ���������
			ListSatResCpu.SetSatellite( it, initX, t );

			// �������� ����������
			//double per = 100.0*it/((double)CU_BlockXYZ);
			//printf( "persent %f %\r", per );
		}

		FILE *Fres = fopen( "DistInterpolate_cpu.txt", "w" );
		int nf = (int)ResultFindMinDist[0];
		for( int it = 0; it < nf; it++ )
		{
			int index = 3*it + 1;
			int iarr = ResultFindMinDist[index];
			int idn = ListSat.GetSatelliteID( iarr );
			fprintf( Fres, "%d\t %d\t %f\t %f\n", iarr, idn, ResultFindMinDist[index+1], ResultFindMinDist[index+2] );  
		}
		fclose( Fres );

		delete OrbitPointsArray_S0;
		delete OrbitPointsArray_S1;
		delete ResultFindMinDist;

		delete Orb0;
		delete Orb1;

		int t2 = clock();
		printf("TIME for all array CPU %f  ms\n", (double)(t2-t1)/CLOCKS_PER_SEC*1000.0 );
*/