//==============================================================================//
// Andrianov N.G.
// opbit predict module
//==============================================================================//
#pragma once
#include <stdio.h>
#include "PredictOrbitType.h"
#include "InfluenceForce.h"

#include "DefineParam.h"

#include <vector>

//#define _DEBUG
//#include "C:\Program Files (x86)\Visual Leak Detector\include\vld.h"


// name
namespace Orbit
{
	//================================================//
	// ��������� ��������������
	// ��� ��������� � ������� �� �����
	//================================================//
	struct CurrentIntegrateParam
	{
		double hh;
		int NY;	
		int NP;
		int IP9;

		int l;
		int l1;
		int l2;
		int ltek;

		int kkbeg;
		double delt, rr1, rr2, rotn, dd1;
		bool InvertStart;

		void init( double initStep )
		{
			hh = initStep;
			NY = 6;	
			NP = 0;
			IP9 = 2;

			l = 0;
			l1 = 0;
			l2 = 0;
			ltek = 0;

			kkbeg = 1;
			InvertStart = false;
		}
	};

	//================================================//
	// ����� � �������� ��������������
	//================================================//
	class PredictOrbitSat
	{
	private:

		//================================================//
		// ���������� � �����������
		int Nxyz;
		double initStepH;
		double *Tolerance;

		double SetStepIntegration( double h );
		void SetTolerance();
		//================================================//

		//================================================//
		// PredictOrbitConst.cpp
		const double *IC_a;
		const double *IC_a1;
		const double *IC_w;
		const double *IC_g;
		const double *IC_c;

		// ������� ��������� �������� ��� ���������� ��������������
		// �������� ������, ���������� ��������� �� ������ ��������
		// ���������� ������� ����� �������������
		double *InitParamA( );	// 120
		double *InitParamA1( ); // 40
		double *InitParamW( );	// 76
		double *InitParamG( );	// 6
		double *InitParamC( );	// 16
		// ������������� ��������
		void InitIntegrationConstant();
		// delete constant
		void DeleteIntegrationConstant();
		//================================================//

		//================================================//
		// ��������� �� ����� � �������������
		Force::InfluenceForce *IF;
		//================================================//

		//================================================//
		// PredictOrbitIntegrator.cpp
		// ���������� ��������� �������� � 8�� ������ ��� ��������������
		// ����� RK
		void GetStartOrbitPoints( OrbitPoints & OP, double h );
		// �������������� ��������� ��������
		void OrbitIntegration( OrbitPoints & OP, double t, CurrentIntegrateParam &Iparam, double *OrbitPointsArray = NULL );
		void OrbitIntegration_step( OrbitPoints & OP, double t, CurrentIntegrateParam &Iparam, double *OrbitPointsArray = NULL );

		// ������� � ������ �����
		void RightFxyzv( double t, double *xt, double *fx );
		// ������ ���������� ��������������
		void RunIntegration( double t0, double *x0, double t, double *OrbitPointsArray = NULL );
		//================================================//

#ifdef GPUCOMPILE
		//================================================//
		// GPU
		// cuda event handles
		cudaEvent_t start, stop;
		float gpuTime;

		void InitIntegrationGPU();
		void DeleteIntegrationGPU();
		// ������ ���������� �������������� �� GPU
		void RunIntegrationGpu( double t_s, double t, double h, double e, SatelliteArray &ListSat, SatelliteArray &ListSatResGpu );
		void RunIntegrationGpu_Approach( double t_s, double t, double h, double e, SatelliteArray &ListSat, SatelliteArray &ListSatVerify, std::vector< PointWithMinDist > &result );
		//================================================//
#endif

	public:
		// ������������
		PredictOrbitSat();
		~PredictOrbitSat();
		// ������� �������� � ��������
		int Init_CPU( );
		int DeInit_CPU();

		// ������� ��������
		void CalcNewArrayPosition( double t0, double *x0, Force::InfluenceForce *inIF, double *Tarr, double *resArr, int size );
		void CalcNewPosition( double t0, double *x0, double t, Force::InfluenceForce *inIF, double *OrbitPointsArray = NULL );
		double GetDistTwoPoints( double x1, double y1, double z1, double x2, double y2, double z2 );

		// ������ ���������
		void SetListSatellite_T1(  double ts, double te, SatelliteArray &ListSat );
		void SetListSatelliteFromFile(  double ts, double te, SatelliteArray &ListSat );
		void CheckListSatellite_T1( SatelliteArray &ListSatResCpu );
		void CheckListSatelliteCPU_GPU( SatelliteArray &ListSatResCpu, SatelliteArray &ListSatResGpu );
	
		// ������� ������
		void CalcNewPositionArray( double t, Force::InfluenceForce *inIF, SatelliteArray &ListSat, SatelliteArray &ListSatResCpu );
#ifdef GPUCOMPILE
		void CalcNewPositionArrayGPU( double t_s,  double t, Force::InfluenceForce *inIF, SatelliteArray &ListSat, SatelliteArray &ListSatResGpu );
#endif
		//================================================//
		// SatCloseApproach.cpp
		// �������� ������� ���������
		void CalcSatCloseApproach_CPU( double t_s, double t, Force::InfluenceForce *inIF, SatelliteArray &ListSat, SatelliteArray &ListSatResCpu );

		void CalcSatCloseApproach_CPU_list( double t_s, double t, Force::InfluenceForce *inIF, SatelliteArray &ListSat, SatelliteArray &ListSatVerify );
		void CalcSatCloseApproach_CPU_onetoall( double t_s, double t, double *Xverify, double atmv, double sunv, Force::InfluenceForce *inIF, SatelliteArray &ListSat, std::vector< PointWithMinDist > &result, double StepSec );
		void FindSmallDistant_cpu( double *OrbitPointsArray_S0, double *OrbitPointsArray_S1, int iln, int MaxPointNS,  std::vector< PointWithMinDist > &result );
		void CalcSatCloseApproach_CPU_verify( double t_s, double t, double *Xverify, double atmv, double sunv, Force::InfluenceForce *inIF, SatelliteArray &ListSat, std::vector< PointWithMinDist > &result );

#ifdef GPUCOMPILE
		void CalcSatCloseApproach_GPU( double t_s,  double t, Force::InfluenceForce *inIF, SatelliteArray &ListSat, SatelliteArray &ListSatResGpu );

		void CalcSatCloseApproach_GPU_list( double t_s, double t, Force::InfluenceForce *inIF, SatelliteArray &ListSat, SatelliteArray &ListSatVerify );
#endif

		int FindSmallDistant( double *OrbitPointsArray_S0, double *OrbitPointsArray_S1, int ns, double *ResultFindMinDist );
		double GetInterpolatePoint( double *x, double *y, int np, double xf );

		int ReCalculateStep( double *OrbitPointsArray_S0, double *OrbitPointsArray_new, double Tstart, double Tend, double Tstep );
		int FindSmallDistantStep( double *OrbitPointsArray_S0, double *OrbitPointsArray_S1, int ns, double *ResultFindMinDist );
		//================================================//
	};
}
//==============================================================================//