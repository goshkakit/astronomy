//==============================================================================//
// Andrianov N.G.
// opbit predict 
// module find Influence Force

// files ftp://maia.usno.navy.mil/ser7
//==============================================================================//
#pragma once
#ifndef _InfluenceForce_H_
#define _InfluenceForce_H_

#include "DefineParam.h"
#include <vector>

#ifdef GPUCOMPILE
#include "cuda_runtime.h"
#endif

// name
namespace Force
{
	struct atmIndex
	{
		int DD, MM, YYYY;
		double F107;
		double F81;
		double aKp;
		int data;
	};
	struct tmp_loadfinals
	{
		int y, m, d;
		double mjd_cur;
		double px, spx, py, spy, dUT1, sdUT1;

		int ip1, ip2;
	};
	class InfluenceForce
	{
	private:
	
		//================================================//
		// InfluencePlanet.cpp
		// ������ ������
		int DE403_clu;
		int DE403_step;
		int DE403_size;
		double *FileMemDE403;
		// �������� ��������
		void DE403LoadFileToMemory();
		void DE403LoadFileToMemory(std::string path_DE403);
		// ������� ������
		void DE403DeleteMemory();
		// Planet Position
		double PLCOORD[11*3];

		// ���������� �������
		void pln_coords( int i, double *pos );
		// ����������� ��� �������
		double mu_plan( int i );
		// ����� ������� ���������
		bool pln_flag( int i );
		// ������� ��������
		void cheb2( int mode, double x, int degree, double *coeff, double *w, double *dw, double *ddw );
		// ���������� ������� �� �������� ��������
		void DE403( double t, int n_pl, int mode, double* x_pl, double ajd0, double delt0 );
		// ���������� ��������� ������ ������������ �����
		void planets_update_geo( double t, double ajd0, double delt0 );
		// ���������� ��������� ������ ������������ ������
		void planets_update_sol( double t, double ajd0, double delt0 );
		// ���������� �� ������
		void planets_grav( double *x, double *f_gr );
		//================================================//

		//================================================//
		//InfluenceSun.cpp
		// ���������� ������ ������
		double shadow_geo( double *r_ka, double *r_sun );
		// ����
		void sp_cannonball( double *x, double *f_sp, double *pln_coords11, double mu_plan11, double sp_q );
		// �������� ������
		void sp_cannonballForce( double *x, double sp_q, double *f_sp );
		//================================================//

		//================================================//
		//InfluenceEGM96.cpp
		double *h_egm96;

		// ��������� �����
		void GetF_Harm_egm96( double *x_b, int n_harm, double *f_harm );
		// ���� ��-�� ��������
		void GetHarmForce( double *x, double *Fharm );
		// init egm
		void InitHarmInCPU();
		//================================================//

		//================================================//
		// InfluenceAtmosphere.cpp
		// ��������� �����
		void Atm_drag( double *x, double t, double *f, double sigma_up, double ajd0, double delt0 );
		// ���� ���������
		//double Roa2004_2( double time, double *x, double ajd0, double delt0 );
		void InitAtm();
		void InitAtm(std::string path_atm_conf);
		std::vector<atmIndex> ListAtmIndex;
		//================================================//

		//================================================//
		// InfluenceTime.cpp
		
		// ������� ������� � ���������� �� �����
		int InitTAU_UTCcorrect();
		int InitTAU_UTCcorrect(std::string path_TAI_UTC);
		// �������� �������
		void DeInitTAU_UTCcorrect();
		// ��������� �������� ������� - ����� �����
		int GetTAU_UTCkmax();
		// ������ �������
		int GetTAU_UTC_size();
		//��������� ����� TDB � UTC
		double tdb_utc( double ajd, double delt, double t );
		// YYYYMMDD.0d0 � HHMMSS.SSSS... � ������� �������������
		void jddeltt( double dt, double tm, double *ajd, double *delt, double *t );
		// �� �������� ������������� � YYYYMMDD.0d0 � HHMMSS.SSSS
		void dttm( double ajd, double delt, double t, double *dt, double *tm );
		//================================================//

		//================================================//
		//InfluenceNutationEarth.cpp
		// ��������� � ������� �����
		//������ �������� ������� ���������
		double E2000( double E1, double E2);
		void FA2000( double AED, double *FA );
		//���������� �������� �������
		void N2000( int N, double AJD, double *HYT );
		//����� ������� �������
		void NM2000( double E, double *HUT );
		//================================================//

		//================================================//
		// InfluiencePrecessionEarth.cpp
		//������ ������� ��������� 
		void PM2000( double E1, double E2, double *P );
		//================================================//

		//================================================//
		// InfluencePoleEarth.cpp
		// ������� ��������
		int finals_n;
		double *finals_tab;

		// ������������� ������
		void iers_init();
		void iers_init(std::string path_finals);
		// �������� ��������
		void iers_delete();
		// ��������� ��������� �������� ��������� ������ ����� � ������
		void get_xyt ( double t, double *xyt, double ajd0, double delt0 );
		// ������� �������� ��������������� �������� �������
		void PM_mat ( double x, double y, double *A_pole );
		//================================================//

		//================================================//
		// InfluenceEarthRotation.cpp
		// ��������� �������� ������� ��������� �� EME2000 � ����������� ��� �������� ��
		void state_to_itrf( double t, double *x_2000, double *x_g, double *A_rot );
		// ��������� ��������� ������� ��������������� ��������� �������� �����
		void ER_mat ( double t, double d_UT1, double *A_rotat, double ajd0, double delt0 );
		// ��������� ��������� ������� �������� �� ITRF � ICRF
		void iers_mat( double t, double *A, double ajd0, double delt0 );
		//================================================//

		//================================================//
		// AdditionalFunction.cpp
		// �������������� �������
		void matMul( double *in1, double *in2, double *out );
		void matVecMul_V6( double *inMat, double *inVec, double *outVec );
		double DMOD( double X, double Y );
		double DDIM( double X, double Y  );
		double DSIGN( double X, double Y );
		double Dmax( double X, double Y );
		double Dmin( double X, double Y );
		//================================================//

		//================================================//
#ifdef GPUCOMPILE
		// GPU init
		void cuDE403LoadFileToGPU();
		void cuSetupAndLoadN2000DataToGPU();
		void cuiers_init_gpu( double *taiutc, int sizeTai, double *table, int sizeTable, int infinals_n );
		void cuLoadegm96ToGPU();
#endif
		//================================================//

	public:
		InfluenceForce();
		~InfluenceForce();

		double current_px_ = 0; //sec.of arc
		double current_py_ = 0; //sec.of arc
		double current_dUt1_ = 0; //sec. of time

		// ������������
		int ST;
		double SIGMA_ATM;
		double SIGMA_SUN;
		// �������� �������
		double S_ajd0;
		double S_delt0;
		// ������������� ������ �� CPU
		void Init_CPU();
		void Init_CPU(std::string path_TAI_UTC, std::string path_DE403, std::string path_finals, std::string path_atmconf);
		void DeInit_CPU();
		// Init koeff
		void SetSigmaAtm( double sg ) { SIGMA_ATM = sg; };
		void SetSigmaSun( double sg ) { SIGMA_SUN = sg; };
		// ������ ��������� �������
		std::vector< tmp_loadfinals > pole_offset;
		
		// �������� �������, ������� ������� � ������ �������
		// ������ � ���������� �� �����
		int TAUUTC_size;
		int TAUUTC_kmax;
		double *TAIUTCCorrect;


		//================================================//
#ifdef GPUCOMPILE
		// GPU
		// ��������� ��� ��������
		double *d_egm96;
				// ��������� ��� �������� �������
		double *d_AMPL;
		double *d_ARG;
		// �������� �������
		double *d_TAIUTCDATA;
		double *d_finals_tab;
		int *d_finals_n;
		// ��������� ������ �� GPU
		double *d_FileMemDE403;
		// ������������� ������ �� CPU
		void Init_GPU();
		void DeInit_GPU();
#endif
		//================================================//

		//================================================//
		// InfluenceTime.cpp
		void set_time( double date, double time, double *ajd0, double *delt0, double *t );
		double get_time( double date, double time, double ajd0, double delt0 );
		double Ajd_dt ( double ajd);
		double dt_ajd ( double dt );
		//================================================//

		//================================================//
		// InfluenceEarthRotation.cpp
		void GetTemeMatrix( double t, double *A, double ajd0, double delt0 );
		void iers_update_matrix( double t, double *A_rot, double ajd0, double delt0 );
		//================================================//

		//================================================//
		// AdditionalFunction.cpp
		void matVecMul( double *inMat, double *inVec, double *outVec );
		void transpose( double *A, double *B );
		//================================================//

		//================================================//
		// InfluenceForce.cpp
		// ���������� ������ ������
		int rh_fun( double t, double *x, double *f );
		int rh_fun_kepler(double t, double *xt, double *fx );
		int rh_fun_grav( double t, double *x, double *f );
		//================================================//

		//================================================//
		// Test
		int TestRotation();
		//================================================//

		//================================================//
		// convert coordinate
		void ConvertXYZtoRADEC(double *posICRF, double *TelICRF, double *Ra, double *Dec);
		void ITRFToICRF(double jd, double *posITRF, double *posICRF);
		void ICRFToITRF(double jd, double *posICRF, double *posITRF);
		//================================================//
	};
};
//==============================================================================//
#endif