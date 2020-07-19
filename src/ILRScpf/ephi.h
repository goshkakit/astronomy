#ifndef EPHIHEADER_H 
#define EPHIHEADER_H

#define MAXENT 36000

//---------------------------------------------------------------------------------//
// заголовок
//---------------------------------------------------------------------------------//
struct eh {
	double prf;
	double xpond_xmit_delay;
	double xpond_utc_off;
	double xpond_osc_drift;
	double cofm_corr;
	int cospar_id;
	int sic;
	int norad;
	int eph_seq;
	int sttpyear;
	int sttpmonth;
	int sttpday;
	int sttphour;
	int sttpmin;
	int sttpsec;
	int endpyear;
	int endpmonth;
	int endpday;
	int endphour;
	int endpmin;
	int endpsec;
	int ephsep;
	int compat;
	int version;
	int epyear;
	int epmon;
	int epday;
	int ephour;
	int tar_type;
	int ref_frame;
	int rot_type;
	int cofm_app;
	int ltcorr;
	int atrk_ro_0;
	int xtrk_ro_0;
	int rtrk_ro_0;
	int atrk_ro_1;
	int xtrk_ro_1;
	int rtrk_ro_1;
	int atrk_ro_2;
	int xtrk_ro_2;
	int rtrk_ro_2;
	char eph_source[4];
	char tar_name[10];
	char eph_notes[10];
};
//---------------------------------------------------------------------------------//
// эфемерида
//---------------------------------------------------------------------------------//
struct ei {
	double jdi[MAXENT];
	double jdf[MAXENT];
	double outv[6][MAXENT];
	double backv[6][MAXENT];
	double c12[MAXENT];
	double c23[MAXENT];
	double aberoutv[3][MAXENT];
	double aberbackv[3][MAXENT];
	int leapflag[MAXENT];
	double jds[MAXENT];
	double oscrelcorr[MAXENT];
	double offsetjdi[MAXENT];
	double offsetjdf[MAXENT];
	char offsetname[11][MAXENT];
	double offsetoutv[3][MAXENT];
	double offsetbackv[3][MAXENT];
	double offsetjds[MAXENT];
	double rotjdi[MAXENT];
	double rotjdf[MAXENT];
	double rotangv[3][MAXENT];
	double gast[MAXENT];
	double rotjds[MAXENT];
	double eopjdi[MAXENT];
	double eopjdf[MAXENT];
	double eopxy[2][MAXENT];
	double eopdutc[MAXENT];
	double eopjds[MAXENT];
	int	nv;
	int  nvoff;
	int  nvrot;
	int  nveop;
} ;

//---------------------------------------------------------------------------------//
//
//---------------------------------------------------------------------------------//
//struct ei {
//	double* jdi;//[MAXENT];
//	double* jdf;//[MAXENT];
//
//	double* outv;//[6][MAXENT];
//	double* backv;//[6][MAXENT];
//
//	double* c12;//[MAXENT];
//	double* c23;//[MAXENT];
//
//	double* aberoutv;//[3][MAXENT];
//	double* aberbackv;//[3][MAXENT];
//
//	int* leapflag;//[MAXENT];
//	double* jds;//[MAXENT];
//	double* oscrelcorr;//[MAXENT];
//	double* offsetjdi;//[MAXENT];
//	double* offsetjdf;//[MAXENT];
//
//	char* offsetname;//[11][MAXENT];
//	double* offsetoutv;//[3][MAXENT];
//	double* offsetbackv;//[3][MAXENT];
//
//	double* offsetjds;//[MAXENT];
//	double* rotjdi;//[MAXENT];
//	double* rotjdf;//[MAXENT];
//
//	double* rotangv;//[3][MAXENT];
//
//	double* gast;//[MAXENT];
//	double* rotjds;//[MAXENT];
//	double* eopjdi;//[MAXENT];
//	double* eopjdf;//[MAXENT];
//
//	double* eopxy;//[2][MAXENT];
//
//	double* eopdutc;//[MAXENT];
//	double* eopjds;//[MAXENT];
//
//	int	nv;
//	int  nvoff;
//	int  nvrot;
//	int  nveop;
//
//	void InitMem()
//	{
//		jdi = new double[MAXENT];
//		jdf = new double[MAXENT];
//		// 2
//		outv = new double[6*MAXENT];
//		backv = new double[6*MAXENT];
//		
//		c12 = new double[MAXENT];
//		c23 = new double[MAXENT];
//
//		// 2
//		aberoutv = new double[3*MAXENT];
//		aberbackv = new double[3*MAXENT];
//
//		leapflag = new int [MAXENT];
//		jds = new double[MAXENT];
//		oscrelcorr = new double[MAXENT];
//		offsetjdi = new double[MAXENT];
//		offsetjdf = new double[MAXENT];
//
//		// 3
//		offsetname = new char[11*MAXENT];
//		offsetoutv = new double[3*MAXENT];
//		offsetbackv = new double[3*MAXENT];
//
//		offsetjds = new double[MAXENT];
//		rotjdi = new double[MAXENT];
//		rotjdf = new double[MAXENT];
//
//		// 1
//		rotangv = new double[3*MAXENT];
//
//		gast = new double[MAXENT];
//		rotjds = new double[MAXENT];
//		eopjdi = new double[MAXENT];
//		eopjdf = new double[MAXENT];
//
//		// 1
//		eopxy = new double[2*MAXENT];
//
//		eopdutc = new double[MAXENT];
//		eopjds = new double[MAXENT];
//	}
//};
//---------------------------------------------------------------------------------//
#endif
