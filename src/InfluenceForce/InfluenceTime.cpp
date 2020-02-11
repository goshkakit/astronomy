//==============================================================================//
// Andrianov N.G.
// opbit predict 
// module find Influence Force
// Time
//==============================================================================//
#include <math.h>
#include <cmath>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <sstream>

#include "InfluenceForce.h"

namespace Force
{
	struct tmp_loadtauutc
	{
		double jd;
		double dt;
		double k1;
		double k2;
	};

	// ���������� �� ���������� ������
	int NINTC( double A )
	{
		int ia = (int)A;
		double da = ia;

		if( A > 0 )
		{
			double delta = A-da;
			if( delta >= 0.5 )
				ia = ia+1;
		}
		else
		{
			double delta = A-da;
			if( delta <= -0.5 )
				ia = ia-1;
		}
		return ia;
	}
	//==============================================================================//
	// TDB_UTC : JD, delt, t || TDB - UTC
	// ������� ������� ��������� ����� TDB � UTC,
	// ������������� ������� � �������� TAI � UTC ��-�� ���������� �������� �����.
	//==============================================================================//
	int InfluenceForce::GetTAU_UTCkmax(){
		return TAUUTC_kmax;
	}
	int InfluenceForce::GetTAU_UTC_size(){
		return TAUUTC_size;
	}
	double StdToDouble2( std::string str )
	{
		double res = 0;
		std::istringstream os( str );
		os >> res;
		return res;
	}
	int InfluenceForce::InitTAU_UTCcorrect()
	{
		// list load
		std::vector< tmp_loadtauutc > list;

		// load file
		char *TD_path = "data/TAI-UTC.DAT";
		std::ifstream rfs( TD_path, std::ios::in );
		printf( "LOAD TAU-UTC.DAT\n" );

		while( !rfs.eof() )
		{
			std::string line;
			std::getline( rfs, line );
			if( line.length() < 10 )
				break;

			tmp_loadtauutc pt;

			// ��������� ����
			std::string s_jd = line.substr( 17, 9 );
			pt.jd = StdToDouble2( s_jd );		//!!!!linux

			// ��������
			std::string s_dt = line.substr( 38, 10 );
			pt.dt = StdToDouble2( s_dt );

			// k1
			std::string s_k1 = line.substr( 60, 5 );
			pt.k1 = StdToDouble2( s_k1 );
			// k2
			std::string s_k2 = line.substr( 70, 9 );
			pt.k2 = StdToDouble2( s_k2 );

			list.push_back( pt );
		}
		rfs.close();

		printf("Size TAU-UTC.DAT = %zu\n", list.size() );

		TAUUTC_size = list.size()*4; // 39*4
		TAUUTC_kmax = list.size()-1; // 38

		TAIUTCCorrect = new double[TAUUTC_size];
		for( unsigned int kl = 0; kl < list.size(); kl++ )
		{
			TAIUTCCorrect[ 4*kl + 0 ] = list[kl].jd;
			TAIUTCCorrect[ 4*kl + 1 ] = list[kl].dt;
			TAIUTCCorrect[ 4*kl + 2 ] = list[kl].k1;
			TAIUTCCorrect[ 4*kl + 3 ] = list[kl].k2;
		}
		return 0;
	}
	int InfluenceForce::InitTAU_UTCcorrect(std::string path_TAI_UTC)
	{
		// list load
		std::vector< tmp_loadtauutc > list;

		// load file
		const char *TD_path = path_TAI_UTC.data();
		std::ifstream rfs(TD_path, std::ios::in);
		printf("LOAD TAU-UTC.DAT\n");

		while (!rfs.eof())
		{
			std::string line;
			std::getline(rfs, line);
			if (line.length() < 10)
				break;

			tmp_loadtauutc pt;

			// ��������� ����
			std::string s_jd = line.substr(17, 9);
			pt.jd = StdToDouble2(s_jd);		//!!!!linux

											// ��������
			std::string s_dt = line.substr(38, 10);
			pt.dt = StdToDouble2(s_dt);

			// k1
			std::string s_k1 = line.substr(60, 5);
			pt.k1 = StdToDouble2(s_k1);
			// k2
			std::string s_k2 = line.substr(70, 9);
			pt.k2 = StdToDouble2(s_k2);

			list.push_back(pt);
		}
		rfs.close();

		printf("Size TAU-UTC.DAT = %zu\n", list.size());

		TAUUTC_size = list.size() * 4; // 39*4
		TAUUTC_kmax = list.size() - 1; // 38

		TAIUTCCorrect = new double[TAUUTC_size];
		for (unsigned int kl = 0; kl < list.size(); kl++)
		{
			TAIUTCCorrect[4 * kl + 0] = list[kl].jd;
			TAIUTCCorrect[4 * kl + 1] = list[kl].dt;
			TAIUTCCorrect[4 * kl + 2] = list[kl].k1;
			TAIUTCCorrect[4 * kl + 3] = list[kl].k2;
		}
		return 0;
	}
	//==============================================================================//
	// �������� �������
	//==============================================================================//
	void InfluenceForce::DeInitTAU_UTCcorrect()
	{
		delete TAIUTCCorrect;
	}
	//==============================================================================//
	// TDB_UTC : JD, delt, t || TDB - UTC
	// ������� ������� ��������� ����� TDB � UTC, ���������� �������-
	// ��� �������� � 32.184 ������� � ���������� �������, ��������-
	// ����� � ������� TAI � UTC ��-�� ���������� �������� �����.
	// real*8, intent(in) :: ajd, delt, t
	//==============================================================================//
	double InfluenceForce::tdb_utc( double ajd, double delt, double t )
	{
		double ajdm, jd;
		double tdb_utc_res = 0;

		double *tab = TAIUTCCorrect;
		int kmax = GetTAU_UTCkmax();
		// 1965 MAR   1 = JD 2438820.5  
		// TAI-UTC = 3.6401300 S + (MJD - 38761.) X 0.001296 S
		// ������ ���������� � �������� ������ � ����� TAI-UTC.dat
		//  1961 JAN  1 =JD 2437300.5  TAI-UTC=   1.4228180 S + (MJD - 37300.) X 0.001296 S
		ajdm = ajd+(delt+t)/86.40-2400000.0;
		jd = ajd+(delt+t)/86.40;

		//n4 = 4*kmax-4;
		for( int i = 0; i <= kmax; i++ )
		{
			int i4 = 4*(kmax-i);
			// if (tab[i4+1].le.jd.or.i.eq.kmax)
			if( tab[i4] <= jd || i == kmax)
			{
				tdb_utc_res = tab[i4 + 1] + (ajdm - tab[i4+2])*tab[i4+3] + 32.184;
				//printf( "tab[i4+1] = %f\n", tab[i4+1] );
				return tdb_utc_res;
			}
		}
		tdb_utc_res = tab[1] + (ajd - tab[2])*tab[3] + 32.184;
		return tdb_utc_res;
	}
	//==============================================================================//
	// dt_ajd : date || JD
	// ������� ��������� ��������� ���� �� ����, ���������� � �������
	// YYYYMMDD.0000d0, �.�. ������� ���� �� �������� ����� ������.
	// ���������� ��������� ���� �������� ������� ������.
	//==============================================================================//
	double InfluenceForce::dt_ajd ( double dt )
	{
		double dt_ajd;

		int year, month, day, century, year_c;
		year = floor(dt/10000.0);
		month = floor((dt-year*10000)/100.0);
		day = ((int)(dt))-year*10000-month*100;
		if (month > 2)
		{
			month = month - 3;
		}
		else
		{
			month = month + 9;
			year = year - 1;
		}

		century = year/100;
		year_c = year - 100*century;
		dt_ajd = (146097*century)/4+(1461*year_c)/4+ (153*month+2)/5+day+1721119;

		return  dt_ajd;
	} 
	//==============================================================================//
	// ��������� �������� ���� � ������� ������� YYYYMMDD.0d0 � 
	// HHMMSS.SSSS... � ������� �������������: ��������� ���� ������� JD, 
	// ����� ������������ TDB � �������� delt, � ����� �� ������ ����� �� 
	// ��� � ���. ������ t.
	//
	// double *ajd - ��������� ���� �� ������ �����
	// double *delt - �������� ������� �������� ������� � TDB � ����� ������������ �������� 
	// double *t - ����� �� ������ ����� �� ����������� �������
	//==============================================================================//
	void InfluenceForce::jddeltt( double dt, double tm, double *ajd, double *delt, double *t )
	{
		int h, m, tm_i;
		double s, dutc, tm_frac;

		*ajd = dt_ajd(dt) - 0.5;
		tm_frac = tm - floor(tm);
		tm_i = floor(tm);

		//if (abs(tm_frac-1.d0).le.1.d-8)
		if (fabs(tm_frac-1.0) <= 1.0e-8) //!!!linux
		{
			tm_frac = 0.0;
			tm_i = tm_i+1;
		}
		h = tm_i/10000;
		m = (tm_i - h*10000)/100;
		s = (double)(tm_i - h*10000.0 - m*100.0) + tm_frac;

		double tmp1	= ((double)(h))*3.6 + ((double)(m))*0.06 + s*1.0e-3;

		dutc = tdb_utc( *ajd, -10.8, tmp1 );
		//printf("tdb_utc return %f\n", dutc );

		double tmp2 = -10.8 + dutc/1000.0;

		*delt = tmp2;
		*t = tmp1;
	}
	//==============================================================================//
	// ��������� �������
	// ajd0 - ��������� ���� ������ �������.
	// delt0 - �������� �� 3 ���� ���������� ����������� ������� � ������� TDB � UTC.
	// t - ����� �� ������ ����� � ������� ������
	//
	// �������� ����� �������
	//------------------------------------------------------------------------------//
	// TDB
	// ��� �������������� ���������� � ����� ����� ������� TDB, ��� �� ����������� �����
	// TDB - ����� ������� � ������� ������� � ���������� ��������� �������
	// ���������� �� ����� ������� ����� �� 32.183 �������
	// TDB = TAI + 32.183 - ���������� �������������� ����������
	//------------------------------------------------------------------------------//

	//------------------------------------------------------------------------------//
	// UTC
	// ��� ������ ������� ���������� � UTC
	// ����������� ���������, ������� ���������
	// ������� ������ ����� TAI, �� ������� ��� ������� �� ����� ����� ������ �������
	// � 1968 ����. � 2012 ���� ��� ���������� ��������� 35 ������
	// UTC = TAI + 35.0;
	//------------------------------------------------------------------------------//

	//------------------------------------------------------------------------------//
	// MDB
	// ��������� � � ������ ���������� MDB - ��� �� UTC + 3:00
	// � ������ ������� ���������� ����� MDB
	// ������� ��� � ����������� ����� Tet = t + 32.184 + 35 - 3*3600.
	//------------------------------------------------------------------------------//

	//------------------------------------------------------------------------------//
	// JD
	// ��������� ���� Julian Date
	// J2000 ������������� 12:00 TDB 1 ������ 2000 ����
	// J2000 = 2451545.0
	// ���� ����� ��������� ���� �������� 86400 ������
	// ������ ������������ ����� �� ������ ����� 
	// JD0 = JD - 0.5;
	// ��������� ���� ������� - ����� ����� ����� �� ������� ��� �� ��������
	//------------------------------------------------------------------------------//
	//
	// ajd0 - ��������� ���� �� ������ �����
	// delt0 - �������� ������� �������� ������� � TDB � ����� ������������ �������� � ������� ������
	// t - ����� �� ������ ����� �� ����������� ������� � ������� ������
	// ����� ����������� ����� Tet = ajd0*86.40 + delt0 + t;
	// ��������� ���� �� ����� jd = ajd0 + ( delt0 + t )/86.40 - ���������� � ������
	//==============================================================================//
	void InfluenceForce::set_time( double date, double time, double *ajd0, double *delt0, double *t )
	{
		jddeltt(date, time, ajd0, delt0, t );

		double jd = *ajd0;
		double d = *delt0;
		double mt = *t;
		//printf("set_time: Data = %f, Time = %f\nt = %f,\t jd = %f delta %f, CH = %f\n", date, time, mt, jd, d, mt/3.6 );
	}

	//==============================================================================//
	// ��������� �������
	// ����� �� ������ ����� �� ����������� ������� �������� � double ajd0, double delt0
	//==============================================================================//
	double InfluenceForce::get_time( double date, double time, double ajd0, double delt0 )
	{
		double t;
		double jd, delt;
		bool time_rdy = true;

		if( time_rdy )
		{
			jddeltt( date, time, &jd, &delt, &t );
			t = (jd-ajd0)*86.40 + (delt-delt0) + t;

			//printf("get_time: Data = %f, Time = %f\nt = %f,\t jd = %f delta %f, CH = %f from SetTime\n", date, time, t, jd, delt, t/3.6 );
		}
		else
		{
			printf( " Warnint time module is not initialized! Set t = 0\n" );
			t = 0.0;
		}
		return t;
	}
	//==============================================================================//
	// ajd_dt : JD || date
	// ������� ���������� ����, ���������� � ������� YYYYMMDD.0000d0, 
	// �� ��������� ����.
	//==============================================================================//
	double InfluenceForce::Ajd_dt ( double ajd )
	{
		int y, d, j, m;

		j = NINTC(ajd-1721119.0);
		y = (4*j-1)/146097;
		d = (4*j-1-146097*y)/4;
		j = (4*d+3)/1461;
		d = (4*d+7-1461*j)/4;
		m = (5*d-3)/153;
		d = (5*d+2-153*m)/5;
		y = 100*y+j;
		if (m < 10)
		{
			m=m+3;
		}
		else
		{
			m=m-9;
			y=y+1;
		}
		double ajd_dt=y*10000+m*100+d;
		return ajd_dt;
	}
	//==============================================================================//
	// dttm : JD, delt, t || date, time
	// ��������� �������� �� �������� �������������: ��������� ���� �������
	// JD, ������ ������������ TDB � �������� delt, � ������� �� ������ ��-
	// ��� �� ��� � ���. ������ t -- � ���� � ������� ������� YYYYMMDD.0d0
	// � HHMMSS.SSSS...
	//==============================================================================//
	void  InfluenceForce::dttm( double ajd, double delt, double t, double *dt, double *tm )
	{
		double dutc, res, ajd_i, tt, dtt, ajd_ii, tt_i;
		double tth, hh, am, amm, sec, delt1;
		dutc = tdb_utc(ajd, delt, t);
		res = DMOD(ajd, 1.0);
		ajd_i = ajd - res;

		tt = t + delt + 10.80 + (res + 0.50)*86.40 - dutc/1000.0;
		dtt = DMOD(tt,86.40);

		tt_i = tt - dtt;
		ajd_ii = ajd_i + tt_i/86.4;

		if (dtt< 0.0 ) 
		{
			dtt = 86.40-fabs(dtt); //!!!linux
			ajd_ii = ajd_ii - 1.0;
		}
		else if( fabs(dtt-86.40) < 1.0E-9)  //!!!linux
		{
			dtt = 0.0;
			ajd_ii = ajd_ii+1;
		}

		*dt = Ajd_dt(ajd_ii);

		tth = dtt/3.60;
		int ihh = (int)(tth+1.0E-12);
		hh = ihh;

		am = (tth-hh)*60.0;
		int iamm = (int)(am+1.0E-10);
		amm = iamm;

		sec = (am-amm)*60.0;

		*tm = hh*10000.0 + amm*100.0 + sec;

		delt1 = delt-dutc/1000.0;
	}
	//==============================================================================//
};

//########################################################################//
//double tab[244] = { 2437300.5,  1.422818, 37300.,  .0012960, //  1
//					2437512.5,  1.372818, 37300.,  .0012960, //  2
//					2437665.5,  1.845858, 37665.,  .0011232, //  3
//					2438334.5,  1.945858, 37665.,  .0011232, //  4
//					2438395.5,  3.240130, 38761.,  .0012960, //  5
//					2438486.5,  3.340130, 38761.,  .0012960, //  6
//					2438639.5,  3.440130, 38761.,  .0012960, //  7
//					2438761.5,  3.540130, 38761.,  .0012960, //  8
//					2438820.5,  3.640130, 38761.,  .0012960, //  9
//					2438942.5,  3.740130, 38761.,  .0012960, // 10
//					2439004.5,  3.840130, 38761.,  .0012960, // 11
//					2439126.5,  4.313170, 39126.,  .0025920, // 12
//					2439887.5,  4.213170, 39126.,  .0025920, // 13
//					2441317.5, 10.000000, 41317.,  .0000000, // 14
//					2441499.5, 11.000000, 41317.,  .0000000, // 15
//					2441683.5, 12.000000, 41317.,  .0000000, // 16
//					2442048.5, 13.000000, 41317.,  .0000000, // 17
//					2442413.5, 14.000000, 41317.,  .0000000, // 18
//					2442778.5, 15.000000, 41317.,  .0000000, // 19
//					2443144.5, 16.000000, 41317.,  .0000000, // 20
//					2443509.5, 17.000000, 41317.,  .0000000, // 21
//					2443874.5, 18.000000, 41317.,  .0000000, // 22
//					2444239.5, 19.000000, 41317.,  .0000000, // 23
//					2444786.5, 20.000000, 41317.,  .0000000, // 24
//					2445151.5, 21.000000, 41317.,  .0000000, // 25
//					2445516.5, 22.000000, 41317.,  .0000000, // 26
//					2446247.5, 23.000000, 41317.,  .0000000, // 27
//					2447161.5, 24.000000, 41317.,  .0000000, // 28
//					2447892.5, 25.000000, 41317.,  .0000000, // 29
//					2448257.5, 26.000000, 41317.,  .0000000, // 30
//					2448804.5, 27.000000, 41317.,  .0000000, // 31
//					2449169.5, 28.000000, 41317.,  .0000000, // 32
//					2449534.5, 29.000000, 41317.,  .0000000, // 33
//					2450083.5, 30.000000, 41317.,  .0000000, // 34
//					2450630.5, 31.000000, 41317.,  .0000000, // 35
//					2451179.5, 32.000000, 41317.,  .0000000, // 36
//					2453736.5, 33.000000, 41317.,  .0000000, // 37
//					2454832.5, 34.000000, 41317.,  .0000000, // 38
//					2456109.5, 35.000000, 41317.,  .0000000, // 39
//					.0,   .000000,     0.,  .0000000, // 40
//					.0,   .000000,     0.,  .0000000, // 41
//					.0,   .000000,     0.,  .0000000, // 42
//					.0,   .000000,     0.,  .0000000, // 43
//					.0,   .000000,     0.,  .0000000, // 44
//					.0,   .000000,     0.,  .0000000, // 45
//					.0,   .000000,     0.,  .0000000, // 46
//					.0,   .000000,     0.,  .0000000, // 47
//					.0,   .000000,     0.,  .0000000, // 48
//					.0,   .000000,     0.,  .0000000, // 49
//					.0,   .000000,     0.,  .0000000, // 50
//					.0,   .000000,     0.,  .0000000, // 51
//					.0,   .000000,     0.,  .0000000, // 52
//					.0,   .000000,     0.,  .0000000, // 53
//					.0,   .000000,     0.,  .0000000, // 54
//					.0,   .000000,     0.,  .0000000, // 55
//					.0,   .000000,     0.,  .0000000, // 56
//					.0,   .000000,     0.,  .0000000, // 57
//					.0,   .000000,     0.,  .0000000, // 58
//					.0,   .000000,     0.,  .0000000, // 59
//					.0,   .000000,     0.,  .0000000, // 60
//					.0,   .000000,     0.,  .0000000};// 61

//TAIUTCCorrect = new double[ntu];
//for( int it = 0; it < ntu; it++ )
//{
//	TAIUTCCorrect[it] = tab[it];
//}