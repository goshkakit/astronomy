//==============================================================================//
// Andrianov N.G.
// opbit predict 
// module find Influence Force
// Motion Pole Earth
//==============================================================================//
#include <math.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "InfluenceForce.h"

namespace Force
{
	//==============================================================================//
	// удаление ресурсов
	//==============================================================================//
	void InfluenceForce::iers_delete()
	{
		delete finals_tab;
	};
	//==============================================================================//
	// ѕрцедура инициализации модул€ iers. ¬ пам€ть записываетс€ массив 
	// координат полюса «емли дл€ последующего расчета матрицы поворота из
	// гринвичской —  в инерциальную (EME2000).
	// !!!!
	// Ќадо переделать данный модуль!!!
	// http://data.iers.org/products/7/1378/orig/finals.all
	//==============================================================================//
	//The format of the finals.data, finals.daily, and finals.all files is:
	//
	//Col.#    Format  Quantity
	//-------  ------  -------------------------------------------------------------
	//1-2      I2      year (to get true calendar year, add 1900 for MJD<=51543 or add 2000 for MJD>=51544)
	//3-4      I2      month number
	//5-6      I2      day of month
	//7        X       [blank]
	//8-15     F8.2    fractional Modified Julian Date (MJD UTC)
	//16       X       [blank]
	//17       A1      IERS (I) or Prediction (P) flag for Bull. A polar motion values
	//18       X       [blank]
	//19-27    F9.6    Bull. A PM-x (sec. of arc)
	//28-36    F9.6    error in PM-x (sec. of arc)
	//37       X       [blank]
	//38-46    F9.6    Bull. A PM-y (sec. of arc)
	//47-55    F9.6    error in PM-y (sec. of arc)
	//56-57    2X      [blanks]
	//58       A1      IERS (I) or Prediction (P) flag for Bull. A UT1-UTC values
	//59-68    F10.7   Bull. A UT1-UTC (sec. of time)
	//69-78    F10.7   error in UT1-UTC (sec. of time)
	//79       X       [blank]
	//80-86    F7.4    Bull. A LOD (msec. of time) -- NOT ALWAYS FILLED
	//87-93    F7.4    error in LOD (msec. of time) -- NOT ALWAYS FILLED
	//94-95    2X      [blanks]
	//96       A1      IERS (I) or Prediction (P) flag for Bull. A nutation values
	//97       X       [blank]
	//98-106   F9.3    Bull. A dPSI (msec. of arc)
	//107-115  F9.3    error in dPSI (msec. of arc)
	//116      X       [blank]
	//117-125  F9.3    Bull. A dEPSILON (msec. of arc)
	//126-134  F9.3    error in dEPSILON (msec. of arc)
	//135-144  F10.6   Bull. B PM-x (sec. of arc)
	//145-154  F10.6   Bull. B PM-y (sec. of arc)
	//155-165  F11.7   Bull. B UT1-UTC (sec. of time)
	//166-175  F10.3   Bull. B dPSI (msec. of arc)
	//176-185  F10.3   Bull. B dEPSILON (msec. of arc)
	double StdToDouble( std::string str )
	{
		double res = 0;
		std::istringstream os( str );
		os >> res;
		return res;
	}
	int StdToInt( std::string str )
	{
		int res = 0;
		std::istringstream os( str );
		os >> res;
		return res;
	}

	void InfluenceForce::iers_init()
	{
		//------------------------------------//
		// load file
		char *TD_path = "data/finals.all";
		std::ifstream rfs( TD_path, std::ios::in );
		printf( "LOAD finals.all\n" );

		
		int kl = 0;
		while( !rfs.eof() )
		{
			kl++;
			//printf( "k = %d\r", kl );
			std::string line;
			std::getline( rfs, line );


			if( line.length() < 70 )
				break;
			try
			{
				tmp_loadfinals pt;

				// data
				std::string s_y = line.substr( 0, 2 );
				pt.y = StdToInt( s_y );

				std::string s_m = line.substr( 2, 2 );
				pt.m = StdToInt( s_m );
				
				std::string s_d = line.substr( 4, 2 );
				pt.d = StdToInt( s_d );

				// j data
				std::string s_j = line.substr( 7, 8 );
				pt.mjd_cur = StdToDouble( s_j );
	
				// ip1
				std::string s_ip1 = line.substr( 16, 1 );
				if( s_ip1 == "I" )
					pt.ip1 = 1;
				else
					pt.ip1 = 2;
				// ip2
				std::string s_ip2 = line.substr( 57, 1 );
				if( s_ip2 == "I" )
					pt.ip2 = 1;
				else
					pt.ip2 = 2;

				// 28 8; 37 9;  47 8;  58 10;  69 9
				// x
				std::string s_px = line.substr( 18, 9 ); 
				pt.px = StdToDouble( s_px );
				std::string s_spx = line.substr( 28, 8 ); 
				pt.spx = StdToDouble( s_spx );

				// y
				std::string s_py = line.substr( 37, 9 ); 
				pt.py = StdToDouble( s_py );
				std::string s_spy = line.substr( 47, 8 ); 
				pt.spy = StdToDouble( s_spy );

				// t
				std::string s_t1 = line.substr( 58, 10 ); 
				pt.dUT1 = StdToDouble( s_t1 );
				std::string s_t2 = line.substr( 69, 9 ); 
				pt.sdUT1 = StdToDouble( s_t2 );

				pole_offset.push_back( pt );
			}
			catch( ... )
			{
				printf("End, error convert\n" );
				break;
			}
		}
		rfs.close();

		printf( "load size = %zu\n", pole_offset.size() );
		//FILE *flist = fopen( "data/tmp_fall.txt", "w" );
		//for( unsigned int it = 0; it < list.size(); it++ )
		//{
		//	int ip1 = 1;
		//	int ip2 = 1;
		//	fprintf( flist, "%d %d %d %f %d %f %f %f %f %d %f %f\n", list[it].y, list[it].m, list[it].d, list[it].mjd_cur, ip1, list[it].px, list[it].spx, list[it].py, list[it].spy, ip2, list[it].dUT1, list[it].sdUT1 );
		//}
		//fclose( flist );
		//------------------------------------//

		// указатель на массив поправки времени
		double *tab = TAIUTCCorrect;
		// временный массив
		double *utc = new double[40000];
		// индекс максимальной записи 0 1 2 .. kmax
		int tsize = GetTAU_UTC_size();
		// массив с поправками
		finals_tab = new double[260000];
		
		// вычислени€
		double ajd50 = 2433282.50;
		double mjd0 = 2400000.50;
		double radsec = 206264.8062470970;
		int ierr, i, j, jend, j_fin, k, un;
		double jd_cur;
		double b, t0, ak, jd_b, jd_e, tdt_UTC, UT1_tdt;

		j = 1;
		for( int it = 0; it < tsize; it++ )
		{
			utc[it+1] = tab[it];
			j++;

		}
		j = j/4;

		jend = j;
		jd_b = utc[1];
		jd_e = utc[(jend-1)*4+1];
		printf("j %d  jd_b %f  jd_e   %f\n (jend-1)*4+1 = %d\n", j, jd_b, jd_e, (jend-1)*4+1);

		double mjd_cur;
		double dUT1;

		j = 1;
		for( unsigned int Klinei = 0; Klinei < pole_offset.size(); Klinei++ )
		{
			// get param
			//int y, m, d;
			int ip1 = pole_offset[Klinei].ip1;
			int ip2 = pole_offset[Klinei].ip2;

			mjd_cur = pole_offset[Klinei].mjd_cur;
			dUT1 = pole_offset[Klinei].dUT1;

			double px = pole_offset[Klinei].px;
			double spx = pole_offset[Klinei].spx;
			double py = pole_offset[Klinei].py;
			double spy = pole_offset[Klinei].spy;
			double sdUT1 = pole_offset[Klinei].sdUT1;
			
			i = (j-1)*5;
			jd_cur = mjd_cur + mjd0;
			if (jd_cur < jd_b)
			{
				b  = utc[2];
				t0 = utc[3];
				ak = utc[4];
			}
			else if (jd_cur > jd_e) 
			{
				b  = utc[(jend-1)*4+2];
				t0 = utc[(jend-1)*4+3];
				ak = utc[(jend-1)*4+4];
			}
			else
			{
				for( k = 1; k <= jend-1; k++ )
				{
					if ((jd_cur >= utc[(k-1)*4+1]) && (jd_cur < utc[k*4+1]) ) 
					{
						b  = utc[(k-1)*4+2];
						t0 = utc[(k-1)*4+3];
						ak = utc[(k-1)*4+4];
					}
				}
			}

			double tdt_utc = b + (jd_cur-t0-mjd0)*ak + 32.1840;
			finals_tab[i+1] = mjd_cur;

			if ((ip1 != 1)&&(ip1 != 2 )) 
			{
				px   = finals_tab[i-3];
				py   = finals_tab[i-2];
				dUT1 = finals_tab[i-1];
			}

			finals_tab[i+2] = px/radsec;
			finals_tab[i+3] = py/radsec;
			finals_tab[i+4] = dUT1/1000.0;

			UT1_tdt = (dUT1-tdt_utc)/1000.0;

			finals_tab[i+5] = -UT1_tdt*1000.0;

			j = j+1;
		}

		// дополнение таблицы нул€ми 
		finals_n = j;
		j_fin = j;
		for( j = j_fin; j <= 20000; j++ )
		{
			i = (j-1)*5;
			finals_tab[i+1] = mjd_cur+(double)(j_fin-j+1);
			finals_tab[i+2] = 0.0;
			finals_tab[i+3] = 0.0;
			finals_tab[i+4] = dUT1/1000.0;
			finals_tab[i+4] = UT1_tdt*1000.0;
		}
		printf("finals_n = %d\n", finals_n );
		delete utc;
	}
	void InfluenceForce::iers_init(std::string path_finals)
	{
		//------------------------------------//
		// load file
		const char *TD_path = path_finals.data();
		std::ifstream rfs(TD_path, std::ios::in);
		printf("LOAD finals.all\n");


		int kl = 0;
		while (!rfs.eof())
		{
			kl++;
			//printf( "k = %d\r", kl );
			std::string line;
			std::getline(rfs, line);


			if (line.length() < 70)
				break;
			try
			{
				tmp_loadfinals pt;

				// data
				std::string s_y = line.substr(0, 2);
				pt.y = StdToInt(s_y);

				std::string s_m = line.substr(2, 2);
				pt.m = StdToInt(s_m);

				std::string s_d = line.substr(4, 2);
				pt.d = StdToInt(s_d);

				// j data
				std::string s_j = line.substr(7, 8);
				pt.mjd_cur = StdToDouble(s_j);

				// ip1
				std::string s_ip1 = line.substr(16, 1);
				if (s_ip1 == "I")
					pt.ip1 = 1;
				else
					pt.ip1 = 2;
				// ip2
				std::string s_ip2 = line.substr(57, 1);
				if (s_ip2 == "I")
					pt.ip2 = 1;
				else
					pt.ip2 = 2;

				// 28 8; 37 9;  47 8;  58 10;  69 9
				// x
				std::string s_px = line.substr(18, 9);
				pt.px = StdToDouble(s_px);
				std::string s_spx = line.substr(28, 8);
				pt.spx = StdToDouble(s_spx);

				// y
				std::string s_py = line.substr(37, 9);
				pt.py = StdToDouble(s_py);
				std::string s_spy = line.substr(47, 8);
				pt.spy = StdToDouble(s_spy);

				// t
				std::string s_t1 = line.substr(58, 10);
				pt.dUT1 = StdToDouble(s_t1);
				std::string s_t2 = line.substr(69, 9);
				pt.sdUT1 = StdToDouble(s_t2);

				pole_offset.push_back(pt);
			}
			catch (...)
			{
				printf("End, error convert\n");
				break;
			}
		}
		rfs.close();

		printf("load size = %zu\n", pole_offset.size());
		//FILE *flist = fopen( "data/tmp_fall.txt", "w" );
		//for( unsigned int it = 0; it < list.size(); it++ )
		//{
		//	int ip1 = 1;
		//	int ip2 = 1;
		//	fprintf( flist, "%d %d %d %f %d %f %f %f %f %d %f %f\n", list[it].y, list[it].m, list[it].d, list[it].mjd_cur, ip1, list[it].px, list[it].spx, list[it].py, list[it].spy, ip2, list[it].dUT1, list[it].sdUT1 );
		//}
		//fclose( flist );
		//------------------------------------//

		// указатель на массив поправки времени
		double *tab = TAIUTCCorrect;
		// временный массив
		double *utc = new double[40000];
		// индекс максимальной записи 0 1 2 .. kmax
		int tsize = GetTAU_UTC_size();
		// массив с поправками
		finals_tab = new double[260000];

		// вычислени€
		double ajd50 = 2433282.50;
		double mjd0 = 2400000.50;
		double radsec = 206264.8062470970;
		int ierr, i, j, jend, j_fin, k, un;
		double jd_cur;
		double b, t0, ak, jd_b, jd_e, tdt_UTC, UT1_tdt;

		j = 1;
		for (int it = 0; it < tsize; it++)
		{
			utc[it + 1] = tab[it];
			j++;

		}
		j = j / 4;

		jend = j;
		jd_b = utc[1];
		jd_e = utc[(jend - 1) * 4 + 1];
		printf("j %d  jd_b %f  jd_e   %f\n (jend-1)*4+1 = %d\n", j, jd_b, jd_e, (jend - 1) * 4 + 1);

		double mjd_cur;
		double dUT1;

		j = 1;
		for (unsigned int Klinei = 0; Klinei < pole_offset.size(); Klinei++)
		{
			// get param
			//int y, m, d;
			int ip1 = pole_offset[Klinei].ip1;
			int ip2 = pole_offset[Klinei].ip2;

			mjd_cur = pole_offset[Klinei].mjd_cur;
			dUT1 = pole_offset[Klinei].dUT1;

			double px = pole_offset[Klinei].px;
			double spx = pole_offset[Klinei].spx;
			double py = pole_offset[Klinei].py;
			double spy = pole_offset[Klinei].spy;
			double sdUT1 = pole_offset[Klinei].sdUT1;

			i = (j - 1) * 5;
			jd_cur = mjd_cur + mjd0;
			if (jd_cur < jd_b)
			{
				b = utc[2];
				t0 = utc[3];
				ak = utc[4];
			}
			else if (jd_cur > jd_e)
			{
				b = utc[(jend - 1) * 4 + 2];
				t0 = utc[(jend - 1) * 4 + 3];
				ak = utc[(jend - 1) * 4 + 4];
			}
			else
			{
				for (k = 1; k <= jend - 1; k++)
				{
					if ((jd_cur >= utc[(k - 1) * 4 + 1]) && (jd_cur < utc[k * 4 + 1]))
					{
						b = utc[(k - 1) * 4 + 2];
						t0 = utc[(k - 1) * 4 + 3];
						ak = utc[(k - 1) * 4 + 4];
					}
				}
			}

			double tdt_utc = b + (jd_cur - t0 - mjd0)*ak + 32.1840;
			finals_tab[i + 1] = mjd_cur;

			if ((ip1 != 1) && (ip1 != 2))
			{
				px = finals_tab[i - 3];
				py = finals_tab[i - 2];
				dUT1 = finals_tab[i - 1];
			}

			finals_tab[i + 2] = px / radsec;
			finals_tab[i + 3] = py / radsec;
			finals_tab[i + 4] = dUT1 / 1000.0;

			UT1_tdt = (dUT1 - tdt_utc) / 1000.0;

			finals_tab[i + 5] = -UT1_tdt*1000.0;

			j = j + 1;
		}

		// дополнение таблицы нул€ми 
		finals_n = j;
		j_fin = j;
		for (j = j_fin; j <= 20000; j++)
		{
			i = (j - 1) * 5;
			finals_tab[i + 1] = mjd_cur + (double)(j_fin - j + 1);
			finals_tab[i + 2] = 0.0;
			finals_tab[i + 3] = 0.0;
			finals_tab[i + 4] = dUT1 / 1000.0;
			finals_tab[i + 4] = UT1_tdt*1000.0;
		}
		printf("finals_n = %d\n", finals_n);
		delete utc;
	}
	//==============================================================================//
	// ѕроцедура получени€ значени€ координат полюса «емли и сдвига по вре-
	// мени в момент t.
	//==============================================================================//
	void InfluenceForce::get_xyt ( double t, double *xyt, double ajd0, double delt0 )
	{
		double jd0000 = 2400000.0;
		double xyt4[4];
		double jd, jd00, res;
		int k, k5;

		xyt[0] = 0.0;
		xyt[1] = 0.0;
		xyt[2] = 0.0;

		jd   = ajd0+(delt0+t)/86.40;
		jd00 = finals_tab[1] + jd0000;

		if ((jd < 23.0E+5) || (jd > 25.0E+5)) 
		{
			printf( "Module IERS -> get_xyt : unexpected Julian date. Return\n" );
			return;
		}
		k = int(jd-jd00)+1;
		if (k <= 0) 
		{
			printf( "jd = %f jd00 = %f\n", jd, jd00 );
			printf( "Module IERS -> get_xtr : Time is less than first record of finals.dat table. Return\n" );
		}
		if (k > finals_n-2) 
		{
			k = finals_n-2;
		}
		else if (k < 1) 
		{
			k = 1;
		}
		res = jd - ( jd00 + (double)(k-1) );
		//res = jd - ( jd00 + (double)(k) );
		k5 = (k-1)*5;
		//k5 = (k)*5;

		//xyt4 = finals_tab(k5+2:k5+5)*(1.d0-res)+finals_tab(k5+7:k5+10)*res
		xyt4[0] = finals_tab[k5+2]*(1.0-res) + finals_tab[k5+7]*res;
		xyt4[1] = finals_tab[k5+3]*(1.0-res) + finals_tab[k5+8]*res;
		xyt4[2] = finals_tab[k5+4]*(1.0-res) + finals_tab[k5+9]*res;
		xyt4[3] = finals_tab[k5+5]*(1.0-res) + finals_tab[k5+10]*res;

		//xyt = xyt4([1, 2, 4])
		xyt[0] = xyt4[0];
		xyt[1] = xyt4[1];
		xyt[2] = xyt4[3];

		//printf( "get_xyt = %f %f %f\n", xyt[0], xyt[1], xyt[2] );
	}
	//==============================================================================//
	// ѕроцедура расчета матрицы поворота соответствующей смещению полюсов
	//==============================================================================//
	void InfluenceForce::PM_mat ( double x, double y, double *A_pole )
	{
		double a, qx, qy, gamma, x2y2;
		qx   = x*x;
		qy   = y*y;
		x2y2 = qx+qy;
		if( (x2y2) < 1.0E-18)
		{
			// —мещение полюсов нулевое или почти нулевое:
			// искома€ матрица - единична€
			for( int it = 0; it < 9; it++ )
				A_pole[it] = 0.0;
			A_pole[0] = 1.0;
			A_pole[4] = 1.0;
			A_pole[8] = 1.0;
		}
		else
		{
			// «адание матрицы, соответствующей смещению полюсов
			a     = sqrt(1.0-x2y2);
			gamma = -x*y*(a-1.0)/(x2y2);
			A_pole[0] = (qy+qx*a)/(x2y2);
			A_pole[4] = (qx+qy*a)/(x2y2);
			A_pole[1] = gamma;
			A_pole[3] = gamma;
			A_pole[8] = a;

			A_pole[2] = x;
			A_pole[6] = -x;
			A_pole[5] = -y;
			A_pole[7] = y;
		}
	}
	//==============================================================================//
};