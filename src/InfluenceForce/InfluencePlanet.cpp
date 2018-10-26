//==============================================================================//
// Andrianov N.G.
// opbit predict 
// module find Influence Force
// find force from planet
//==============================================================================//
#include <math.h>
#include <fstream>
#include <stdio.h>
#include <string.h>

#include "InfluenceForce.h"

namespace Force
{
	//==============================================================================//
	// координаты планет
	//==============================================================================//
	void InfluenceForce::pln_coords( int i, double *pos )
	{
		pos[0] = PLCOORD[i*3 + 0];
		pos[1] = PLCOORD[i*3 + 1];
		pos[2] = PLCOORD[i*3 + 2];
	}
	//==============================================================================//
	// гравитационные параметры планет
	//==============================================================================//
	double InfluenceForce::mu_plan( int i )
	{
		double mu_plan[11] = {
			22.03208047245,		// Меркурий
			324.8587656142,		// Венера
			398.6004415,		// Земля
			42.828287,			// Марс
			126712.59708,		// Юпитер
			37939.51971,		// Сатурн
			5780.158533417,		// Уран
			6871.307771094,		// Нептун
			1.02086,			// Плутон
			4.90279914,			// Луна
			132712439.935 };	// Солнце

			return mu_plan[i];
	}
	//==============================================================================//
	// учитывать ли гравитацию от планет
	// Меркурий
	// Венера
	// Земля
	// Марс
	// Юпитер
	// Сатурн
	// Уран
	// Нептун
	// Плутон
	// Луна
	// Солнце
	//==============================================================================//
	bool InfluenceForce::pln_flag( int i )
	{
		int FLFlag[11] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
		return 1;
	}
	//==============================================================================//
	// координаты планет из полинома чебышева
	//==============================================================================//
	void InfluenceForce::cheb2( int mode, double x, int degree, double *coeff, double *w, double *dw, double *ddw )
	{
		double bk, bk1, bk2, dbk, dbk1, dbk2, ddbk, ddbk1, ddbk2;

		switch( mode ) 
		{
		case 1:
			bk1 = 0;
			bk2 = 0;
			//do i = degree,1,-1
			for( int i = degree; i > 0; i-- )
			{
				bk = coeff[i] - bk2 + 2.0*bk1*x;
				bk2 = bk1;
				bk1 = bk;
			}
			*w = coeff[0] - bk2 + bk1*x;
			//*dw = 0.0;
			//*ddw = 0.0;
			break; 

		case 2:
			bk1 = 0.0;
			bk2 = 0.0;
			dbk1 = 0.0;
			dbk2 = 0.0;
			//do i = degree,1,-1
			for( int i = degree; i > 0; i-- )
			{
				bk = coeff[i]-bk2 + 2.0*bk1*x;
				bk2 = bk1;
				bk1 = bk;
				dbk = 2.0*(dbk1*x + bk2) - dbk2;
				dbk2 = dbk1;
				dbk1 = dbk;
			}
			*w = coeff[0] - bk2 + bk1*x;
			*dw = bk1 + dbk1*x - dbk2;
			//*ddw = 0.0;
			break; 

		case 3:
			bk1 = 0.0;
			bk2 = 0.0;
			dbk1 = 0.0;
			dbk2 = 0.0;
			ddbk1 = 0.0;
			ddbk2 = 0.0;
			//do i = degree,1,-1
			for( int i = degree; i > 0; i-- )
			{
				bk = coeff[i] - bk2 + 2.0*bk1*x;
				bk2 = bk1;
				bk1 = bk;
				dbk = 2.0*(dbk1*x + bk2) - dbk2;
				dbk2 = dbk1;
				dbk1 = dbk;
				ddbk = 2.0*(ddbk1*x + 2.0*dbk2) - ddbk2;

				ddbk2 = ddbk1;
				ddbk1 = ddbk;
			}
			*w = coeff[0] - bk2 + bk1*x;
			*dw = bk1 + dbk1*x - dbk2;
			*ddw = 2.0*dbk1 + ddbk1*x - ddbk2;
			break; 
		}
	};
	//==============================================================================//
	// интерполяция орбиты планет и вычисление координат на заданный момент времени
	//==============================================================================//
	void InfluenceForce::DE403( double t, int n_pl, int mode, double* x_pl, double ajd0, double delt0 )
	{
		double jd50 = 2433282.5;
		double day_begin = -54770.0;
		double day_end = 59726.0;
		double span = 32.0;
		double daysec = 2.314814814815e-5;
		double buffer[746];

		int c_pl;
		int j;
		int knot;
		int degree;
		double dif;
		double t_span;
		double arg;
		double day;

		double intp[11] = {8.0, 32.0, 16.0, 32.0, 32.0, 32.0, 32.0, 32.0, 32.0, 4.0, 32.0 };
		int p[11] = { 3, 147, 183, 273, 303, 330, 354, 378, 396, 414, 702 };
		int c[11] = { 12, 12, 15, 10, 9, 8, 8, 6, 6, 12, 15 };

		int first_rec = 1;
		// время в юлианской дате в системе TDB - эфемиридное время
		day = ajd0-jd50 + (t + delt0)/86.4;
		//printf("input day = %f\n", day );
		// span - шаг блока даных об эфемеридах
		int rec_num = (int)((day-day_begin)/span);
		// номер блока
		int clu = rec_num + first_rec;
		//printf( "clu = %d\n", clu );

		//char *DE_path = "data\\DE403.bin";
		//std::ifstream rfs( DE_path, std::ios::in | std::ios::binary);
		//for( int it = 0; it < clu; it++ )
		//	rfs.read( (char *)buffer, 746*sizeof( double ) );
		//rfs.close();
		memcpy( buffer, &FileMemDE403[746*(clu-1)], 746*sizeof( double ) );

		//nint
		dif = (int)(2.0*(ajd0-jd50-buffer[1])/2.0)+(t+delt0)/86.40;
		// число коэффициентов для планеты
		c_pl = c[n_pl];
		// степень полинома
		degree = c_pl-1;
		// длительность полинома
		t_span = intp[n_pl];
		// номер подблока
		knot = (int)(dif/t_span);
		// перевод аргумента к интервалу интерполяции [-1 1 ] 
		arg = 2.0*(dif/t_span-knot)-1.0;
		// индекс начала подблока
		// 3* - три координаты
		j = p[n_pl] - 1 + 3*knot*c_pl;
		//printf("Pos index J = %d\n", j );

		t_span = daysec/t_span;
		switch( mode ) 
		{
		case 1:
			//do i = 1,3
			for( int i = 0; i < 3; i++ )
			{
				cheb2(1, arg, degree, &buffer[j+i*c_pl], &x_pl[i], NULL, NULL );
				x_pl[i] = x_pl[i]*1.e-3;
			}
			break; 
		case 2:
			//do i = 1,3
			for( int i = 0; i < 3; i++ )
			{
				cheb2(2, arg, degree, &buffer[j+i*c_pl], &x_pl[i], &x_pl[i+3], NULL);
				x_pl[i] = x_pl[i]*1.0e-3;
				x_pl[i+3] = x_pl[i+3]*t_span;
			}
			break; 
		case 3:
			//do i = 1,3
			for( int i = 0; i < 3; i++ )
			{
				cheb2(3, arg, degree, &buffer[j+i*c_pl], &x_pl[i], &x_pl[i+3], &x_pl[i+6] );
				x_pl[i] = x_pl[i]*1.0e-3;
				x_pl[i+3] = x_pl[i+3]*t_span;
				x_pl[i+6] = x_pl[i+6]*t_span*t_span*1.0e+3;
			}
			break; 
		}
	}
	//==============================================================================//
	// положение планет
	// важен размер файла
	// необходимо сдедить за изменениями
	//==============================================================================//
	void InfluenceForce::DE403LoadFileToMemory()
	{
		printf("Load DE403\n");
		// 21 347 536 byte
		// 2668442 double
		// 3577 block
		DE403_clu = 3577;
		DE403_step = 746;
		DE403_size = 2668442;
		FileMemDE403 = new double[ DE403_size + 10000 ];

		char *DE_path = "data/DE403.bin";
		std::ifstream rfs( DE_path, std::ios::in | std::ios::binary);
		for( int it = 0; it < DE403_clu; it++ )
		{
			double *ptr = &FileMemDE403[ DE403_step*it];
			rfs.read( (char *)ptr, DE403_step*sizeof( double ) );
		}
		rfs.close();
	}
	//==============================================================================//
	// clear all
	//==============================================================================//
	void InfluenceForce::DE403DeleteMemory()
	{
		printf("Delete FileMemDE403\n");
		delete FileMemDE403;
	}
	//==============================================================================//
	// Процедура обновления массива координат планет.
	// Координаты записываются в ГЕОЦЕНТРИЧЕСКОЙ СК
	//==============================================================================//
	void InfluenceForce::planets_update_geo( double t, double ajd0, double delt0 )
	{
		double EM_r[3];
		double M_r[3];

		DE403( t, 2, 1, EM_r, ajd0, delt0 );	// Земля+Луна
		DE403( t, 9, 1, M_r, ajd0, delt0 );	// Луна (геоцентр)

		//pln_coords(10,:) = M_r
		for( int it = 0; it < 3; it++ )				// Луна
			PLCOORD[9*3+it]= M_r[it];

		//pln_coords(3,:) = EM_r-M_r/M_ME			// Земля 
		double M_ME = 82.300578;
		double POS_E[3];
		for( int it = 0; it < 3; it++ )
		{
			PLCOORD[2*3+it]= EM_r[it] - M_r[it]/M_ME;
			// положение земли относительно солнца
			POS_E[it] = PLCOORD[2*3+it];
		}

		DE403(t, 0, 1, &PLCOORD[0*3], ajd0, delt0 ); // Меркурий
		DE403(t, 1, 1, &PLCOORD[1*3], ajd0, delt0 ); // Венера
		DE403(t, 3, 1, &PLCOORD[3*3], ajd0, delt0 ); // Марс
		DE403(t, 4, 1, &PLCOORD[4*3], ajd0, delt0 ); // Юпитер
		DE403(t, 5, 1, &PLCOORD[5*3], ajd0, delt0 ); // Сатурн
		DE403(t, 6, 1, &PLCOORD[6*3], ajd0, delt0 ); // Уран
		DE403(t, 7, 1, &PLCOORD[7*3], ajd0, delt0 ); // Нептун
		DE403(t, 8, 1, &PLCOORD[8*3], ajd0, delt0 ); // Плутон

		// положение солнца это - положение земли
		for( int it = 0; it < 3; it++ )
			PLCOORD[10*3+it] = -PLCOORD[2*3+it];		// Солнце

		// переводим все планеты состему координат относительно земли
		// 9 - луна
		// 10 - солнце
		for( int i = 0; i < 9; i++ )
		{
			PLCOORD[i*3+0] = PLCOORD[i*3+0] - POS_E[0];
			PLCOORD[i*3+1] = PLCOORD[i*3+1] - POS_E[1];
			PLCOORD[i*3+2] = PLCOORD[i*3+2] - POS_E[2];
		}
	}
	//==============================================================================//
	// Процедура обновления массива координат планет.
	// Координаты записываются в ГЕЛИОЦЕНТРИЧЕСКОЙ СК
	//==============================================================================//
	void InfluenceForce::planets_update_sol( double t, double ajd0, double delt0 )
	{
		double EM_r[3];
		double M_r[3];

		DE403( t, 2, 1, EM_r, ajd0, delt0);	// Земля+Луна
		DE403( t, 9, 1, M_r, ajd0, delt0);	// Луна (геоцентр)

		//pln_coords(10,:) = M_r
		for( int it = 0; it < 3; it++ )				// Луна
			PLCOORD[9*3+it]= M_r[it];

		//pln_coords(3,:) = EM_r-M_r/M_ME ! Земля 
		double M_ME = 82.300578;
		for( int it = 0; it < 3; it++ )
			PLCOORD[2*3+it]= EM_r[it] - M_r[it]/M_ME;

		DE403(t, 0, 1, &PLCOORD[0*3], ajd0, delt0 ); // Меркурий
		DE403(t, 1, 1, &PLCOORD[1*3], ajd0, delt0 ); // Венера
		DE403(t, 3, 1, &PLCOORD[3*3], ajd0, delt0 ); // Марс
		DE403(t, 4, 1, &PLCOORD[4*3], ajd0, delt0 ); // Юпитер
		DE403(t, 5, 1, &PLCOORD[5*3], ajd0, delt0 ); // Сатурн
		DE403(t, 6, 1, &PLCOORD[6*3], ajd0, delt0 ); // Уран
		DE403(t, 7, 1, &PLCOORD[7*3], ajd0, delt0 ); // Нептун
		DE403(t, 8, 1, &PLCOORD[8*3], ajd0, delt0 ); // Плутон

		//pln_coords(10,:) = pln_coords(10,:) + pln_coords(3,:)
		for( int it = 0; it < 3; it++ )				// Луна
			PLCOORD[9*3+it] += PLCOORD[2*3+it];

		for( int it = 0; it < 3; it++ )
			PLCOORD[10*3+it] = 0;		// Солнце
		
	}
	//==============================================================================//
	// гравитация от планет
	//==============================================================================//
	void InfluenceForce::planets_grav( double *x, double *f_gr )
	{
		double BC[3];
		double BO[3];
		double BC_r;
		double BO_r;

		f_gr[0] = 0;
		f_gr[1] = 0;
		f_gr[2] = 0;

		for( int i = 0; i <= 10; i++ ) 
		{
			if( pln_flag(i) )
			{
				// body - center
				pln_coords( i, BC );	

				//(BC,BC)
				BC_r = BC[0]*BC[0] + BC[1]*BC[1] + BC[2]*BC[2];	

				//BC_r.LT.1.d-12
				if ( BC_r < 1.0e-12 )	
				{
					//printf("BC_r < 1.0e-12 \n");
					// body - object
					BO[0] = BC[0]-x[0]; 
					BO[1] = BC[1]-x[1];
					BO[2] = BC[2]-x[2];
					// 1/(BO*BO)
					BO_r = 1.0/( BO[0]*BO[0] + BO[1]*BO[1] + BO[2]*BO[2] );
					BO_r = BO_r*sqrt(BO_r);

					f_gr[0] = f_gr[0] + mu_plan(i)*BO[0]*BO_r;
					f_gr[1] = f_gr[1] + mu_plan(i)*BO[1]*BO_r;
					f_gr[2] = f_gr[2] + mu_plan(i)*BO[2]*BO_r;
				}
				else
				{
					// body - object
					BO[0] = BC[0]-x[0]; 
					BO[1] = BC[1]-x[1];
					BO[2] = BC[2]-x[2];

					BC_r = 1.0/BC_r;
					BC_r = BC_r*sqrt(BC_r);

					// 1/(BO*BO)
					BO_r = 1.0/( BO[0]*BO[0] + BO[1]*BO[1] + BO[2]*BO[2] );
					BO_r = BO_r*sqrt(BO_r);

					f_gr[0] = f_gr[0] + mu_plan(i)*( BO[0]*BO_r - BC[0]*BC_r );
					f_gr[1] = f_gr[1] + mu_plan(i)*( BO[1]*BO_r - BC[1]*BC_r );
					f_gr[2] = f_gr[2] + mu_plan(i)*( BO[2]*BO_r - BC[2]*BC_r );
				}
			}
		}
	};
	//==============================================================================//
};