//==============================================================================//
// Andrianov N.G.
// opbit predict module
// Integration motion
// RK method
//==============================================================================//

#include <stdio.h>
#include <time.h>
#include "PredictOrbitSat.h" 

// name
namespace Orbit
{
	//==============================================================================//
	// функция в правой части
	// можем выбирать разные варианты
	//==============================================================================//
	void PredictOrbitSat::RightFxyzv( double t, double *xt, double *fx )
	{
		// учет только гравитации, гармоник, давления солнца (sun off)
		// похожая реализация есть для GPU
		//IF->rh_fun_grav( t, xt, fx );

		// Учет всех воздействий
		IF->rh_fun( t, xt, fx );

		// куплеровское движение
		//IF->rh_fun_kepler( t, xt, fx );
	};

	//==============================================================================//
	// вычисление начальных точек методом Рунге-Кутты
	//==============================================================================//
	void PredictOrbitSat::GetStartOrbitPoints( OrbitPoints & OP, double h )
	{
		//printf("Get Start Position by RK\n");
		int kkr[12] = { 3,4,5,6,1,2,4,3,2,1,6,5 };
		double g[6] = {
			0.069431844202973712388026755553595247452137,
			0.330009478207571867598667120448377656399712,
			0.669990521792428132401332879551622343600287,
			0.930568155797026287611973244446404752547862,
			1.069431844202973712388026755553595247452137,
			1.330009478207571867598667120448377656399712 
		};

		int NY = OP.NY;
		double *X0 = OP.X0;
		double TP0 = OP.TP0;

		//===============================================================//
		//if( kp == 1 || kp == 2 )
		// MAIN cycle of starting procedure
		// TR,HR - current time & current step
		// P  = F1+2*F2+2*F3+F4
		// fun(double t, double *xt, double *fx )
		double TPN;
		double XN[6];	// вектор состояния
		double FN1[6];	// значения функции
		double FN2[6];	// значения функции

		for (int j = 1; j <= 6; ++j)
		{
			int kk = kkr[j - 1];				// nodes following: 3,4,5,6,1,2 
			double hr = h * (1.0 - g[kk - 1]);	// step
			kk = kkr[j + 5];
			double tr = TP0 + hr;				// next time step


			RightFxyzv( TP0, X0, FN1 );
			//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

			// F2=F(T0+1/2H,X0+1/2HF1) 
			for( int it = 0; it < NY; it++ )
				XN[it] = hr * FN1[it]*0.5 + X0[it];
			TPN = TP0 + hr*0.5;
			RightFxyzv( TPN, XN, FN2 );
			//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

			// F3=F(T0+1/2H,X0+1/2HF2)
			for( int it = 0; it < NY; it++ )
			{
				XN[it] = hr * FN2[it]*0.5 + X0[it];
				FN1[it] += FN2[it]*2.0;
			}
			TPN = TP0 + hr*0.5;
			RightFxyzv( TPN, XN, FN2 );
			//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

			// F4=F(T0+H,X0+HF3)
			for( int it = 0; it < NY; it++ )
			{
				XN[it] = hr * FN2[it] + X0[it];
				FN1[it] += FN2[it]*2.0;
			}
			TPN = TP0 + hr;
			RightFxyzv( TPN, XN, FN2 );
			//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

			// X1=X0+1/6H(F1+2F2+2F3+F4)
			for( int it = 0; it < NY; it++ )
				XN[it] = X0[it] +  hr / 6.0 *( FN1[it] + FN2[it] );

			// вычисление точек
			// 4,3,2,1,6,5 []
			if( j == 1)	RightFxyzv(tr, XN, OP.FP4 );
			if( j == 2)	RightFxyzv(tr, XN, OP.FP3 );
			if( j == 3)	RightFxyzv(tr, XN, OP.FP2 );
			if( j == 4)	RightFxyzv(tr, XN, OP.FP1 );
			if( j == 5)	RightFxyzv(tr, XN, OP.FP6 );
			if( j == 6)	RightFxyzv(tr, XN, OP.FP5 );
		}
		//===============================================================//
	};
	//=====================================================================//
	// функция численного интегрирования
	//=====================================================================//
	void PredictOrbitSat::OrbitIntegration( OrbitPoints & OP, double t, CurrentIntegrateParam &Iparam, double *OrbitPointsArray )
	{
		//double e = 1.0e-8;
		int la[5] = { 0,24,48,72,96 };
		int la1[5] = { 0,8,16,24,32 };
		int la2[5] = { 0,24,24,24,48 };
		double alim[2] = { 0.0403536069, 1.0 };

		int isinv[5] = { 5,2,2,2,1 };
		int isge[5] = { 4,5,5,5,3 };
		int isle[5] = { 4,1,1,1,3 };
		int iseq[5] = { 4,2,2,2,3 };

		// получение констант
		const double *a = IC_a;
		const double *a1 = IC_a1;
		const double *w = IC_w;
		const double *g = IC_g;
		const double *c__ = IC_c;

		// начальные точки для старта интегратора
		double TP0 = OP.TP0;
		double TPN = OP.TPN;
		double *X0 = OP.X0;
		double *FP0 = OP.FP0;
		double *FP1 = OP.FP1;
		double *FP2 = OP.FP2;
		double *FP3 = OP.FP3;
		double *FP4 = OP.FP4;
		double *FP5 = OP.FP5;
		double *FP6 = OP.FP6;
		double *FP7 = OP.FP7;
		double Rst = t;

		// new t_steppt
		double t_steppt = TP0; // время старта записи точек

		//===============================================================//
		//double hh = h;
		//int NY = 6;	
		//int NP = 0;
		//int IP9 = 2;

		//int l = 0;
		//int l1 = 0;
		//int l2 = 0;
		//int ltek = 0;

		//int kkbeg = 1;
		//double delt, rr1, rr2, rotn, dd1;
		//bool InvertStart = false;

		double hh = Iparam.hh;
		
		int NY = Iparam.NY;	
		int NP = Iparam.NP;
		int IP9 = Iparam.IP9;

		int l = Iparam.l;
		int l1 = Iparam.l1;
		int l2 = Iparam.l2;
		int ltek = Iparam.ltek;

		int kkbeg = Iparam.kkbeg;
		double delt = Iparam.delt;
		double rr1 = Iparam.rr1;
		double rr2 = Iparam.rr2;
		double rotn = Iparam.rotn;
		double dd1 = Iparam.dd1;
		bool InvertStart = Iparam.InvertStart;

		//------------------------------------------------//
		// запись начального положения
		if( OrbitPointsArray != NULL )
		{
			// позиция для записи плюс 1, в самое начало
			int itwrite = 1;

			// запись значений
			double w_t = TP0;
			double w_x = X0[0];
			double w_y = X0[1];
			double w_z = X0[2];

			OrbitPointsArray[ itwrite ] = w_t;
			OrbitPointsArray[ itwrite + 1] = w_x;
			OrbitPointsArray[ itwrite + 2] = w_y;
			OrbitPointsArray[ itwrite + 3] = w_z;
		}
		//------------------------------------------------//

		// изменение направления интегрирования
		if (hh * (t - TP0) < 0.0)
		{ 
			// direct change.
			printf("direct change\n");
			for( int it = 0; it < NY; it++ )
			{
				double tmp = FP1[it];
				FP1[it] = FP6[it];
				FP6[it] = tmp;

				tmp = FP2[it];
				FP2[it] = FP5[it];
				FP5[it] = tmp;

				tmp = FP3[it];
				FP3[it] = FP4[it];
				FP4[it] = tmp;
			}

			// step change
			hh = -hh;
			if (IP9 == 1) {	hh /= 0.7;	}
			if (IP9 == 5) {	hh *= 0.7;	}

			// case nuber cange 
			IP9 = isinv[IP9 - 1];

			//  2 steps of extrapolation
			kkbeg = 3;

			//goto L1400;
			InvertStart = true;
		}
		//===============================================================//

		int isstep = 0;
		double minstep = 1000;
		// main cyrcle
		while( 1 )
		{
			if( InvertStart == false )
			{
				//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
				// Add new 
				//while(hh * ( t_steppt - TP0 - hh) <= 0.0 )
				//{
				//	// A8: interpolatin at the destination time 
				//	double delts = t_steppt - TP0;
				//	double Xres[6];
				//	for ( int ii = 0; ii < NY; ii++ ) 
				//	{
				//		rr1 = 0.0;
				//		for (int j = 1; j <= 13; j += 4) 
				//		{
				//			rr1 = delts / hh * (c__[j - 1] * FP3[ii] + c__[j] * FP4[ii] + c__[j + 1] * FP5[ii] + c__[j + 2] * FP6[ii] + rr1);
				//			Xres[ii] = FP0[ii] + rr1 * hh;
				//		}
				//	}

				//	FILE *fre = fopen( "pt.log", "at" );
				//	for ( int ii = 0; ii < NY; ii++ ) 
				//		fprintf( fre, "%f\t", Xres[ii] );
				//	fprintf( fre, "%f\n", t_steppt );
				//	fclose( fre );
				//	// next pt
				//	t_steppt += 0.001;
				//}
				//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//

				//===============================================================//
				//  S4: time overflow 
				if(hh * (t - TP0 - hh) <= 0.0 )
					break; 

				if( isstep >= 9990 )
				{
					printf( "bad satt\n" );
					break;
				}
				// if overflow absent => step of implicit RK-method 
				// A5: step of implicit RK-method + error estimation
				TP0 +=hh;

				// DELT = maximum error
				// RR1,RR2 - working values
				// ROTH - error for i-th equation
				delt = 0.0;
				for (int ii = 0; ii < NY; ii++ )
				{
					rr2 = FP3[ii]*w[72] + FP4[ii]*w[73] + FP5[ii]*w[74] + FP6[ii]*w[75];

					//step counter starting setting NP = 0;
					if (NP == 0 || kkbeg == 3 ){	
						FP0[ii] +=  hh * rr2;
						continue;
					}

					rr1 = FP7[ii] + FP2[ii]*a1[l1 + 4] + FP3[ii]*a1[l1 + 5] + FP4[ii]*a1[l1 + 6] + FP5[ii]*a1[l1 + 7];
					dd1 = hh * (rr1 - rr2)/Tolerance[ii];
					rotn = fabs(dd1);
					if (rotn >= delt){	delt = rotn;}

					FP0[ii] +=  hh * rr2;
				}
				NP++;
				//===============================================================//

				//===============================================================//
				//STEP: increase, decrease or do not change? 
				if ( NP == 1 || kkbeg == 3 ){ 
					IP9 = iseq[IP9 - 1];
					kkbeg = 1;
				}
				else if(alim[0] >= delt){	
					// A6:step increase 
					if (IP9 == 1) {	hh = hh; }
					if (IP9 == 2) {	hh /= 0.7;	}
					if (IP9 == 3) {	hh /= 0.7;	}
					if (IP9 == 4) {	hh /= 0.7;	}
					if (IP9 == 5) {	hh = hh;	}
					IP9 = isge[IP9 - 1];
				}
				else if (alim[1] <= delt) {	
					// A7:step decrease 
					if (IP9 == 1) {	hh = hh; }
					if (IP9 == 2) {	hh *= 0.7;	}
					if (IP9 == 3) {	hh *= 0.7;	}
					if (IP9 == 4) {	hh *= 0.7;	}
					if (IP9 == 5) {	hh = hh;	}
					IP9 = isle[IP9 - 1];
				}
				else
				{
					IP9 = iseq[IP9 - 1];
					kkbeg = 1;
				}
				//===============================================================//
			}
			InvertStart = false;
			//===============================================================//
			// A4: extrapolation 2|4 points + interpolation 
			l = la[IP9 - 1];
			l1 = la1[IP9 - 1];
			l2 = la2[IP9 - 1];
			//===============================================================//

			//===============================================================//
			for(int ii = 0; ii < NY; ii++ ) 
				FP7[ii] = 0.0;
			
			//  Main cycle of extrapolation 
			// FP0 + hh( a1*FP1 + a2*FP2 + a3*FP3 + a4*FP4 + a5*FP5  + a6*FP6 )
			for (int kk = kkbeg; kk <= 4; ++kk) 
			{
				ltek = l + kk * 6 - 6;
				// время в первой точке FP1
				TPN = TP0 + g[kk - 1] * hh;

				for(int ii = 0; ii < NY; ii++ ) 
				{
					X0[ii] = FP0[ii] + hh * (a[ltek] * FP1[ii] + a[ltek + 1] * FP2[ii] + a[ltek + 2] * FP3[ii] + a[ltek + 3] * FP4[ii] + a[ltek + 4]* FP5[ii] + a[ltek + 5] *FP6[ii]);
					// accumulate values for error estimation
					FP7[ii] += FP2[ii] * a1[l1 + kk - 1];
				}
				RightFxyzv( TPN, X0, FP1 );

				// offset point
				for( int it = 0; it < NY; it++ )
				{
					double tmp = FP1[it];
					FP1[it] = FP2[it];
					FP2[it] = FP3[it];
					FP3[it] = FP4[it];
					FP4[it] = FP5[it];
					FP5[it] = FP6[it];
					FP6[it] = tmp;
				}
			}
			//===============================================================//

			//===============================================================//
			// Main cycle of interpolation
			// X0 = FP0 + hh*( w1*FP1 + w2*FP2 + w3*FP3 + w4*FP4 + w5*FP5 + w6*FP6 )
			for (int kk = 1; kk <= 4; ++kk)
			{
				for ( int ii = 0; ii < NY; ii++ )
					X0[ii] =  FP0[ii] + hh * (w[l2]*FP1[ii] + w[l2 + 1]*FP2[ii] + w[l2 + 2]*FP3[ii] + w[l2 + 3]*FP4[ii] + w[l2 + 4]*FP5[ii] + w[l2 + 5]*FP6[ii] );

				// new point for FP3 .... FP6
				TPN = TP0 + g[kk - 1] * hh;
				if( kk == 1) RightFxyzv( TPN, X0, FP3 );
				if( kk == 2) RightFxyzv( TPN, X0, FP4 );
				if( kk == 3) RightFxyzv( TPN, X0, FP5 );
				if( kk == 4) RightFxyzv( TPN, X0, FP6 );
				l2 += 6;
			}
			//===============================================================//

			//printf("%f  ", hh );
			if( hh < minstep )
				minstep = hh;

			//------------------------------------------------//
			isstep++;
			if( OrbitPointsArray != NULL )
			{
				// позиция для записи начиная со второго места плюс 1
				int itwrite = isstep*4 + 1;

				// запись значений
				double w_t = TP0;
				double w_x = FP0[0];
				double w_y = FP0[1];
				double w_z = FP0[2];

				OrbitPointsArray[ itwrite ] = w_t;
				OrbitPointsArray[ itwrite + 1] = w_x;
				OrbitPointsArray[ itwrite + 2] = w_y;
				OrbitPointsArray[ itwrite + 3] = w_z;
			}
			//------------------------------------------------//

			//FILE *fre = fopen( "pt_st.log", "at" );
			//for ( int ii = 0; ii < NY; ii++ ) 
			//	fprintf( fre, "%f\t", FP0[ii] );
			//fprintf( fre, "%f\n", TP0 );
			//fclose( fre );
		}
		//===============================================================//
		// A8: interpolatin at the destination time 
		delt = t - TP0;
		for ( int ii = 0; ii < NY; ii++ ) 
		{
			rr1 = 0.0;
			for (int j = 1; j <= 13; j += 4) 
			{
				rr1 = delt / hh * (c__[j - 1] * FP3[ii] + c__[j] * FP4[ii] + c__[j + 1] * FP5[ii] + c__[j + 2] * FP6[ii] + rr1);
				X0[ii] = FP0[ii] + rr1 * hh;
			}
		}
		//===============================================================//
		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
		// Add new 
		//while(hh * ( t_steppt - t) <= 0.0 )
		//{
		//	// A8: interpolatin at the destination time 
		//	double delts = t_steppt - TP0;
		//	double Xres[6];
		//	for ( int ii = 0; ii < NY; ii++ ) 
		//	{
		//		rr1 = 0.0;
		//		for (int j = 1; j <= 13; j += 4) 
		//		{
		//			rr1 = delts / hh * (c__[j - 1] * FP3[ii] + c__[j] * FP4[ii] + c__[j + 1] * FP5[ii] + c__[j + 2] * FP6[ii] + rr1);
		//			Xres[ii] = FP0[ii] + rr1 * hh;
		//		}
		//	}

		//	FILE *fre = fopen( "pt.log", "at" );
		//	for ( int ii = 0; ii < NY; ii++ ) 
		//		fprintf( fre, "%f\t", Xres[ii] );
		//	fprintf( fre, "%f\n", t_steppt );
		//	fclose( fre );
		//	// next pt
		//	t_steppt += 0.001;
		//}
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//

		// возвращаем шаг интегрирования
		Iparam.hh = hh;
		Iparam.NY = NY;	
		Iparam.NP = NP;
		Iparam.IP9 = IP9;

		Iparam.l = l;
		Iparam.l1 = l1;
		Iparam.l2 = l2;
		Iparam.ltek = ltek;

		Iparam.kkbeg = kkbeg;
		Iparam.delt = delt;
		Iparam.rr1 = rr1;
		Iparam.rr2 = rr2;
		Iparam.rotn = rotn;
		Iparam.dd1 = dd1;
		Iparam.InvertStart = InvertStart;

		// начальное время, это текущее законченное
		OP.TP0 = TP0;

		//------------------------------------------------//
		// число точек
		if( OrbitPointsArray != NULL )
		{
			OrbitPointsArray[ 0 ] = isstep + 1;
		}
		//------------------------------------------------//

		printf( "\nN step = %d\nminstep = %f\n", isstep, minstep );

	}
	//##############################################################################//
	//							NEW FUNC FOR FIND DIST								//
	//==============================================================================//
	// save step by 1 sec
	//==============================================================================//
	void PredictOrbitSat::OrbitIntegration_step( OrbitPoints & OP, double t, CurrentIntegrateParam &Iparam, double *OrbitPointsArray )
	{
		//double e = 1.0e-8;
		int la[5] = { 0,24,48,72,96 };
		int la1[5] = { 0,8,16,24,32 };
		int la2[5] = { 0,24,24,24,48 };
		double alim[2] = { 0.0403536069, 1.0 };

		int isinv[5] = { 5,2,2,2,1 };
		int isge[5] = { 4,5,5,5,3 };
		int isle[5] = { 4,1,1,1,3 };
		int iseq[5] = { 4,2,2,2,3 };

		// получение констант
		const double *a = IC_a;
		const double *a1 = IC_a1;
		const double *w = IC_w;
		const double *g = IC_g;
		const double *c__ = IC_c;

		// начальные точки для старта интегратора
		double TP0 = OP.TP0;
		double TPN = OP.TPN;
		double *X0 = OP.X0;
		double *FP0 = OP.FP0;
		double *FP1 = OP.FP1;
		double *FP2 = OP.FP2;
		double *FP3 = OP.FP3;
		double *FP4 = OP.FP4;
		double *FP5 = OP.FP5;
		double *FP6 = OP.FP6;
		double *FP7 = OP.FP7;
		double Rst = t;

		// new t_steppt
		double t_steppt = TP0; // время старта записи точек
		int itwrites = 0;
		double PointTStep = STEPDTIME/1000.0;

		int isstep = 0;
		double minstep = 1000;
		//===============================================================//
		//double hh = h;
		//int NY = 6;	
		//int NP = 0;
		//int IP9 = 2;

		//int l = 0;
		//int l1 = 0;
		//int l2 = 0;
		//int ltek = 0;

		//int kkbeg = 1;
		//double delt, rr1, rr2, rotn, dd1;
		//bool InvertStart = false;

		double hh = Iparam.hh;
		int NY = Iparam.NY;	
		int NP = Iparam.NP;
		int IP9 = Iparam.IP9;

		int l = Iparam.l;
		int l1 = Iparam.l1;
		int l2 = Iparam.l2;
		int ltek = Iparam.ltek;

		int kkbeg = Iparam.kkbeg;
		double delt = Iparam.delt;
		double rr1 = Iparam.rr1;
		double rr2 = Iparam.rr2;
		double rotn = Iparam.rotn;
		double dd1 = Iparam.dd1;
		bool InvertStart = Iparam.InvertStart;

		// изменение направления интегрирования
		if (hh * (t - TP0) < 0.0)
		{ 
			// direct change.
			//printf("direct change\n");
			for( int it = 0; it < NY; it++ )
			{
				double tmp = FP1[it];
				FP1[it] = FP6[it];
				FP6[it] = tmp;

				tmp = FP2[it];
				FP2[it] = FP5[it];
				FP5[it] = tmp;

				tmp = FP3[it];
				FP3[it] = FP4[it];
				FP4[it] = tmp;
			}

			// step change
			hh = -hh;
			if (IP9 == 1) {	hh /= 0.7;	}
			if (IP9 == 5) {	hh *= 0.7;	}

			// case nuber cange 
			IP9 = isinv[IP9 - 1];

			//  2 steps of extrapolation
			kkbeg = 3;

			//goto L1400;
			InvertStart = true;
		}
		//===============================================================//
		// main cyrcle
		while( 1 )
		{
			if( InvertStart == false )
			{
				//===============================================================//
				//  S4: time overflow 
				if(hh * (t - TP0 - hh) <= 0.0 )
					break; 

				//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
				// Add new 
				if( OrbitPointsArray != NULL )
				{
					while(hh * ( t_steppt - TP0 - hh) <= 0.0 )
					{
						// A8: interpolatin at the destination time 
						double delts = t_steppt - TP0;
						double Xres[6];
						for ( int ii = 0; ii < NY; ii++ ) 
						{
							rr1 = 0.0;
							for (int j = 1; j <= 13; j += 4) 
							{
								rr1 = delts / hh * (c__[j - 1] * FP3[ii] + c__[j] * FP4[ii] + c__[j + 1] * FP5[ii] + c__[j + 2] * FP6[ii] + rr1);
								Xres[ii] = FP0[ii] + rr1 * hh;
							}
						}
						//FILE *fre = fopen( "pt.log", "at" );
						//for ( int ii = 0; ii < NY; ii++ ) 
						//	fprintf( fre, "%f\t", Xres[ii] );
						//fprintf( fre, "%f\n", t_steppt );
						//fclose( fre );

						int iw = itwrites*4 + 1;
						OrbitPointsArray[ iw  ] = t_steppt;
						OrbitPointsArray[ iw  + 1] = Xres[0];
						OrbitPointsArray[ iw  + 2] = Xres[1];
						OrbitPointsArray[ iw  + 3] = Xres[2];
						itwrites++;
					
						// next pt
						t_steppt +=  PointTStep;
					}
				}
				//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//

				// if overflow absent => step of implicit RK-method 
				// A5: step of implicit RK-method + error estimation
				TP0 +=hh;

				// DELT = maximum error
				// RR1,RR2 - working values
				// ROTH - error for i-th equation
				delt = 0.0;
				for (int ii = 0; ii < NY; ii++ )
				{
					rr2 = FP3[ii]*w[72] + FP4[ii]*w[73] + FP5[ii]*w[74] + FP6[ii]*w[75];

					//step counter starting setting NP = 0;
					if (NP == 0 || kkbeg == 3 ){	
						FP0[ii] +=  hh * rr2;
						continue;
					}

					rr1 = FP7[ii] + FP2[ii]*a1[l1 + 4] + FP3[ii]*a1[l1 + 5] + FP4[ii]*a1[l1 + 6] + FP5[ii]*a1[l1 + 7];
					dd1 = hh * (rr1 - rr2)/Tolerance[ii];
					rotn = fabs(dd1);
					if (rotn >= delt){	delt = rotn;}

					FP0[ii] +=  hh * rr2;
				}
				NP++;
				//===============================================================//

				//===============================================================//
				//STEP: increase, decrease or do not change? 
				if ( NP == 1 || kkbeg == 3 ){ 
					IP9 = iseq[IP9 - 1];
					kkbeg = 1;
				}
				else if(alim[0] >= delt){	
					// A6:step increase 
					if (IP9 == 1) {	hh = hh; }
					if (IP9 == 2) {	hh /= 0.7;	}
					if (IP9 == 3) {	hh /= 0.7;	}
					if (IP9 == 4) {	hh /= 0.7;	}
					if (IP9 == 5) {	hh = hh;	}
					IP9 = isge[IP9 - 1];
				}
				else if (alim[1] <= delt) {	
					// A7:step decrease 
					if (IP9 == 1) {	hh = hh; }
					if (IP9 == 2) {	hh *= 0.7;	}
					if (IP9 == 3) {	hh *= 0.7;	}
					if (IP9 == 4) {	hh *= 0.7;	}
					if (IP9 == 5) {	hh = hh;	}
					IP9 = isle[IP9 - 1];
				}
				else
				{
					IP9 = iseq[IP9 - 1];
					kkbeg = 1;
				}
				//===============================================================//
			}
			InvertStart = false;
			//===============================================================//
			// A4: extrapolation 2|4 points + interpolation 
			l = la[IP9 - 1];
			l1 = la1[IP9 - 1];
			l2 = la2[IP9 - 1];
			//===============================================================//

			//===============================================================//
			for(int ii = 0; ii < NY; ii++ ) 
				FP7[ii] = 0.0;

			//  Main cycle of extrapolation 
			// FP0 + hh( a1*FP1 + a2*FP2 + a3*FP3 + a4*FP4 + a5*FP5  + a6*FP6 )
			for (int kk = kkbeg; kk <= 4; ++kk) 
			{
				ltek = l + kk * 6 - 6;
				// время в первой точке FP1
				TPN = TP0 + g[kk - 1] * hh;

				for(int ii = 0; ii < NY; ii++ ) 
				{
					X0[ii] = FP0[ii] + hh * (a[ltek] * FP1[ii] + a[ltek + 1] * FP2[ii] + a[ltek + 2] * FP3[ii] + a[ltek + 3] * FP4[ii] + a[ltek + 4]* FP5[ii] + a[ltek + 5] *FP6[ii]);
					// accumulate values for error estimation
					FP7[ii] += FP2[ii] * a1[l1 + kk - 1];
				}
				RightFxyzv( TPN, X0, FP1 );

				// offset point
				for( int it = 0; it < NY; it++ )
				{
					double tmp = FP1[it];
					FP1[it] = FP2[it];
					FP2[it] = FP3[it];
					FP3[it] = FP4[it];
					FP4[it] = FP5[it];
					FP5[it] = FP6[it];
					FP6[it] = tmp;
				}
			}
			//===============================================================//

			//===============================================================//
			// Main cycle of interpolation
			// X0 = FP0 + hh*( w1*FP1 + w2*FP2 + w3*FP3 + w4*FP4 + w5*FP5 + w6*FP6 )
			for (int kk = 1; kk <= 4; ++kk)
			{
				for ( int ii = 0; ii < NY; ii++ )
					X0[ii] =  FP0[ii] + hh * (w[l2]*FP1[ii] + w[l2 + 1]*FP2[ii] + w[l2 + 2]*FP3[ii] + w[l2 + 3]*FP4[ii] + w[l2 + 4]*FP5[ii] + w[l2 + 5]*FP6[ii] );

				// new point for FP3 .... FP6
				TPN = TP0 + g[kk - 1] * hh;
				if( kk == 1) RightFxyzv( TPN, X0, FP3 );
				if( kk == 2) RightFxyzv( TPN, X0, FP4 );
				if( kk == 3) RightFxyzv( TPN, X0, FP5 );
				if( kk == 4) RightFxyzv( TPN, X0, FP6 );
				l2 += 6;
			}
			//===============================================================//

			//printf("%f  ", hh );
			if( hh < minstep )
				minstep = hh;

			isstep++;
		}
		//===============================================================//
		// A8: interpolatin at the destination time 
		delt = t - TP0;
		for ( int ii = 0; ii < NY; ii++ ) 
		{
			rr1 = 0.0;
			for (int j = 1; j <= 13; j += 4) 
			{
				rr1 = delt / hh * (c__[j - 1] * FP3[ii] + c__[j] * FP4[ii] + c__[j + 1] * FP5[ii] + c__[j + 2] * FP6[ii] + rr1);
				X0[ii] = FP0[ii] + rr1 * hh;
			}
		}
		//===============================================================//
		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
		// Add new 
		if( OrbitPointsArray != NULL )
		{
			while( hh*( t_steppt - t) <= 0.0 )
			{
				// A8: interpolatin at the destination time 
				double delts = t_steppt - TP0;
				double Xres[6];
				for ( int ii = 0; ii < NY; ii++ ) 
				{
					rr1 = 0.0;
					for (int j = 1; j <= 13; j += 4) 
					{
						rr1 = delts / hh * (c__[j - 1] * FP3[ii] + c__[j] * FP4[ii] + c__[j + 1] * FP5[ii] + c__[j + 2] * FP6[ii] + rr1);
						Xres[ii] = FP0[ii] + rr1 * hh;
					}
				}
				//FILE *fre = fopen( "pt.log", "at" );
				//for ( int ii = 0; ii < NY; ii++ ) 
				//	fprintf( fre, "%f\t", Xres[ii] );
				//fprintf( fre, "%f\n", t_steppt );
				//fclose( fre );

				int iw = itwrites*4 + 1;
				OrbitPointsArray[ iw ] = t_steppt;
				OrbitPointsArray[ iw + 1] = Xres[0];
				OrbitPointsArray[ iw + 2] = Xres[1];
				OrbitPointsArray[ iw + 3] = Xres[2];
				itwrites++;

				// next pt
				t_steppt += PointTStep;
			}
		}
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//

		// возвращаем шаг интегрирования
		Iparam.hh = hh;
		Iparam.NY = NY;	
		Iparam.NP = NP;
		Iparam.IP9 = IP9;

		Iparam.l = l;
		Iparam.l1 = l1;
		Iparam.l2 = l2;
		Iparam.ltek = ltek;

		Iparam.kkbeg = kkbeg;
		Iparam.delt = delt;
		Iparam.rr1 = rr1;
		Iparam.rr2 = rr2;
		Iparam.rotn = rotn;
		Iparam.dd1 = dd1;
		Iparam.InvertStart = InvertStart;

		// начальное время, это текущее законченное
		OP.TP0 = TP0;

		//------------------------------------------------//
		// число точек
		if( OrbitPointsArray != NULL )
		{
			OrbitPointsArray[ 0 ] = itwrites;
			printf("Arr size = %d\n", itwrites );
		}
		//------------------------------------------------//
		//printf( "\nPredict: N step = %d minstep = %f\n", isstep, minstep );
	}
	//==============================================================================//
};