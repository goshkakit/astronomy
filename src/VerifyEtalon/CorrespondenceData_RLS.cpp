#include "CorrespondenceData.h"

extern "C" {
#include "novas.h"
#include "eph_manager.h" /* remove this line for use with solsys version 2 */
}
#include <conio.h>

// структура с временем измерени€ и направлением приема
// дл€ всех поментов записи необходимо будет получить расположение объекта
// координаты и углы
struct ObzorRLS
{
	std::vector< double > tsev_tnc;
	std::vector< double > az0;
	std::vector< double > pla0;
};
//==============================================================================//
// определение точности измерений станции
//==============================================================================//
void CorrespondenceData::RunCorrespondenceDataRLS(  )
{
	printf( "-------------------------------------------------------\n" );

	//-------------------------------------------------------//
	// положение станции  алининград
	double lg = 20.1812187944;	// долгота grad
	double lt = 54.8571174472;	// широта grad
	double H = 98.302;			// высота m

	// ѕоложение телескопа  исловодск
	//double lg = 42.661662;	// долгота grad
	//double lt = 43.72202;	// широта grad
	//double H = 2070.0;			// высота m

	// положение телескопа в москве
	//double lg = 37.54794;
	//double lt = 55.80364;
	//double H = 190.0;

	double Tpos_mos[3] = { 2848.74429, 2189.68052, 5252.33756 };
	//Xs = 2848.74429;
	//Ys = 2189.68052;
	//Zs = 5252.33756;

	double Tpos_kisl[3] = { 3395.61326, 3128.29531, 4388.71898 };
	printf( "Station Position Kisl orig\n" );
	printf( "%f\t %f\t %f\n", Tpos_kisl[0], Tpos_kisl[1], Tpos_kisl[2] );
	printf( "Station Position Mos orig\n" );
	printf( "%f\t %f\t %f\n", Tpos_mos[0], Tpos_mos[1], Tpos_mos[2] );

	double R_mean = 6363.0;		// средний радиус земли km
	printf( "Station position\n" );
	printf( "Lt = %f\n", lt );
	printf( "Lg = %f\n", lg );
	printf( "H = %f\n", H );

	// направление обзора станции
	double Az_st = 240.0205;		// направление - азимут grad

	double Hkm = H/1000.0;
	//double pi = 3.1415926535;
	double torad = pi/180.0;
	double lgr = lg*torad;
	double ltr = lt*torad;

	// перевод в декартову систему координат
	double Xs = ( R_mean + Hkm )*cos( ltr )*cos( lgr );
	double Ys = ( R_mean + Hkm )*cos( ltr )*sin( lgr );
	double Zs = ( R_mean + Hkm )*sin( ltr );

	// координаты станции
	double Tpos[3];
	Tpos[0] = Xs;
	Tpos[1] = Ys;
	Tpos[2] = Zs;
	printf( "Station Position\n" );
	printf( "%f\t %f\t %f\n", Tpos[0], Tpos[1], Tpos[2] );
	//-------------------------------------------------------//

	//-------------------------------------------------------//
	//  аталог
	//char *TLEFile = "TLE20140714.txt";
	char *TLEFile = "catalog-20140711-0000.txt";
	//char *TLEFile = "TLE20131020.txt";

	// «агрузка.  аталог NORAD
	tleload.LoadData( TLEFile, 2 );
	//int id = 39491; //4288;
	//int id = 33105;

	for( int is = 0; is < tleload.NORADList.size(); is++ )
	{
		printf( "-------------------------------------------------------\n" );
		//orbitN = tleload.GetFromID( id );
		orbitN = tleload.NORADList[ is ];
		int id = is;

		printf( "Get Sat = %s\n", orbitN->SatId().c_str() );
		printf( "Get Sat Name = %s\n", orbitN->SatName().c_str() );
		printf( "AP = %f PR = %f E = %f\n", orbitN->Apogee(), orbitN->Perigee(), orbitN->Eccentricity() );
		//-------------------------------------------------------//

		// матрица перехода
		// отстование по времени
		//double Toffset = 0;
		// момент обзора из записей
		//ObzorRLS orls;
		//for( int it = 0; it < orls.tsev_tnc.size(); it++ )
		//{
		//	double tsev = orls.tsev_tnc[it];
		//}

		//  онвертаци€ в московский формат времени
		// ¬рем€ по ћоскве
		// Result JD = 2456853.0022222 = 20140714 160312 - 4UTC
		double Tmdb_d = 20140714.0;	// YYYYMMDD
		double Tmdb_t = 161125.66889953 + 16280E-6; //161116;//150312.0;	// HHMMSS.SS
		double Tutc_d = 20140714.0;	// YYYYMMDD
		double Tutc_t = 131125.66889953 + 16280E-6; //131116;//120312.0;	// HHMMSS.SS

		//201013 013841647
		//double Tmdb_d = 20131020.0;	// YYYYMMDD
		//double Tmdb_t = 043841.647;	// HHMMSS.SS
		//double Tutc_d = 20131020.0;	// YYYYMMDD
		//double Tutc_t = 013841.647;	// HHMMSS.SS

		double jd_orig = 2.45658556853758E+0006;
		double Ra_orig = 6.14682169938829E+0000;
		double Dec_orig = 1.00405374356832E+0000;
		//-------------------------------------------------------//
		// установка времени
		double int1, ajd1, delt1;
		double date1 = Tmdb_d;	
		double time1 = Tmdb_t;		
		IForce->set_time(date1, time1, &ajd1, &delt1, &int1 );

		// матрица перевода в земную систему
		double Arot[9];
		IForce->iers_update_matrix( int1, Arot, ajd1, delt1 );
		// матрица перехода из земной в нормальную систему
		double invArot[9];
		IForce->transpose( Arot, invArot ); 

		// перевод в ICRF координат телескопа
		double Telicrf[3];
		IForce->matVecMul(  invArot, Tpos, Telicrf );

		printf( "Station Position icrf\n" );
		printf( "%f\t %f\t %f\n", Telicrf[0], Telicrf[1], Telicrf[2] );

		double jd = IForce->dt_ajd( Tutc_d ) - 0.5;
		double secday_m = 12.0*3600.0 + 3.0*60.0 + 12.0;
	
		double tm_i = Tutc_t;
		int h = tm_i/10000;
		int m = (tm_i - h*10000)/100;
		double s = (double)(tm_i - h*10000.0 - m*100.0);
		double secday	= ((double)(h))*3600.0 + ((double)(m))*60.0 + s;
		jd = jd + secday/86400.0;

		printf( "JD Evnt = %f secday = %f %f\n" , jd, secday_m, secday );
		printf( "Delta jd orig = %f\n", jd_orig - jd );
		//-------------------------------------------------------//

		//-------------------------------------------------------//
		// прогноз орбиты на момент измерени€
		double Sicrf[3];
		Sicrf[0] = 0;
		Sicrf[1] = 0;
		Sicrf[2] = 0;

		double jdmeg = jd;

		// прогнозируем “Ћ≈ на момент измерени€
		// прогноз по SGP4
		double Tpredict = jdmeg - orbitN->Epoch().Date();
		Tpredict = Tpredict*24.0*60.0;
		printf( "Tpredict = %f min %f hour\n", Tpredict, Tpredict/60.0 );

		double Ptle[3];
		try
		{
			// прогноз положени€ в километрах
			cEciTime eci1 = orbitN->GetPosition( Tpredict );
			Ptle[0] = eci1.Position().m_x;
			Ptle[1] = eci1.Position().m_y;
			Ptle[2] = eci1.Position().m_z;
		}
		catch( ... )
		{
			printf("NORAD ERROR Tpredict = %f min, it = %d\n", Tpredict, id );
		}
		printf( "Object Position TEME\n" );
		printf( "%f\t %f\t %f\n", Ptle[0], Ptle[1], Ptle[2] );

		// нельз€ присваивать, в разных системах координат
		//Sicrf[0] = Ptle[0];
		//Sicrf[1] = Ptle[1];
		//Sicrf[2] = Ptle[2];

		// ConvertTEMEtoICRF
		// переод спрогнозированного значени€ в систему координат ICRF
		// установка времени
		//double date1 = calcTime.gdata;
		//double time1 = calcTime.gtime; 
		//double int1, ajd1, delt1;
		//IForce->set_time(date1, time1, &ajd1, &delt1, &int1 );

		// матрица перехода из ICRF в TEME
		double A_Teme[9];
		IForce->GetTemeMatrix( int1, A_Teme, ajd1, delt1 );
		// обратна€ матрица перехода из TEME в ICRF
		double invA_Teme[9];
		IForce->transpose( A_Teme, invA_Teme ); 

		// перевод вектора состо€ни€ в ICRF 
		IForce->matVecMul( invA_Teme, Ptle, Sicrf );

		printf( "Object Position Sicrf\n" );
		printf( "%f\t %f\t %f\n", Sicrf[0], Sicrf[1], Sicrf[2] );
		//-------------------------------------------------------//

		//-------------------------------------------------------//
		// вычисление просчитанных измерений
		double Ra;
		double Dec;
		// вектор направлени€
		double pos[3];
		pos[0] = Sicrf[0] - Telicrf[0];
		pos[1] = Sicrf[1] - Telicrf[1];
		pos[2] = Sicrf[2] - Telicrf[2];

		// ”гол между вектором вверх и направлением на объект
		double cosa = pos[0]*Telicrf[0] + pos[1]*Telicrf[1] + pos[2]*Telicrf[2];
		double abs1 =  sqrt( pos[0]*pos[0] +  pos[1]*pos[1] +  pos[2]*pos[2] );
		double abs2 =  sqrt( Telicrf[0]*Telicrf[0] +  Telicrf[1]*Telicrf[1] +  Telicrf[2]*Telicrf[2] );
		cosa = cosa/abs1/abs2;
		double az1 = acos( cosa );
		az1 = az1/torad;
		printf( "az1 = %f\n", az1 );

		ConvertXYZtoRADEC( Sicrf, Telicrf, &Ra, &Dec );
		printf( "Ra = %f\t Dec = %f\n", Ra, Dec );

		double delta_ra = Ra_orig - Ra;
		double delta_dec = Dec_orig - Dec;
		printf( "Detla angle = %f %f rad\n", delta_ra, delta_dec );
		printf( "Detla angle = %f %f degree\n", delta_ra/torad, delta_dec/torad );

		// дальность до объекта
		double satDist = sqrt( pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2] );
		double rDist = sqrt(  Sicrf[0]*Sicrf[0] + Sicrf[1]*Sicrf[1] + Sicrf[2]*Sicrf[2] );
		printf("dist = %f\t %f\n", satDist, rDist );

		vector2radec( pos, &Ra, &Dec );
		double Ra_h = Ra;
		double Dec_deg = Dec;

		Ra = Ra/12.0*pi;
		Dec = Dec/180.0*pi;
		printf( "Ra = %f\t Dec = %f\n", Ra, Dec );
		printf( "Ra_o = %f\t Dec_o = %f\n", Ra_orig, Dec_orig );
		//-------------------------------------------------------//

		//-------------------------------------------------------//
		// ѕеревод в азимут и угол места
		double xp = 0.082004;
		double yp = 0.298188;
		double ut1_utc =  0.2857338;

		double leap_secs = 35.0;
		double jd_ut1 = jdmeg + ut1_utc/86400.0;
		double delta_t = 32.184 + leap_secs - ut1_utc;
		short int accuracy = 0;

		on_surface location;
		location.height = H;
		location.latitude = lt;
		location.longitude = lg;

		short int ref_option = 0;

		double zd;
		double az;
		double rar;
		double decr;
		equ2hor ( jd_ut1, delta_t, accuracy, xp, yp, &location, Ra_h, Dec_deg, ref_option, &zd, &az, &rar, &decr );
		printf("\n View position satellite\n" );
		printf( "zd = %f\naz = %f\nrar = %f\ndecr = %f\n", zd, az, rar, decr );


		double zenangle = GetZenitAngle( Ra, Dec, Telicrf );
		zenangle = zenangle/pi*180.0;
		printf( "zenangle = %f\n", zenangle );
		//-------------------------------------------------------//

		printf( "Delta Az = %f\n", az-Az_st );

		// Good satt
		bool goodSatt = false;

		double res_z = 90.0 - zenangle;
		double res_az = az - Az_st;

		double view_z = 3.0;
		double view_az = -7.0;
		double delta_z = 2.0;
		double delta_az = 2.0;

		if( abs( res_z - view_z ) < delta_z && abs( res_az - view_az ) < delta_az )
			goodSatt = true;

		if( goodSatt )
		{
			//getch();
			FILE *fre = fopen( "listSatt.txt", "at" );
			fprintf( fre, "%d\t %s\t %f\t %f\t %f\t %f\t %f\n", id, orbitN->SatId().c_str(), zenangle, az, res_z, res_az, satDist );
			fclose( fre );
		}
		printf( "-------------------------------------------------------\n" );
	}

	return;
};
//==============================================================================//