#include "CorrespondenceData.h"

extern "C" {
#include "novas.h"
#include "eph_manager.h" /* remove this line for use with solsys version 2 */
}

//################################### NOVAS ####################################//
//==============================================================================//
// земныек координата в неподвижные координаты
//==============================================================================//
void CorrespondenceData::Ter2Cel( double jday, double jdmeg, double *vec1, double *vec2 )
{
	// поправки полюса
	//13 1 1 56293.00 P  0.082004 0.005865  0.298188 0.007111  P 0.2857338 0.0055637                 P   -75.505     .600    -8.527     .600                                                     
	//13 1 2 56294.00 P  0.080977 0.005927  0.298532 0.007205  P 0.2847340 0.0056502                 P   -75.505     .600    -8.480     .600                                                     
	//13 1 3 56295.00 P  0.079960 0.005988  0.298894 0.007299  P 0.2836315 0.0057362                 P   -75.519     .600    -8.293     .600	
	double xp = 0.082004;
	double yp = 0.298188;
	double ut1_utc =  0.2857338;

	// 13 3 7 56358.00 P  0.049385 0.009202  0.348086 0.012560  P 0.2133315 0.0105444 
	//double xp = 0.049385;
	//double yp = 0.348086;
	//double ut1_utc =  0.2133315;


	double leap_secs = 35.0;
	double jd_ut_high = jday;
	double jd_ut_low = jdmeg - jday + ut1_utc/86400.0;
	double delta_t = 32.184 + leap_secs - ut1_utc;
	short int method = 0;
	short int accuracy = 0;
	short int option = 0;
	short int res = ter2cel( jd_ut_high, jd_ut_low, delta_t, method, accuracy, option, xp, yp, vec1, vec2 );
	//printf("res = %d\n", res );

}
//==============================================================================//
// доработать
// перевод координат
//==============================================================================//
void CorrespondenceData::EquToHor( double jd, double *pos, double *zd, double *az, double *rar, double *decr )
{
	double ra, dec;
	vector2radec( pos, &ra, &dec );

	double xp = 0.082004;
	double yp = 0.298188;
	double ut1_utc =  0.2857338;

	double leap_secs = 35.0;
	double jd_ut1 = jd + ut1_utc/86400.0;
	double delta_t = 32.184 + leap_secs - ut1_utc;
	short int accuracy = 0;

	on_surface location;
	location.height = 190.0;
	location.latitude = 55.80364;
	location.longitude = 37.54794;

	short int ref_option = 0;

	equ2hor ( jd_ut1, delta_t, accuracy, xp, yp, &location, ra, dec, ref_option, zd, az, rar, decr );
}; 
//==============================================================================//
// функция идентификации измерений.
// поиск среди измерений по заданному номеру спутника и загруженному TLE 
// поиск заданного спутника в измерениях
//==============================================================================//
//double CorrespondenceData::TLEidentify( char *fname, double *Tel_ITRF, bool printrep )
//{
//	// загрузка файлов
//	optLoad.LoadData( fname, false );
//
//	// получеие списка измерений
//	OpticArray = optLoad.GetOpticList();
//
//	// проверяем весь список измерений
//	double distAngle = 0;
//	for( int it = 0; it < OpticArray.size(); it++ )
//	{
//		double jdmeg = OpticArray[it].jd;
//
//		// прогнозируем ТЛЕ на момент измерения м
//		// прогноз по SGP4
//		double Tpredict = jdmeg - orbitN->Epoch().Date();
//		Tpredict = Tpredict*24.0*60.0;
//		cEciTime eci1 = orbitN->GetPosition( Tpredict );
//
//		double Ptle[3];
//		// прогноз положения в километрах
//		Ptle[0] = eci1.Position().m_x;
//		Ptle[1] = eci1.Position().m_y;
//		Ptle[2] = eci1.Position().m_z;
//
//		// нужно перевести прогноз из TEME в ICRF
//		double PICRF[3];
//
//		TimePoint calcTime;
//		calcTime.jd = OpticArray[it].jd;
//		calcTime.gdata = OpticArray[it].dataMDB;
//		calcTime.gtime = OpticArray[it].timeMDB;
//
//		ConvertTEMEtoICRF( Ptle, PICRF, calcTime );
//
//		// перевести в RA и DEC
//		// положение телескопа на момент измерения в ICRF в км
//		double Tel_ICRF[3];
//		GetTelescopePosInICRF( Tel_ITRF, Tel_ICRF, calcTime );
//
//		// вычисление просчитанных измерений
//		double Ra;
//		double Dec;
//		ConvertXYZtoRADEC( PICRF, Tel_ICRF, &Ra, &Dec );
//
//		double pi = 3.1415926535;
//		double tog = 180.0/pi;
//		if( Ra < 0 )
//			Ra = 2.0*pi + Ra;
//
//		double optRa =  OpticArray[it].Ra;
//		double optDec =  OpticArray[it].Dec;
//
//		double dRa = Ra - optRa;
//		double dDec = Dec - optDec;
//
//		distAngle += abs( dRa ) + abs( dDec );
//
//		if( printrep )
//		{
//			printf( "data: %.2f %.2f\tOPT: %f\t%f\tTLE: %f\t%f\n", OpticArray[it].dataMDB, OpticArray[it].timeMDB, optRa*tog, optDec*tog, Ra*tog, Dec*tog );
//
//			printf( "%d:rad:\t dRa = %.7f\t  dDec = %.7f\n", it, dRa, dDec );
//
//			dRa = dRa/pi*180.0*3600.0;
//			dDec = dDec/pi*180.0*3600.0;
//
//			printf( "%d:sec:\t dRa = %.7f\t  dDec = %.7f\n-----\n", it, dRa, dDec );
//		}
//
//	}
//	distAngle = distAngle/OpticArray.size();
//
//	return distAngle;
//}
//==============================================================================//
