#include "CorrespondenceData.h"
#include "common/DataConverter.h"

extern "C" {
#include "novac/novas.h"
#include "novac/eph_manager.h" /* remove this line for use with solsys version 2 */
}

//==============================================================================//
// инициализация
//==============================================================================//
void CorrespondenceData::InitModyle()
{
	// модуль воздействий
	IForce = new Force::InfluenceForce();
	IForce->Init_CPU();

	//-----------------------------------------------------------------------//
	// эфемериды
	short int error = 0;
	short int de_num = 0;
	double jd_beg, jd_end;
	if ((error = ephem_open ("data\\eph\\lnx1900.405", &jd_beg,&jd_end,&de_num)) != 0)
	{
		if (error == 1)
			printf ("JPL ephemeris file not found.\n");
		else
			printf ("Error reading JPL ephemeris file header.\n");
	}
	else
	{
		printf ("JPL ephemeris DE%d open. Start JD = %10.2f  End JD = %10.2f\n",
			de_num, jd_beg, jd_end);
		printf ("\n");
	}
	//-----------------------------------------------------------------------//
};
//==============================================================================//
// тест соответствия с лазерными измерениями
//==============================================================================//
void CorrespondenceData::RunCorrespondenceData( char *optfname, char *htsName, char *tleName, double *Tpos )
{
	printf("\n**********************************************\n");

	std::cout << htsName << endl;
	std::cout << tleName << endl;
	std::cout << optfname << endl;

	// лазерные измерения
	cpfe.LoadCpffile( htsName );
	int id = cpfe.GetNoradId();

	// Каталог NORAD
	tleload.LoadData( tleName, 2 );
	orbitN = tleload.GetFromID( id );
	printf( "Get Sat = %s\n", orbitN->SatId().c_str() );
	printf( "Get Sat Name = %s\n", orbitN->SatName().c_str() );
	printf( "AP = %f PR = %f E = %f\n", orbitN->Apogee(), orbitN->Perigee(), orbitN->Eccentricity() );

	// измерения
	optLoad.LoadData( optfname, false );
	OpticArray = optLoad.GetOpticList();
	printf("**********************************************\n");

	// идем по точкам орбиты и выбираем видимые
	printf("----------------------------------------------------------\n");
	printf("\nNO CORRECT:\ndRa\tdDec\tPro\tNor\tZ\n");
	for( int it = 0; it < OpticArray.size(); it ++ )
	{
		CorrespondenceCPF( it, false, Tpos, false, false );
	}
	printf("----------------------------------------------------------\n");
	printf("\nABBERATION CORRECT:\ndRa\tdDec\tPro\tNor\tZ\n");
	for( int it = 0; it < OpticArray.size(); it ++ )
	{
		CorrespondenceCPF( it, false, Tpos, true, false );
	}
	
	printf("\nABBERATION AND TIME CORRECT:\ndRa\tdDec\tPro\tNor\tZ\n");
	for( int it = 0; it < OpticArray.size(); it ++ )
	{
		CorrespondenceCPF( it, true, Tpos, true, false );
	}

	printf("----------------------------------------------------------\n");
	printf("\nABBERATION, REFRACTION AND TIME CORRECT:\ndRa\tdDec\tPro\tNor\tZ\n");
	for( int it = 0; it < OpticArray.size(); it ++ )
	{
		CorrespondenceCPF( it, true, Tpos, true, true );
	}
	printf("----------------------------------------------------------\n");

	//TestCorrectRefract();

	return;
};
void CorrespondenceData::ITRFToICRF( double jd, double *posITRF, double *posICRF)
{
	DataConverter Dconv;
	// MDB
	jd = jd +0.125;
	double dataMDB = Dconv.JDtoYYYYMMDD(jd);
	double timeMDB = Dconv.SECtoHHMMSS(dataMDB, jd);

	// установка времени
	double int1, ajd1, delt1;
	IForce->set_time(dataMDB, timeMDB, &ajd1, &delt1, &int1);

	// матрица перевода в земную систему
	double Arot[9];
	IForce->iers_update_matrix(int1, Arot, ajd1, delt1);
	// матрица перехода из земной в нормальную систему
	double invArot[9];
	IForce->transpose(Arot, invArot);

	IForce->matVecMul(invArot, posITRF, posICRF);
}
//==============================================================================//
// проверка соответствия
//==============================================================================//
void CorrespondenceData::CorrespondenceCPF( int it, bool Tcorr, double *Tpos, bool abbcorr, bool refcorr )
{
	// отстование по времени
	double Toffset = 0;

	double pi = 3.1415926535;
	// установка времени
	double int1, ajd1, delt1;
	double date1 = OpticArray[it].dataMDB;
	double time1 = OpticArray[it].timeMDB;
	IForce->set_time(date1, time1, &ajd1, &delt1, &int1 );

	// матрица перевода в земную систему
	double Arot[9];
	IForce->iers_update_matrix( int1, Arot, ajd1, delt1 );
	// матрица перехода из земной в нормальную систему
	double invArot[9];
	IForce->transpose( Arot, invArot ); 

	// перевод в ICRF координат телескопа
	double Telicrf[3];
	//IForce->matVecMul(  invArot, Tpos, Telicrf );
	//Ter2Cel( jd, OpticArray[it].jd, Tpos, Telicrf );

	ITRFToICRF(OpticArray[it].jd, Tpos, Telicrf);

	// поправка на расстояние до спутника
	if( Tcorr )
	{
		double jdmeg = OpticArray[it].jd;
		// прогнозируем ТЛЕ на момент измерения м
		// прогноз по SGP4
		double Tpredict = jdmeg - orbitN->Epoch().Date();
		Tpredict = Tpredict*24.0*60.0;
		cEciTime eci1 = orbitN->GetPosition( Tpredict );

		double Ptle[3];
		// прогноз положения в километрах
		Ptle[0] = eci1.Position().m_x;
		Ptle[1] = eci1.Position().m_y;
		Ptle[2] = eci1.Position().m_z;

		TimePoint calcTime;
		calcTime.jd = OpticArray[it].jd;
		calcTime.gdata = OpticArray[it].dataMDB;
		calcTime.gtime = OpticArray[it].timeMDB;

		// положение спутника по TLE в момент измерения
		double PICRF[3];
		ConvertTEMEtoICRF( Ptle, PICRF, calcTime );

		double lx = PICRF[0] - Telicrf[0];
		double ly = PICRF[1] - Telicrf[1];
		double lz = PICRF[2] - Telicrf[2];
		double satDistTle = sqrt( lx*lx + ly*ly + lz*lz );

		Toffset = satDistTle/300.0; // задержка в мисисекундах
		//printf("\n\nTIME CORRECT = %.2f\tdist = %f\n", Toffset, satDistTle );
		Toffset = Toffset/1000.0;	// задержка в секундах
	}

	//-------------------------------------------------------//
	// матрица поворота в момент отражения света от спутника
	double int2, ajd2, delt2;
	double date2 = OpticArray[it].dataMDB;
	double time2 = OpticArray[it].timeMDB - Toffset;
	IForce->set_time(date2, time2, &ajd2, &delt2, &int2 );

	// матрица перевода в земную систему
	double Arot_m[9];
	IForce->iers_update_matrix( int2, Arot_m, ajd2, delt2 );
	// матрица перехода из земной в нормальную систему
	double invArot_m[9];
	IForce->transpose( Arot_m, invArot_m ); 
	//-------------------------------------------------------//

	//-------------------------------------------------------//
	// вычисление положения на орбите через интерполяцию
	// time
	double jd = IForce->dt_ajd( OpticArray[it].dataUTC ) - 0.5;
	double secday = (OpticArray[it].jd - jd )*86400.0;
	secday = secday - Toffset;

	double Sicrf[3];
	double Gpos[6];
	cpfe.GetPos( jd, secday, Gpos );
	IForce->matVecMul(  invArot_m, Gpos, Sicrf );
	//Ter2Cel( jd, OpticArray[it].jd, Gpos, Sicrf );

	double Sicrf2[3];
	double Gpos2[6];
	cpfe.GetPos( jd, secday+0.2, Gpos2 );
	IForce->matVecMul(  invArot_m, Gpos2, Sicrf2 );
	//Ter2Cel( jd, OpticArray[it].jd, Gpos2, Sicrf2 );
	//-------------------------------------------------------//

	// вычисление просчитанных измерений
	double Ra;
	double Dec;
	ConvertXYZtoRADEC( Sicrf, Telicrf, &Ra, &Dec );
	// вектор направления
	double pos[3];
	pos[0] = Sicrf[0] - Telicrf[0];
	pos[1] = Sicrf[1] - Telicrf[1];
	pos[2] = Sicrf[2] - Telicrf[2];
	//vector2radec( pos, &Ra, &Dec );
	//Ra = Ra/12.0*pi;
	//Dec = Dec/180.0*pi;

	// дальность до объекта
	double satDist = sqrt(pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]);
	double rDist = sqrt(Sicrf[0] * Sicrf[0] + Sicrf[1] * Sicrf[1] + Sicrf[2] * Sicrf[2]);
	//printf("dist = %f\t %f\n", satDist, rDist );

	// положение на секунду позже
	double Ra1;
	double Dec1;
	ConvertXYZtoRADEC( Sicrf2, Telicrf, &Ra1, &Dec1 );
	// вектор направления
	//pos[0] = Sicrf2[0] - Telicrf[0];
	//pos[1] = Sicrf2[1] - Telicrf[1];
	//pos[2] = Sicrf2[2] - Telicrf[2];
	//vector2radec( pos, &Ra1, &Dec1 );
	//Ra1 = Ra1/12.0*pi;
	//Dec1 = Dec1/180.0*pi;

	double optRa =  OpticArray[it].Ra;
	double optDec =  OpticArray[it].Dec;

	if( abbcorr == false )
	{
		//printf("ABS ERROR\n");
		CalcError( Ra, Dec, Ra1, Dec1, optRa, optDec, Telicrf );
	}
	else
	{
		//printf("\nABBERATION CORRECT:");
		double jd_utc = OpticArray[it].jd;
		double leap_secs = 35.0;
		double jd_tt = jd_utc + ((double) leap_secs + 32.184) / 86400.0;
		cat_entry sat;
		make_cat_entry ("DUMMY","xxx",0, optRa/pi*12.0, optDec/pi*180.0, 0.0, 0.0, 0.0, 0.0, &sat );
		double cra, cdec;
		short int res = virtual_star ( jd_tt, &sat, 0, &cra, &cdec);
		cra = cra/12.0*pi;
		cdec = cdec/180.0*pi;

		//double dra = optRa - cra;
		//double ddec = optDec - cdec;
		//double dd = sqrt( dra*dra + ddec*ddec );
		//printf( "%.2f sec\n", dd/pi*180.0*3600.0 );

		if( refcorr )
		{
			double outRa;
			double outDec;
			refractcorrect( cra, cdec, Telicrf, &outRa, &outDec, satDist, rDist );
			CalcError( Ra, Dec, Ra1, Dec1, outRa, outDec, Telicrf );
		}
		else
		{
			CalcError( Ra, Dec, Ra1, Dec1, cra, cdec, Telicrf );
		}
	}
	//-----------------------------------------------------------------------//
}
//==============================================================================//
// вычисление ошибок
//==============================================================================//
void CorrespondenceData::CalcError( double Ra, double Dec, double Ra1, double Dec1, double optRa, double optDec, double *Telicrf )
{
	double pi = 3.1415926535;
	double az = GetZenitAngle( optRa, optDec, Telicrf );

	Ra = Ra*cos(Dec);
	Ra1 = Ra1*cos(Dec1);
	optRa = optRa*cos(optDec);

	// вектор орбиты
	double AA1 = Ra1 - Ra;
	double AA2 = Dec1 - Dec;
	double AAabs = sqrt( AA1*AA1 + AA2*AA2 );
	// вектор от точки к измерению
	double BB1 = optRa - Ra;
	double BB2 = optDec - Dec;
	double omabs = sqrt( BB1*BB1 + BB2*BB2 );

	S3DCoordinate VO = S3DCoordinate( AA1, AA2, 0 );
	S3DCoordinate VM = S3DCoordinate( BB1, BB2, 0 );
	S3DCoordinate VZ = S3DCoordinate( 0, 0, 1 );

	// продольная ошибка
	double cosa = ( AA1*BB1 + AA2*BB2 )/AAabs/omabs;
	double perr = cosa*omabs;

	// поперечная ошибка
	double nerr = sqrt( omabs*omabs - perr*perr );

	// определение знака поперечной ошибки
	S3DCoordinate vecz = VO^VM;
	double sig = vecz*VZ/VZ.norm()/vecz.norm();
	if( sig < 0 ) sig = -1;
	if( sig > 0 ) sig = 1;
	nerr = nerr*sig;

	double dRa = Ra - optRa;
	double dDec = Dec - optDec;
	dRa = dRa/pi*180.0*3600.0;
	dDec = dDec/pi*180.0*3600.0;
	perr = perr/pi*180.0*3600.0;
	nerr = nerr/pi*180.0*3600.0;

	printf( "%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n", dRa, dDec, perr, nerr, az/pi*180.0 );
}
//==============================================================================//
// вычисление угла к зениту
//==============================================================================//
double CorrespondenceData::GetZenitAngle( double inRa, double inDec, double *tel_icrf)
{
	double pi = 3.1415926535;
	double tog = 180.0/pi;
	// измеренные значения координат телескопов, они искажены
	// double inRa, double inDec

	// вектор телескопа
	S3DCoordinate Rtel;
	Rtel.x = tel_icrf[0];
	Rtel.y = tel_icrf[1];
	Rtel.z = tel_icrf[2];

	// вектор наблюдения
	S3DCoordinate Rn;
	Rn.x = cos(inDec)*cos(inRa);
	Rn.y = cos(inDec)*sin(inRa);
	Rn.z = sin(inDec);

	double cosa = (Rtel*Rn )/( Rtel.norm()*Rn.norm() );
	double zenit = acos( cosa );

	//printf( "zenit Angle = %f, In Grad = %f \n", zenit, zenit/pi*180.0 );

	return zenit;
}
//==============================================================================//
// перевод спрогнозированного положения в координаты RA DEC
//==============================================================================//
void CorrespondenceData::ConvertXYZtoRADEC( double *resultPosition, double *inTelescopePosition, double *Ra, double *Dec )
{
	// входные данные задаются в ICRF
	// вектор направления в системе ICRF
	double x = resultPosition[0] - inTelescopePosition[0];
	double y = resultPosition[1] - inTelescopePosition[1];
	double z = resultPosition[2] - inTelescopePosition[2];

	double r = atan2( y, x );
	double d = atan2( z, sqrt( x*x + y*y ) );

	double pi = 3.1415926535;
	
	if( r < 0 )
		r = 2.0*pi + r;

	*Ra = r;
	*Dec = d;
}
//==============================================================================//
// получение положения телескопа в момент измерения в системе ICRF
//==============================================================================//
void CorrespondenceData::ConvertTEMEtoICRF( double *inPTEME, double *outPICRF, TimePoint calcTime )
{
	// установка времени
	double date1 = calcTime.gdata;
	double time1 = calcTime.gtime; 
	double int1, ajd1, delt1;
	IForce->set_time(date1, time1, &ajd1, &delt1, &int1 );

	// матрица перехода из ICRF в TEME
	double A_Teme[9];
	IForce->GetTemeMatrix( int1, A_Teme, ajd1, delt1 );
	// обратная матрица перехода из TEME в ICRF
	double invA_Teme[9];
	IForce->transpose( A_Teme, invA_Teme ); 

	// перевод вектора состояния в ICRF 
	IForce->matVecMul( invA_Teme, inPTEME, outPICRF );
}
//==============================================================================//