#include "CorrespondenceData.h"
#include "common/DataConverter.h"

extern "C" {
#include "novac/novas.h"
#include "novac/eph_manager.h" /* remove this line for use with solsys version 2 */
}

//==============================================================================//
// �������������
//==============================================================================//
void CorrespondenceData::InitModyle()
{
	// ������ �����������
	IForce = new Force::InfluenceForce();
	IForce->Init_CPU();

	//-----------------------------------------------------------------------//
	// ���������
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
// ���� ������������ � ��������� �����������
//==============================================================================//
void CorrespondenceData::RunCorrespondenceData( char *optfname, char *htsName, char *tleName, double *Tpos )
{
	printf("\n**********************************************\n");

	std::cout << htsName << endl;
	std::cout << tleName << endl;
	std::cout << optfname << endl;

	// �������� ���������
	cpfe.LoadCpffile( htsName );
	int id = cpfe.GetNoradId();

	// ������� NORAD
	tleload.LoadData( tleName, 2 );
	orbitN = tleload.GetFromID( id );
	printf( "Get Sat = %s\n", orbitN->SatId().c_str() );
	printf( "Get Sat Name = %s\n", orbitN->SatName().c_str() );
	printf( "AP = %f PR = %f E = %f\n", orbitN->Apogee(), orbitN->Perigee(), orbitN->Eccentricity() );

	// ���������
	optLoad.LoadData( optfname, false );
	OpticArray = optLoad.GetOpticList();
	printf("**********************************************\n");

	// ���� �� ������ ������ � �������� �������
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

	// ��������� �������
	double int1, ajd1, delt1;
	IForce->set_time(dataMDB, timeMDB, &ajd1, &delt1, &int1);

	// ������� �������� � ������ �������
	double Arot[9];
	IForce->iers_update_matrix(int1, Arot, ajd1, delt1);
	// ������� �������� �� ������ � ���������� �������
	double invArot[9];
	IForce->transpose(Arot, invArot);

	IForce->matVecMul(invArot, posITRF, posICRF);
}
//==============================================================================//
// �������� ������������
//==============================================================================//
void CorrespondenceData::CorrespondenceCPF( int it, bool Tcorr, double *Tpos, bool abbcorr, bool refcorr )
{
	// ���������� �� �������
	double Toffset = 0;

	double pi = 3.1415926535;
	// ��������� �������
	double int1, ajd1, delt1;
	double date1 = OpticArray[it].dataMDB;
	double time1 = OpticArray[it].timeMDB;
	IForce->set_time(date1, time1, &ajd1, &delt1, &int1 );

	// ������� �������� � ������ �������
	double Arot[9];
	IForce->iers_update_matrix( int1, Arot, ajd1, delt1 );
	// ������� �������� �� ������ � ���������� �������
	double invArot[9];
	IForce->transpose( Arot, invArot ); 

	// ������� � ICRF ��������� ���������
	double Telicrf[3];
	//IForce->matVecMul(  invArot, Tpos, Telicrf );
	//Ter2Cel( jd, OpticArray[it].jd, Tpos, Telicrf );

	ITRFToICRF(OpticArray[it].jd, Tpos, Telicrf);

	// �������� �� ���������� �� ��������
	if( Tcorr )
	{
		double jdmeg = OpticArray[it].jd;
		// ������������ ��� �� ������ ��������� �
		// ������� �� SGP4
		double Tpredict = jdmeg - orbitN->Epoch().Date();
		Tpredict = Tpredict*24.0*60.0;
		cEciTime eci1 = orbitN->GetPosition( Tpredict );

		double Ptle[3];
		// ������� ��������� � ����������
		Ptle[0] = eci1.Position().m_x;
		Ptle[1] = eci1.Position().m_y;
		Ptle[2] = eci1.Position().m_z;

		TimePoint calcTime;
		calcTime.jd = OpticArray[it].jd;
		calcTime.gdata = OpticArray[it].dataMDB;
		calcTime.gtime = OpticArray[it].timeMDB;

		// ��������� �������� �� TLE � ������ ���������
		double PICRF[3];
		ConvertTEMEtoICRF( Ptle, PICRF, calcTime );

		double lx = PICRF[0] - Telicrf[0];
		double ly = PICRF[1] - Telicrf[1];
		double lz = PICRF[2] - Telicrf[2];
		double satDistTle = sqrt( lx*lx + ly*ly + lz*lz );

		Toffset = satDistTle/300.0; // �������� � ������������
		//printf("\n\nTIME CORRECT = %.2f\tdist = %f\n", Toffset, satDistTle );
		Toffset = Toffset/1000.0;	// �������� � ��������
	}

	//-------------------------------------------------------//
	// ������� �������� � ������ ��������� ����� �� ��������
	double int2, ajd2, delt2;
	double date2 = OpticArray[it].dataMDB;
	double time2 = OpticArray[it].timeMDB - Toffset;
	IForce->set_time(date2, time2, &ajd2, &delt2, &int2 );

	// ������� �������� � ������ �������
	double Arot_m[9];
	IForce->iers_update_matrix( int2, Arot_m, ajd2, delt2 );
	// ������� �������� �� ������ � ���������� �������
	double invArot_m[9];
	IForce->transpose( Arot_m, invArot_m ); 
	//-------------------------------------------------------//

	//-------------------------------------------------------//
	// ���������� ��������� �� ������ ����� ������������
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

	// ���������� ������������ ���������
	double Ra;
	double Dec;
	ConvertXYZtoRADEC( Sicrf, Telicrf, &Ra, &Dec );
	// ������ �����������
	double pos[3];
	pos[0] = Sicrf[0] - Telicrf[0];
	pos[1] = Sicrf[1] - Telicrf[1];
	pos[2] = Sicrf[2] - Telicrf[2];
	//vector2radec( pos, &Ra, &Dec );
	//Ra = Ra/12.0*pi;
	//Dec = Dec/180.0*pi;

	// ��������� �� �������
	double satDist = sqrt(pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]);
	double rDist = sqrt(Sicrf[0] * Sicrf[0] + Sicrf[1] * Sicrf[1] + Sicrf[2] * Sicrf[2]);
	//printf("dist = %f\t %f\n", satDist, rDist );

	// ��������� �� ������� �����
	double Ra1;
	double Dec1;
	ConvertXYZtoRADEC( Sicrf2, Telicrf, &Ra1, &Dec1 );
	// ������ �����������
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
// ���������� ������
//==============================================================================//
void CorrespondenceData::CalcError( double Ra, double Dec, double Ra1, double Dec1, double optRa, double optDec, double *Telicrf )
{
	double pi = 3.1415926535;
	double az = GetZenitAngle( optRa, optDec, Telicrf );

	Ra = Ra*cos(Dec);
	Ra1 = Ra1*cos(Dec1);
	optRa = optRa*cos(optDec);

	// ������ ������
	double AA1 = Ra1 - Ra;
	double AA2 = Dec1 - Dec;
	double AAabs = sqrt( AA1*AA1 + AA2*AA2 );
	// ������ �� ����� � ���������
	double BB1 = optRa - Ra;
	double BB2 = optDec - Dec;
	double omabs = sqrt( BB1*BB1 + BB2*BB2 );

	S3DCoordinate VO = S3DCoordinate( AA1, AA2, 0 );
	S3DCoordinate VM = S3DCoordinate( BB1, BB2, 0 );
	S3DCoordinate VZ = S3DCoordinate( 0, 0, 1 );

	// ���������� ������
	double cosa = ( AA1*BB1 + AA2*BB2 )/AAabs/omabs;
	double perr = cosa*omabs;

	// ���������� ������
	double nerr = sqrt( omabs*omabs - perr*perr );

	// ����������� ����� ���������� ������
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
// ���������� ���� � ������
//==============================================================================//
double CorrespondenceData::GetZenitAngle( double inRa, double inDec, double *tel_icrf)
{
	double pi = 3.1415926535;
	double tog = 180.0/pi;
	// ���������� �������� ��������� ����������, ��� ��������
	// double inRa, double inDec

	// ������ ���������
	S3DCoordinate Rtel;
	Rtel.x = tel_icrf[0];
	Rtel.y = tel_icrf[1];
	Rtel.z = tel_icrf[2];

	// ������ ����������
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
// ������� ������������������ ��������� � ���������� RA DEC
//==============================================================================//
void CorrespondenceData::ConvertXYZtoRADEC( double *resultPosition, double *inTelescopePosition, double *Ra, double *Dec )
{
	// ������� ������ �������� � ICRF
	// ������ ����������� � ������� ICRF
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
// ��������� ��������� ��������� � ������ ��������� � ������� ICRF
//==============================================================================//
void CorrespondenceData::ConvertTEMEtoICRF( double *inPTEME, double *outPICRF, TimePoint calcTime )
{
	// ��������� �������
	double date1 = calcTime.gdata;
	double time1 = calcTime.gtime; 
	double int1, ajd1, delt1;
	IForce->set_time(date1, time1, &ajd1, &delt1, &int1 );

	// ������� �������� �� ICRF � TEME
	double A_Teme[9];
	IForce->GetTemeMatrix( int1, A_Teme, ajd1, delt1 );
	// �������� ������� �������� �� TEME � ICRF
	double invA_Teme[9];
	IForce->transpose( A_Teme, invA_Teme ); 

	// ������� ������� ��������� � ICRF 
	IForce->matVecMul( invA_Teme, inPTEME, outPICRF );
}
//==============================================================================//