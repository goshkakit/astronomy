#pragma once

#include <iostream>
#include <chrono>
#include <ctime>  
#include <vector>

// System 
#include "InfluenceForce\InfluenceForce.h"
#include "AstroDataManager.h"
#include "AstroTypes.h"

#include "svd/ap.h"
#include "svd/svd.h"

extern "C" {
#include "novas.h"
#include "eph_manager.h" /* remove this line for use with solsys version 2 */
}

class TelescopeAling
{
private:
	DataConverter Dconv_;
	AstroConstant const_;
	Force::InfluenceForce *IForce;

	//double tel_height = 0.0;
	//double tel_latitude = 55.80364;
	//double tel_longitude = 37.54794;

public:
	TelescopeAling() {}

	void Init()
	{
		// init rotation
		IForce = new Force::InfluenceForce();
		IForce->Init_CPU();

		//-----------------------------------------------------------------------//
		// эфемериды
		short int error = 0;
		short int de_num = 0;
		double jd_beg, jd_end;
		if ((error = ephem_open("data\\eph\\lnx1900.405", &jd_beg, &jd_end, &de_num)) != 0)
		{
			if (error == 1)
				printf("JPL ephemeris file not found.\n");
			else
				printf("Error reading JPL ephemeris file header.\n");
		}
		else
		{
			printf("JPL ephemeris DE%d open. Start JD = %10.2f  End JD = %10.2f\n",
				de_num, jd_beg, jd_end);
			printf("\n");
		}
		//-----------------------------------------------------------------------//
	}

	void getCurrentTime(TelescopeDataTime &dtime)
	{
		auto tm = std::chrono::system_clock::now();
		std::time_t td = std::chrono::system_clock::to_time_t(tm);
		struct tm *aTime = localtime(&td);

		int day = aTime->tm_mday;
		int month = aTime->tm_mon + 1; // Month is 0 - 11, add 1 to get a jan-dec 1-12 concept
		int year = aTime->tm_year + 1900; // Year is # years since 1900

		int hh = aTime->tm_hour;
		int mm = aTime->tm_min;
		int ss = aTime->tm_sec;

		dtime.Y = year;
		dtime.M = month;
		dtime.D = day;

		dtime.hh = hh;
		dtime.mm = mm;
		dtime.ss = ss;
		dtime.ms = 0;
	}

	//void getSystemTime()
	//{
	//	SYSTEMTIME st;
	//	GetSystemTime(&st);
	//	WORD wYear;
	//	WORD wMonth;
	//	WORD wDayOfWeek;
	//	WORD wDay;
	//	WORD wHour;
	//	WORD wMinute;
	//	WORD wSecond;
	//	WORD wMilliseconds;
	//}

	void ITRFToICRF(double jd, double *posITRF, double *posICRF)
	{
		// MDB
		jd = jd + 0.125;
		double dataMDB = Dconv_.JDtoYYYYMMDD(jd);
		double timeMDB = Dconv_.SECtoHHMMSS(dataMDB, jd);
		printf("ITRFToICRF: Data MDB: %f Time MDB: %f", dataMDB, timeMDB);

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


	void ConvertToICRF(double jd,  S3DCoordinate vItrf, S3DCoordinate &vIcrf)
	{
		// MDB
		jd = jd + 0.125;
		double date1 = Dconv_.JDtoYYYYMMDD(jd);
		double time1 = Dconv_.SECtoHHMMSS(date1, jd);

		// установка времени
		double int1, ajd1, delt1;
		IForce->set_time(date1, time1, &ajd1, &delt1, &int1);

		// матрица перевода в земную систему
		double Arot[9];
		IForce->iers_update_matrix(int1, Arot, ajd1, delt1);
		// матрица перехода из земной в нормальную систему
		double invArot[9];
		IForce->transpose(Arot, invArot);

		// перевод в ICRF координат телескопа
		double Telicrf[3];
		double Tpos[3];
		Tpos[0] = vItrf.x;
		Tpos[1] = vItrf.y;
		Tpos[2] = vItrf.z;
		IForce->matVecMul(invArot, Tpos, Telicrf);
		vIcrf.x = Telicrf[0];
		vIcrf.y = Telicrf[1];
		vIcrf.z = Telicrf[2];
	}

	std::vector<SCoordinate> getListPointAzEl(int nAz, int nEl, double minEl )
	{
		std::vector<SCoordinate> list;

		for (int j = 0; j < nEl; j++)
		{
			for (int i = 0; i < nAz; i++)
			{
				double az = 360.0 / (double)nAz*(double(i));
				double el = minEl + (90.0 - minEl) / (double)nEl*(double(j));
				SCoordinate pt = SCoordinate(az, el);
				list.push_back(pt);
			}
		}
		return list;
	}

	void ConvertXYZtoRADEC(S3DCoordinate pos, S3DCoordinate tel, SCoordinate &RA_DEC)
	{
		// входные данные задаются в ICRF
		// вектор направления в системе ICRF
		double x = pos.x - tel.x;
		double y = pos.y - tel.y;
		double z = pos.z - tel.z;

		double r = atan2(y, x);
		double d = atan2(z, sqrt(x*x + y*y));

		if (r < 0)
			r = 2.0*M_PI + r;

		RA_DEC.x = r;
		RA_DEC.y = d;
	}

	void ConvertXYZtoRADEC(S3DCoordinate vec, SCoordinate &RA_DEC)
	{
		// входные данные задаются в ICRF
		// вектор направления в системе ICRF
		double x = vec.x;
		double y = vec.y;
		double z = vec.z;

		double r = atan2(y, x);
		double d = atan2(z, sqrt(x*x + y*y));

		if (r < 0)
			r = 2.0*M_PI + r;

		RA_DEC.x = r;
		RA_DEC.y = d;
	}

	void ConvertTEMEtoICRF(double *inPTEME, double *outPICRF, double date1, double time1)
	{
		// установка времени
		double int1, ajd1, delt1;
		IForce->set_time(date1, time1, &ajd1, &delt1, &int1);

		// матрица перехода из ICRF в TEME
		double A_Teme[9];
		IForce->GetTemeMatrix(int1, A_Teme, ajd1, delt1);
		// обратная матрица перехода из TEME в ICRF
		double invA_Teme[9];
		IForce->transpose(A_Teme, invA_Teme);

		// перевод вектора состояния в ICRF 
		IForce->matVecMul(invA_Teme, inPTEME, outPICRF);
	}

	S3DCoordinate SphericalToXYZ_Deg(double lon, double lat, double r)
	{
		S3DCoordinate pt;
		pt.x = r*cos(lat*const_.GR)*cos(lon*const_.GR);
		pt.y = r*cos(lat*const_.GR)*sin(lon*const_.GR);
		pt.z = r*sin(lat*const_.GR);
		return pt;
	};

	S3DCoordinate SphericalToXYZ_Rad(double lon, double lat, double r)
	{
		S3DCoordinate pt;
		pt.x = r*cos(lat)*cos(lon);
		pt.y = r*cos(lat)*sin(lon);
		pt.z = r*sin(lat);
		return pt;
	};

	int convertXYZtoAzimut(S3DCoordinate vT, S3DCoordinate vO, SCoordinate &AzEl)
	{
		// calc position
		// calc angle

		// вектор от телескопа на объект
		S3DCoordinate vTO = vO - vT;

		// угол места относительно горизонта
		double anglea = vT*vTO / (vT.norm()*vTO.norm());
		if (anglea > 1.0)
		{
			anglea = 1.0;
		}
		double Am = acos(anglea);
		Am = M_PI / 2.0 - Am;

		// вертикальная составляющая вектора от телескопа на объект
		S3DCoordinate vOver = vT*(vTO.norm()*sin(Am)) / vT.norm();

		// горизонтпльная составляющая
		S3DCoordinate vOhor = vTO - vOver;
		double Az = 0.0;
		if(vOhor.norm() > 0.0)
		{
			// нормально к плоскости проходящей через ось  Z и вектор на телескоп
			S3DCoordinate vZ = S3DCoordinate(0, 0, 1);
			S3DCoordinate vN = vT^vZ;

			// вектор направленный на север
			S3DCoordinate vA = vN^vT;

			// азимут
			double cosaAz = vA*vOhor / (vA.norm()*vOhor.norm());
			if (cosaAz > 1.0)
			{
				cosaAz = 1.0;
			}

			Az = acos(cosaAz);

			// учет полушария
			double signAz = vN*vOhor / (vN.norm()*vOhor.norm());
			if (signAz > 0 && Az > 0.0 )
				Az = 2.0*M_PI - Az;
		}

		AzEl.x = Az;
		AzEl.y = Am;

		return 0;
	}

	S3DCoordinate getVectorAzEl(S3DCoordinate tel, SCoordinate AzEl)
	{
		double Az = AzEl.x*const_.GR;
		double El = AzEl.y*const_.GR;

		// положение старта
		S3DCoordinate startpos = tel;

		// радиус вектор, перпендикулр к поверхности
		S3DCoordinate vec_n = startpos / startpos.norm();

		// вектор вдоль оси Z
		S3DCoordinate vec_z = S3DCoordinate(0.0, 0.0, 1.0);

		// вектор вдоль широты
		S3DCoordinate vec_lat = vec_z^vec_n;
		vec_lat = vec_lat / vec_lat.norm();

		// вектор вдоль долготы
		S3DCoordinate vec_lon = vec_n^vec_lat;
		vec_lon = vec_lon / vec_lon.norm();

		// вектор в направлении азимута
		S3DCoordinate AzimutVector = vec_lon*cos(Az) + vec_lat*sin(Az);
		AzimutVector = AzimutVector / AzimutVector.norm();

		// вектор в направлении угламеста
		S3DCoordinate ElVector = AzimutVector*cos(El) + vec_n*sin(El);
		ElVector = ElVector / ElVector.norm();

		SCoordinate AzEl_Test;
		convertXYZtoAzimut(startpos, startpos + ElVector, AzEl_Test);
		printf("Input Az= %f\t El= %f\t Test: Az= %f\t Pla= %f\n", Az*const_.RG, El*const_.RG, AzEl_Test.x*const_.RG, AzEl_Test.y*const_.RG);

		return ElVector;
	}

	double GetElAngle(double inRa, double inDec, S3DCoordinate telIcrf)
	{
		// вектор наблюдения
		S3DCoordinate Rn;
		Rn.x = cos(inDec)*cos(inRa);
		Rn.y = cos(inDec)*sin(inRa);
		Rn.z = sin(inDec);

		double cosa = (telIcrf*Rn) / (telIcrf.norm()*Rn.norm());
		if (cosa > 1.0)
		{
			cosa = 1.0;
		}
		double zenit = acos(cosa);
		double El = M_PI / 2.0 - zenit;
		return El;
	}

	void Ter2Cel(double jd_utc, S3DCoordinate PxydUt, S3DCoordinate &vec1, S3DCoordinate &vec2, SCoordinate &posRaDec )
	{
		// поправки полюса
		//13 1 1 56293.00 P  0.082004 0.005865  0.298188 0.007111  P 0.2857338 0.0055637                 P   -75.505     .600    -8.527     .600                                                     
		//double xp = 0.082004;
		//double yp = 0.298188;
		//double ut1_utc = 0.2857338;

		double xp = PxydUt.x;
		double yp = PxydUt.y;
		double ut1_utc = PxydUt.z;
		double leap_secs = 37.0;
		printf("xp= %f\t yp= %f\t dUt1= %f\n", xp, yp, ut1_utc);

		double vter[3], vcel[3];
		// input
		vter[0] = vec1.x;
		vter[1] = vec1.y;
		vter[2] = vec1.z;

		// ∆AT is an integer representing the total count of leap seconds in UTC
		// TT = UTC + ∆AT + 32.184s
		// TT = TAI + 32.184s
		// TAI is International Atomic Time
		double jd_tt = jd_utc + ((double)leap_secs + 32.184) / 86400.0;
		// UT1 = UTC + (UT1–UTC)
		double jd_ut1 = jd_utc + ut1_utc / 86400.0;
		// UT1 = TT – ∆T
		// ∆T = TT - UT1 = UTC + ∆AT + 32.184s - UTC - (UT1–UTC)
		// ∆T = ∆AT + 32.184s - (UT1–UTC)
		double delta_t = 32.184 + leap_secs - ut1_utc;

		short int method = 0;
		short int accuracy = 0;
		short int option = 0;

		//method(short int)
		//Selection for method
		//= 0 ... CIO - based method
		//= 1 ... equinox - based method
		//accuracy(short int)
		//Selection for accuracy
		//= 0 ... full accuracy
		//= 1 ... reduced accuracy
		//option(short int)
		//= 0 ... The output vector is referred to GCRS axes.
		//= 1 ... The output vector is produced with respect to the
		//equator and equinox of date.
		//The 'option' flag only works for the equinox - based method.

		// error = ter2cel (jd_ut1,0.0,delta_t,1,accuracy,0,x_pole,y_pole,vter,vcel)
		short int res = ter2cel(jd_ut1, 0, delta_t, method, accuracy, option, xp, yp, vter, vcel);
		if (res != 0)
		{
			printf("ter2cel res= %d\n", res );
		}

		// result
		vec2.x = vcel[0];
		vec2.y = vcel[1];
		vec2.z = vcel[2];

		double ra, dec;
		short int error = 0;
		if ((error = vector2radec(vcel, &ra, &dec)) != 0)
		{
			printf("Error %d from vector2radec.", error);
		}
		posRaDec.x = ra / 12.0*M_PI;
		posRaDec.y = dec / 180.0*M_PI;
	}

	// input
	//	ra,dec in rad
	// output
	//	az.el in deg
	void EquToHor(double jd_utc, double ra, double dec, S3DCoordinate tel, S3DCoordinate PxydUT, SCoordinate &AzEl )
	{
		ra = ra / M_PI*12.0;
		dec = dec / M_PI*180.0;

		//double ra, dec;
		//double pos[3];
		//// input
		//pos[0] = vec.x;
		//pos[1] = vec.y;
		//pos[2] = vec.z;
		//vector2radec(pos, &ra, &dec);

		double xp = PxydUT.x;
		double yp = PxydUT.y;
		double ut1_utc = PxydUT.z;
		double leap_secs = 37.0;

		double jd_ut1 = jd_utc + ut1_utc / 86400.0;
		double jd_tt = jd_utc + ((double)leap_secs + 32.184) / 86400.0;
		double delta_t = 32.184 + leap_secs - ut1_utc;
		short int accuracy = 0;

		SCoordinate vangle;
		ConvertXYZtoRADEC(tel, vangle);
		vangle = vangle*const_.RG;

		on_surface location;
		location.height = 0.0;
		location.latitude = vangle.y;
		location.longitude = vangle.x;
		//printf("Tel Position calc: lon= %f lat= %f\n", vangle.x, vangle.y);
		//location.latitude = tel_latitude;
		//location.longitude = tel_longitude;
		printf("Tel Position: lon= %f lat= %f\n", location.longitude, location.latitude);

		short int ref_option = 0;

		//short int gcrs2equ(double jd_tt, short int coord_sys,short int accuracy, double rag, double decg, double *ra, double *dec);
		short int coord_sys = 1;
		gcrs2equ(jd_tt, coord_sys, accuracy, ra, dec, &ra, &dec);

		double zd, az, rar, decr;
		equ2hor(jd_ut1, delta_t, accuracy, xp, yp, &location, ra, dec, ref_option, &zd, &az, &rar, &decr);
		//equ2hor(jd_ut1, delta_t, accuracy, 0, 0, &location, ra, dec, ref_option, &zd, &az, &rar, &decr);

		AzEl.x = az;
		AzEl.y = 90.0 - zd;
	};

	std::vector<SCoordinate> getPointRaDec(S3DCoordinate tel, std::vector<SCoordinate> &listAzEl, std::vector<double> &listJD, std::vector<S3DCoordinate> &listVectorAzEl, std::vector<S3DCoordinate> &listPxydUT)
	{
		std::vector<SCoordinate> listRaDec;
		for (int i = 0; i < listAzEl.size(); i++)
		{
			printf("%d\n", i);
			// Вычисление вектора в направлении азимута и угла места
			S3DCoordinate VectorAzEl = getVectorAzEl(tel, listAzEl[i]);
			listVectorAzEl.push_back(VectorAzEl);

			// перевод в ICRF
			double JD = listJD[i];
			S3DCoordinate vItrf = VectorAzEl;
			S3DCoordinate vIcrf;
			ConvertToICRF(JD, vItrf, vIcrf);

			S3DCoordinate PxydUT;
			PxydUT.x = IForce->current_px_;
			PxydUT.y = IForce->current_py_;
			PxydUT.z = IForce->current_dUt1_;
			listPxydUT.push_back(PxydUT);	

			S3DCoordinate vTelIcrf;
			ConvertToICRF(JD, tel, vTelIcrf);

			// other convert
			// Az:El to Ra:Dec
			S3DCoordinate vIcrf_n;
			SCoordinate RaDec_n;
			Ter2Cel(JD, PxydUT, vItrf, vIcrf_n, RaDec_n);

			vIcrf_n = vIcrf_n / vIcrf_n.norm();
			vIcrf = vIcrf / vIcrf.norm();

			S3DCoordinate dPos = (vIcrf_n - vIcrf);
			printf("Pos %f %f %f dPos: %e %e %e\n", vIcrf.x, vIcrf.y, vIcrf.z, dPos.x, dPos.y, dPos.z);

			// Tel_ter to Tel_cel
			S3DCoordinate vTelIcrf_n;
			SCoordinate vTelRaDec_n;
			Ter2Cel(JD, PxydUT, tel, vTelIcrf_n, vTelRaDec_n);
			dPos = (vTelIcrf_n - vTelIcrf)*1000.0;
			printf("Pos [km] %f %f %f dPos [m]: %e %e %e\n", tel.x, tel.y, tel.z, dPos.x, dPos.y, dPos.z);

			// вычисление Ra Dec
			SCoordinate RaDec;
			ConvertXYZtoRADEC(vIcrf, RaDec);
			listRaDec.push_back(RaDec);

			SCoordinate dRaDec = (RaDec_n - RaDec)*const_.RG*3600.0;

			double el = GetElAngle(RaDec.x, RaDec.y, vTelIcrf);

			SCoordinate AzEl_n;
			EquToHor(JD, RaDec_n.x, RaDec_n.y, tel, PxydUT, AzEl_n);

			SCoordinate dAzEl = AzEl_n - listAzEl[i];

			printf("Ra= %f\t [%e asec]\t Dec= %f\t [%e asec]\t El_verify= %f\n", RaDec.x*const_.RG, dRaDec.x, RaDec.y*const_.RG, dRaDec.y, el*const_.RG);
			printf("Input: %f %f\t NOVAS: %f %f\t Err: %f %f\n", listAzEl[i].x, listAzEl[i].y, AzEl_n.x, AzEl_n.y, dAzEl.x, dAzEl.y);
		}
		return listRaDec;
	}

	// ICP
	int shiftMaskToCloud(std::vector< S3DCoordinate> &PointsSrc, std::vector< S3DCoordinate> &PointsDst, S3DCoordinate &totalT, S3DMatrix &totalR)
	{
		for (unsigned int it = 0; it < PointsSrc.size(); it++)
		{
			S3DCoordinate pt;
			pt = PointsSrc[it];
			pt = totalR*pt + totalT;
			PointsDst[it] = pt;
		}
		return 0;
	}

	float calculateDistanse(std::vector< S3DCoordinate> &PointsSrc, std::vector< S3DCoordinate> &PointsDst)
	{
		float dist = 0;
		int count = 0;
		for (unsigned int it = 0; it < PointsSrc.size(); it++)
		{
			S3DCoordinate pt1 = PointsSrc[it];
			S3DCoordinate pt2 = PointsDst[it];
			dist += (pt1 - pt2).norm();
			count++;
		}
		dist = dist / count;
		return dist;
	}

	float getShiftRotationFromICP(std::vector< S3DCoordinate> &PointsSrc, std::vector< S3DCoordinate> &PointsDst, S3DCoordinate &totalTn, S3DMatrix &totalRn)
	{
		unsigned int numSrc = PointsSrc.size();

		using namespace ap;
		S3DCoordinate fcenter(0, 0, 0);
		S3DCoordinate ccenter(0, 0, 0);

		float TotalDist = 0;
		S3DCoordinate Pm;
		S3DCoordinate Pc;
		int closestPointsN = 0;

		for (unsigned int i = 0; i < numSrc; i++)
		{
			Pm = PointsSrc[i];
			Pc = PointsDst[i];
			float pdist = (Pm - Pc).norm();
			closestPointsN++;
			TotalDist += pdist;
			fcenter = fcenter + Pm;
			ccenter = ccenter + Pc;
		}

		if (closestPointsN == 0)
		{
			printf("Error\n");
			return -1;
		}

		fcenter = fcenter / (float)closestPointsN;
		ccenter = ccenter / (float)closestPointsN;
		TotalDist = TotalDist / (float)closestPointsN;

		ap::real_2d_array Cov3;
		ap::real_1d_array S;
		ap::real_2d_array U, Vt;

		Cov3.setbounds(1, 3, 1, 3);
		for (int i = 1; i <= 3; i++)
		{
			for (int j = 1; j <= 3; j++)
			{
				Cov3(i, j) = 0;
			}
		}

		float Yx, Yy, Yz, Xx, Xy, Xz;
		for (unsigned int i = 0; i < numSrc; i++)
		{
			Yx = PointsSrc[i].x - fcenter.x;
			Yy = PointsSrc[i].y - fcenter.y;
			Yz = PointsSrc[i].z - fcenter.z;

			Xx = PointsDst[i].x - ccenter.x;
			Xy = PointsDst[i].y - ccenter.y;
			Xz = PointsDst[i].z - ccenter.z;

			Cov3(1, 1) += Xx*Yx;
			Cov3(2, 2) += Xy*Yy;
			Cov3(3, 3) += Xz*Yz;

			Cov3(1, 2) += Xx*Yy;
			Cov3(2, 1) += Xy*Yx;

			Cov3(1, 3) += Xx*Yz;
			Cov3(3, 1) += Xz*Yx;

			Cov3(2, 3) += Xy*Yz;
			Cov3(3, 2) += Xz*Yy;
		}

		// svd
		svddecomposition(Cov3, 3, 3, 2, 2, 2, S, U, Vt);

		// rotation
		S3DMatrix RCur;
		for (int i = 1; i <= 3; i++) {
			for (int j = 1; j <= 3; j++) {
				for (int k = 1; k <= 3; k++) {
					RCur.a[i - 1][k - 1] += U(i, j)*Vt(j, k);
				}
			}
		}

		// translate
		S3DCoordinate TCur = ccenter - RCur*fcenter;

		// Total R,T
		totalTn = TCur + (RCur*totalTn);
		totalRn = RCur*totalRn;

		return TotalDist;
	}

	void estimateRT(std::vector<S3DCoordinate> &Ref, std::vector<S3DCoordinate> &Mount, S3DMatrix &Rmat)
	{
		std::vector < S3DCoordinate > S1 = Ref;
		std::vector < S3DCoordinate > S1tmp = S1;
		std::vector < S3DCoordinate > S2 = Mount;


		// example rotation and offset
		//S3DCoordinate T = S3DCoordinate(-3, 1, 2);
		//S3DMatrix R = S3DMatrix::E();
		//double angle = 25.0 / 180.0*M_PI;
		//S3DMatrix M = R.R12(angle);
		//shiftMaskToCloud(S2, S2, T, M);
		//for (unsigned int it = 0; it < S2.size(); it++)
		//{
		//	S2[it] = S2[it] + S3DCoordinate(0.1*(float)rand()/(float)RAND_MAX, 0.1*(float)rand() / (float)RAND_MAX, 0.1*(float)rand() / (float)RAND_MAX);
		//}
		//printf("Tset  %f %f %f\n", T.x, T.y, T.z);
		//printf("Mset\n");
		//M.print();

		unsigned int NUMBIT = 3;
		S3DCoordinate Tn = S3DCoordinate(0, 0, 0);
		S3DMatrix Mn = S3DMatrix::E();
		for (unsigned int it = 0; it < NUMBIT; it++)
		{
			float d = calculateDistanse(S1tmp, S2);
			printf("D = %f\n", d);
			getShiftRotationFromICP(S1tmp, S2, Tn, Mn);
			shiftMaskToCloud(S1, S1tmp, Tn, Mn);
		}
		float d = calculateDistanse(S1tmp, S2);
		printf("D = %f\n", d);
		printf("Tfind %f %f %f\n", Tn.x, Tn.y, Tn.z);


		printf("Mfind det = %f\n", Mn.GetDET());
		Mn.print();

		Rmat = Mn;
	}

	S3DMatrix getRMatr(std::vector<SCoordinate> p_listAzEl_Tel, std::vector<SCoordinate> p_listAzEl_Real)
	{
		std::vector<S3DCoordinate> Ref_Real;
		for (int i = 0; i < p_listAzEl_Real.size(); i++)
		{
			SCoordinate pt = p_listAzEl_Real[i];
			S3DCoordinate pt3 = SphericalToXYZ_Deg(pt.x, pt.y, 1.0);
			Ref_Real.push_back(pt3);
		}

		std::vector<S3DCoordinate> Ref_Tel;
		for (int i = 0; i < p_listAzEl_Tel.size(); i++)
		{
			SCoordinate pt = p_listAzEl_Tel[i];
			S3DCoordinate pt3 = SphericalToXYZ_Deg(pt.x, pt.y, 1.0);
			Ref_Tel.push_back(pt3);
		}

		S3DMatrix Rmat;
		estimateRT(Ref_Real, Ref_Tel, Rmat);
		return Rmat;
	}
	
	SCoordinate alignAzEl(SCoordinate p_azel_mount, S3DMatrix p_rmat)
	{
		S3DCoordinate viewInMount = p_rmat.inverse()*SphericalToXYZ_Deg(p_azel_mount.x, p_azel_mount.y, 1.0);
		SCoordinate azel_view;
		ConvertXYZtoRADEC(viewInMount, azel_view);

		return azel_view*const_.RG;
	}

	void calculatePoints()
	{
		//######################//
		// Set Time and location

		TelescopeDataTime dtime;
		getCurrentTime(dtime);
		printf("Data MDB: %f\n", dtime.getYYYYMMDD());
		printf("Time MDB: %f\n", dtime.getHHMMSS());

		// MDB
		double dateMDB = dtime.getYYYYMMDD();
		double timeMDB = dtime.getHHMMSS();

		double JD = Dconv_.YYYYMMDDtoJD(dateMDB) - 0.5;
		JD += dtime.getSecOfDay()/86400.0 - 0.125;
		printf("JD: %f\n", JD);

		double telpos[3];
		telpos[0] = 2848.74429;
		telpos[1] = 2189.68052;
		telpos[2] = 5252.33756;
		S3DCoordinate tel_itrf = S3DCoordinate(telpos[0], telpos[1], telpos[2]);
		//S3DCoordinate tel_itrf = SphericalToXYZ_Deg(tel_longitude, tel_latitude, 6371.0);

		//######################//
		// Calculate pattern

		// point for view [Tel ITRC System Coordinate]
		
		std::vector<SCoordinate> listAzEl_Tel = getListPointAzEl(5, 2, 20.0);
		// Вычисление вектора в направлении азимута и угла места
		std::vector<S3DCoordinate> listVectorAzEl_Tel;
		for (int i = 0; i < listAzEl_Tel.size(); i++)
		{
			S3DCoordinate VectorAzEl_Tel = getVectorAzEl(tel_itrf, listAzEl_Tel[i]);
			listVectorAzEl_Tel.push_back(VectorAzEl_Tel);
		}

		// time for view
		std::vector<double> listJD;
		for (int i = 0; i < listAzEl_Tel.size(); i++)
		{
			listJD.push_back(JD);
		}

		// Simulation
		// offset telescope orientation
		std::vector<SCoordinate> listAzEl_Real;
		for (int i = 0; i < listAzEl_Tel.size(); i++)
		{
			SCoordinate pt = listAzEl_Tel[i];
			pt.x += 10;
			pt.y += 2;
			listAzEl_Real.push_back(pt);
		}

		// ra dec for view or get result from frames [ICRS System coordinate]
		std::vector<S3DCoordinate> listVectorAzEl_Real;
		std::vector<S3DCoordinate> listVectorRaDec;
		std::vector<S3DCoordinate> listPxydUT;
		std::vector<SCoordinate> listRaDEc = getPointRaDec(tel_itrf, listAzEl_Real, listJD, listVectorAzEl_Real, listPxydUT );
		for (int i = 0; i < listRaDEc.size(); i++)
		{
			S3DCoordinate pt = SphericalToXYZ_Rad(listRaDEc[i].x, listRaDEc[i].y, 1);
			listVectorRaDec.push_back(pt);
		}

		// points of real view
		std::vector<S3DCoordinate> Ref_Real;
		for (int i = 0; i < listAzEl_Real.size(); i++)
		{
			SCoordinate pt = listAzEl_Real[i];
			S3DCoordinate pt3 = SphericalToXYZ_Deg(pt.x, pt.y, 1.0);
			Ref_Real.push_back(pt3);
		}
		
		std::vector<S3DCoordinate> Ref_Tel;
		for (int i = 0; i < listAzEl_Tel.size(); i++)
		{
			SCoordinate pt = listAzEl_Tel[i];
			S3DCoordinate pt3 = SphericalToXYZ_Deg(pt.x, pt.y, 1.0);
			Ref_Tel.push_back(pt3);
		}

		S3DMatrix Rmat;
		estimateRT(Ref_Real, Ref_Tel, Rmat);
				
		SCoordinate azel_mount = SCoordinate(40, 50);
		S3DCoordinate viewInMount = Rmat.inverse()*SphericalToXYZ_Deg(azel_mount.x, azel_mount.y, 1.0);
		SCoordinate azel_view;
		ConvertXYZtoRADEC(viewInMount, azel_view);
		printf("azel_mount [az, el]");
		azel_mount.print();
		printf("azel_view [az, el]");
 		azel_view = azel_view*const_.RG;
		azel_view.print();
	}
};