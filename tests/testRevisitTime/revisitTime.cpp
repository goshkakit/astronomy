#include <stdio.h>

#include "Norad\coreLib.h"
#include "Norad\cOrbit.h"

#include "common\DataConverter.h"
#include "common\mytypes.h"


#include "InfluenceForce\InfluenceForce.h"

#include "Norad\sgp4ext.h"
#include "Norad\sgp4unit.h"
#include "Norad\sgp4io.h"

int agiNoradTest()
{
	// ----------------------------  TEME   --------------------------------
	//Force::InfluenceForce *IForce = new Force::InfluenceForce();
	//IForce->Init_CPU();
	//DataConverter Dconv;


	char str[2];
	char infilename[15];
	double ro[3];
	double vo[3];
	char typerun, typeinput, opsmode;
	gravconsttype  whichconst;
	int whichcon;
	FILE *infile, *outfile, *outfilee;

	// ----------------------------  locals  -------------------------------
	double p, a, ecc, incl, node, argp, nu, m, arglat, truelon, lonper;
	double sec, jd, rad, tsince, startmfe, stopmfe, deltamin;
	double tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2;
	int  year; int mon; int day; int hr; int min;
	char longstr1[130];
	typedef char str3[4];
	str3 monstr[13];
	char outname[64];
	char longstr2[130];
	elsetrec satrec;

	rad = 180.0 / pi;
	// ------------------------  implementation   --------------------------
	strcpy(monstr[1], "Jan");
	strcpy(monstr[2], "Feb");
	strcpy(monstr[3], "Mar");
	strcpy(monstr[4], "Apr");
	strcpy(monstr[5], "May");
	strcpy(monstr[6], "Jun");
	strcpy(monstr[7], "Jul");
	strcpy(monstr[8], "Aug");
	strcpy(monstr[9], "Sep");
	strcpy(monstr[10], "Oct");
	strcpy(monstr[11], "Nov");
	strcpy(monstr[12], "Dec");

	printf("%s\n", SGP4Version);

	//opsmode = 'a' best understanding of how afspc code works
	//opsmode = 'i' imporved sgp4 resulting in smoother behavior
	printf("input operation mode a, i \n\n");
	opsmode = getchar();
	fflush(stdin);

	//typerun = 'c' compare 1 year of full satcat data
	//typerun = 'v' verification run, requires modified elm file with
	//              start, stop, and delta times
	//typerun = 'm' maunual operation- either mfe, epoch, or dayof yr also
	printf("input type of run c, v, m \n\n");
	typerun = getchar();
	fflush(stdin);

	//typeinput = 'm' input start stop mfe
	//typeinput = 'e' input start stop ymd hms
	//typeinput = 'd' input start stop yr dayofyr
	if ((typerun != 'v') && (typerun != 'c'))
	{
		printf("input mfe, epoch (YMDHMS), or dayofyr approach, m,e,d \n\n");
		typeinput = getchar();
	}
	else
		typeinput = 'e';

	printf("input which constants 721 72 84 \n");
	scanf("%i", &whichcon);
	if (whichcon == 721) whichconst = wgs72old;
	if (whichcon == 72) whichconst = wgs72;
	if (whichcon == 84) whichconst = wgs84;

	getgravconst(whichconst, tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2);

	// ---------------- setup files for operation ------------------
	// input 2-line element set file
	printf("input elset filename: \n");
	scanf("%s", infilename);
	infile = fopen(infilename, "r");
	if (infile == NULL)
	{
		printf("Failed to open file: %s\n", infilename);
		return 1;
	}

	if (typerun == 'c')
		outfile = fopen("tcppall.out", "w");
	else
	{
		if (typerun == 'v')
			outfile = fopen("tcppver.out", "w");
		else
			outfile = fopen("tcpp.out", "w");
	}

	//        dbgfile = fopen("sgp4test.dbg", "w");
	//        fprintf(dbgfile,"this is the debug output\n\n" );

	// ----------------- test simple propagation -------------------
	while (feof(infile) == 0)
	{
		do
		{
			fgets(longstr1, 130, infile);
			strncpy(str, &longstr1[0], 1);
			str[1] = '\0';
		} while ((strcmp(str, "#") == 0) && (feof(infile) == 0));

		if (feof(infile) == 0)
		{
			fgets(longstr2, 130, infile);
			// convert the char string to sgp4 elements
			// includes initialization of sgp4
			twoline2rv(longstr1, longstr2, typerun, typeinput, opsmode, whichconst,
				startmfe, stopmfe, deltamin, satrec);
			fprintf(outfile, "%ld xx\n", satrec.satnum);
			printf(" %ld\n", satrec.satnum);
			// call the propagator to get the initial state vector value
			sgp4(whichconst, satrec, 0.0, ro, vo);

			// generate .e files for stk
			jd = satrec.jdsatepoch;
			strncpy(outname, &longstr1[2], 5);
			outname[5] = '.';
			outname[6] = 'e';
			outname[7] = '\0';
			invjday(jd, year, mon, day, hr, min, sec);
			outfilee = fopen(outname, "w");
			fprintf(outfilee, "stk.v.4.3 \n"); // must use 4.3...
			fprintf(outfilee, "\n");
			fprintf(outfilee, "BEGIN Ephemeris \n");
			fprintf(outfilee, " \n");
			fprintf(outfilee, "NumberOfEphemerisPoints		146 \n");
			fprintf(outfilee, "ScenarioEpoch	  %3i %3s%5i%3i:%2i:%12.9f \n", day, monstr[mon],
				year, hr, min, sec);
			fprintf(outfilee, "InterpolationMethod		Lagrange \n");
			fprintf(outfilee, "InterpolationOrder		5 \n");
			fprintf(outfilee, "CentralBody				Earth \n");
			fprintf(outfilee, "CoordinateSystem			TEME \n");
			fprintf(outfilee, "CoordinateSystemEpoch	%3i %3s%5i%3i:%2i:%12.9f \n", day,
				monstr[mon], year, hr, min, sec);
			fprintf(outfilee, "DistanceUnit			Kilometers \n");
			fprintf(outfilee, " \n");
			fprintf(outfilee, "EphemerisTimePosVel \n");
			fprintf(outfilee, " \n");
			fprintf(outfilee, " %16.8f %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f\n",
				satrec.t, ro[0], ro[1], ro[2], vo[0], vo[1], vo[2]);

			fprintf(outfile, " %16.8f %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f\n",
				satrec.t, ro[0], ro[1], ro[2], vo[0], vo[1], vo[2]);

			tsince = startmfe;
			// check so the first value isn't written twice
			if (fabs(tsince) > 1.0e-8)
				tsince = tsince - deltamin;

			// ----------------- loop to perform the propagation ----------------
			bool satVisible = false;
			while ((tsince < stopmfe) && (satrec.error == 0))
			{
				tsince = tsince + deltamin;

				if (tsince > stopmfe)
					tsince = stopmfe;

				sgp4(whichconst, satrec, tsince, ro, vo);

				if (satrec.error > 0)
					printf("# *** error: t:= %f *** code = %3d\n",
					satrec.t, satrec.error);


				if (satrec.error == 0)
				{
					if ((typerun != 'v') && (typerun != 'c'))
					{
						jd = satrec.jdsatepoch + tsince / 1440.0;
						invjday(jd, year, mon, day, hr, min, sec);

						/*
						// convert
						double zenit = 0;
						bool writeLine = false;
						if ( 1)
						{
							// MDB
							double jdm = satrec.jdsatepoch + 0.125 - 0.5;
							double date1 = Dconv.JDtoYYYYMMDD(jdm);
							double time1 = Dconv.SECtoHHMMSS(date1, jdm);

							// установка времени
							double int1, ajd1, delt1;
							IForce->set_time(date1, time1, &ajd1, &delt1, &int1);

							// матрица перехода из ICRF в TEME
							double A_Teme[9];
							IForce->GetTemeMatrix(int1, A_Teme, ajd1, delt1);
							// обратная матрица перехода из TEME в ICRF
							double invA_Teme[9];
							IForce->transpose(A_Teme, invA_Teme);

							double Ptle[3];
							Ptle[0] = ro[0];
							Ptle[1] = ro[1];
							Ptle[2] = ro[2];
							double PICRF[3];
							// перевод вектора состояния в ICRF 
							IForce->matVecMul(invA_Teme, Ptle, PICRF);

							jdm = jd + 0.125 - 0.5;
							date1 = Dconv.JDtoYYYYMMDD(jdm);
							time1 = Dconv.SECtoHHMMSS(date1, jdm);

							// матрица перевода в земную систему
							double Arot[9];
							IForce->iers_update_matrix(int1, Arot, ajd1, delt1);

							// матрица перехода из земной в нормальную систему
							double invArot[9];
							IForce->transpose(Arot, invArot);

							double Tpos[3];
							double lon = 37.433736 / rad;
							double lat = 55.858237 / rad;
							double h = 6378.135 + 0.186820250686651;
							Tpos[0] = h*cos(lat)*cos(lon);
							Tpos[1] = h*cos(lat)*sin(lon);
							Tpos[2] = h*sin(lat);

							// перевод в ICRF координат телескопа
							double tel_icrf[3];
							IForce->matVecMul(invArot, Tpos, tel_icrf);

							// вектор телескопа
							S3DCoordinate Rtel;
							Rtel.x = tel_icrf[0];
							Rtel.y = tel_icrf[1];
							Rtel.z = tel_icrf[2];

							// вектор наблюдения
							S3DCoordinate Rn;
							Rn.x = PICRF[0];
							Rn.y = PICRF[1];
							Rn.z = PICRF[2];
							S3DCoordinate Rmeg = Rn - Rtel;

							double cosa = (Rtel*Rmeg) / (Rtel.norm()*Rmeg.norm());
							zenit = acos(cosa)*rad;
						}

						if (zenit >= 80.0 && !satVisible)
						{
							satVisible = true;
							writeLine = true;
						}
						if (zenit < 80.0 && satVisible)
						{
							satVisible = false;
							writeLine = true;
						}
						*/
						if (1)
						{
							fprintf(outfile, " %16.8f %16.8f %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f %5i%3i%3i %2i:%2i:%9.6f\n",
								tsince, jd, ro[0], ro[1], ro[2], vo[0], vo[1], vo[2], year, mon, day, hr, min, sec);
							//                            fprintf(outfile, " %16.8f %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f\n",
							//                                           tsince,ro[0],ro[1],ro[2],vo[0],vo[1],vo[2]);
						}
					}
					else
					{
						jd = satrec.jdsatepoch + tsince / 1440.0;
						invjday(jd, year, mon, day, hr, min, sec);

						fprintf(outfilee, " %16.6f %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f \n",
							tsince*60.0, ro[0], ro[1], ro[2], vo[0], vo[1], vo[2]);

						fprintf(outfile, " %16.8f %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f",
							tsince, ro[0], ro[1], ro[2], vo[0], vo[1], vo[2]);

						rv2coe(ro, vo, mu, p, a, ecc, incl, node, argp, nu, m, arglat, truelon, lonper);
						fprintf(outfile, " %14.6f %8.6f %10.5f %10.5f %10.5f %10.5f %10.5f %5i%3i%3i %2i:%2i:%9.6f\n",
							a, ecc, incl*rad, node*rad, argp*rad, nu*rad,
							m*rad, year, mon, day, hr, min, sec);
					}
				} // if satrec.error == 0

			} // while propagating the orbit

			fprintf(outfilee, " END Ephemeris \n");
			fclose(outfilee);

		} // if not eof

	} // while through the input file


	return 0;
}

void main()
{
	//agiNoradTest();
	//return;
	// 1 90732U 00000A   18206.00000000  .00000000  00000-0  58106+0 0    11
	// 2 90732  19.0914   8.4920 2172097 149.9634 343.3148  1.01981561000013
	//std::string str1 = "90732";
	//std::string str2 = "1 90732U 00000A   18206.00000000  .00000000  00000-0  58106+0 0    11";
	//std::string str3 = "2 90732  19.0914   8.4920 2172097 149.9634 343.3148  1.01981561000013";

	std::string str1 = "28492";
	std::string str2 = "1 28492U 00000A   18271.00000000  .00000000  00000-0  21133-4 0      ";
	std::string str3 = "2 28492  98.1014 205.9586 0000999 137.3407 338.9473 14.63855357000012";

	// init rotation
	//Force::InfluenceForce *FI = new Force::InfluenceForce();
	//FI->Init_CPU();

	// init norad
	cTle tleSGP4(str1, str2, str3);

	// vector
	cOrbit *orbit = new cOrbit(tleSGP4);

	// Now create a site object. Site objects represent a location on the 
	// surface of the earth. Here we arbitrarily select a point on the
	// equator.
	//cSite siteEquator(55.858237, 37.433736, 0.186820); // 0.00 N, 100.00 W, 0 km altitude
	cSite siteEquator(-29.8999999960317, -71.2500000018608, 0);

	double SAT_elv = 10.0;

	//22186.0000000

	double JDtle = orbit->Epoch().Date();
	double jdoffset = 2436203.5;
	double jdm = JDtle - jdoffset;
	printf("JDtle = %f\t jdm = %f\n", JDtle, jdm );

	bool SatVisible = false;
	// Calculate position, velocity
	// mpe = "minutes past epoch
	double stepSec = 0.01;
	double stepMin = stepSec/60.0;
	double dayMinutes = 2.0*24 * 60;
	double JDstart = 0;

	FILE *fre = fopen("orbit-coord-az.txt", "w");

	for (double mpe = 0; mpe <= dayMinutes; mpe += stepMin)
	{
		// Get the position of the satellite at time "mpe"
		cEciTime eci = orbit->GetPosition(mpe + 1440.0);

		// Now get the "look angle" from the site to the satellite. 
		// Note that the ECI object "eciSDP4" contains a time associated
		// with the coordinates it contains; this is the time at which
		// the look angle is valid.
		cTopo topoLook = siteEquator.GetLookAngle(eci);
		//double currenJd = jdm + mpe / 24.0 / 60.0;
		double currenJd = JDtle + (mpe + 1440.0) / 24.0 / 60.0;

		//if (topoLook.ElevationDeg() >= 0)
		//{
		//	fprintf(fre, "%f %.7f %f %f %f %f %f \n", mpe, JDtle + (mpe + 1440.0) / 24.0 / 60.0, eci.Position().m_x, eci.Position().m_y, eci.Position().m_z, topoLook.AzimuthDeg(), topoLook.ElevationDeg());

		//}

		if (1)
		{

			if (topoLook.ElevationDeg() >= SAT_elv && !SatVisible)
			{
				SatVisible = true;
				// Print out the results.
				printf("AZ: %.3f\tEL: %.3f\tJDs: %.7f\t", topoLook.AzimuthDeg(), topoLook.ElevationDeg(), currenJd);
				JDstart = currenJd;
			}
			if (topoLook.ElevationDeg() < SAT_elv && SatVisible)
			{
				SatVisible = false;
				double duration = (currenJd - JDstart)*86400.0;
				// Print out the results.
				printf("JDe: %.7f\tDir: %f\tAZ: %.3f\tEL: %.3f\n", currenJd, duration, topoLook.AzimuthDeg(), topoLook.ElevationDeg());

			}
		}
	}

	fclose(fre);
	
}