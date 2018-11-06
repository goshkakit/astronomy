/*
  Naval Observatory Vector Astrometry Software (NOVAS)
  C Edition, Version 3.1
 
  example.c: Examples of NOVAS calculations 
 
  U. S. Naval Observatory
  Astronomical Applications Dept.
  Washington, DC 
  http://www.usno.navy.mil/USNO/astronomical-applications
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "eph_manager.h" /* remove this line for use with solsys version 2 */
#include "novas.h"

int main()
{
	//2008 April 24, 10:36 : 18.0 UTC
	//2007 11 2     9 h 22 min 12 s   UTC

	// current http://maia.usno.navy.mil/ser7/mark3.out

	//const short int year = 2008;
	//const short int month = 4;
	//const short int day = 24;
	//const int h = 10;
	//const int m = 36;
	//const double s = 18;

	const short int year = 2007;
	const short int month = 11;
	const short int day = 2;
	const int h = 9;
	const int m = 22;
	const double s = 12;

	const short int accuracy = 0;
	double jd_tt, jd_utc, hour;
	double jd_ut1, delta_t;

	hour = (double)h + (double)m / 60.0 + s / 3600.0;
	jd_utc = julian_date(year, month, day, hour);
	printf("jd_utc = %lf, hour = %lf\n", jd_utc, hour);

	//1999 JAN  1 = JD 2451179.5  TAI - UTC = 32.0       S + (MJD - 41317.) X 0.0      S
	//2006 JAN  1 = JD 2453736.5  TAI - UTC = 33.0       S + (MJD - 41317.) X 0.0      S
	//2009 JAN  1 = JD 2454832.5  TAI - UTC = 34.0       S + (MJD - 41317.) X 0.0      S
	//2012 JUL  1 = JD 2456109.5  TAI - UTC = 35.0       S + (MJD - 41317.) X 0.0      S

	const double leap_secs = 33.0;

	//Year, Month, Day, Modified Julian Date
	//PM - x[arcsec], error_PM - x[arcsec], 5
	//PM - y[arcsec], error_PM - y[arcsec], 7
	//UT1 - UTC[seconds], error_UT1 - UTC[seconds], 9
	//LOD[milliseconds], error_LOD[milliseconds],
	//dX[milliarcsec], error_dX[milliarcsec], 13
	//dY[milliarcsec], error_dY[milliarcsec]  15
	//(from Bulletin A)
	//
	//PM - x[arcsec], PM - y[arcsec], UT1 - UTC[seconds], 17,18,19
	//dX[milliarcsec], dY[milliarcsec]
	//(from Bulletin B)

	//8 423 54579.00 I -0.005684 0.000043  0.527458 0.000025  I -0.3869274 0.0000059  0.7107 0.0050  I  0.173    0.046 -0.152    0.340		-.005760   .527510 -.3869340     0.298	-0.317
	//8 424 54580.00 I -0.003111 0.000043  0.528424 0.000024  I -0.3875950 0.0000058  0.6313 0.0038  I  0.187    0.046 -0.170    0.340		-.003100   .528470 -.3875920     0.338	-0.261
	//8 425 54581.00 I  0.000019 0.000043  0.529163 0.000025  I -0.3881989 0.0000048  0.5788 0.0038  I  0.213    0.046 -0.187    0.340		-.000050   .529180 -.3881940     0.343	-0.171
	//const double ut1_utc = 0;
	//const double x_pole = -0.002;
	//const double y_pole = +0.529;
	//double mjd = jd_utc - 2400000.5;
	//double jd1 = 54580.00;
	//double jd2 = 54581.00;
	//double x1 = -0.0031;
	//double x2 = 0.00005;
	//double y1 = 0.528424;
	//double y2 = 0.529163;
	//double ut1 = -0.3875920;
	//double ut2 = -0.3881940;

	//711 1 54405.00 I  0.044863 0.000038  0.191004 0.000042  I - 0.2132016 0.0000037  0.6134 0.0027  I - 64.809     .328 - 5.076     .340   .044820   .191060 - .2132210 - 64.000 - 5.000
	//711 2 54406.00 I  0.042326 0.000038  0.191057 0.000043  I - 0.2138878 0.0000042  0.7657 0.0043  I - 64.490     .328 - 5.038     .340   .042290   .191110 - .2139220 - 63.600 - 4.900
	//711 3 54407.00 I  0.040333 0.000038  0.191442 0.000042  I - 0.2147377 0.0000078  0.9320 0.0032  I - 64.313     .473 - 4.674     .340   .040330   .191490 - .2147400 - 63.500 - 4.600

	// A2000
	//711 1 54405.00 I  0.044850 0.000038  0.191027 0.000042  I -0.2131998 0.0000037  0.6134 0.0027  I -0.087    0.130 -0.114    0.340		.044820   .191060 -.2132210     0.235	-0.039
	//711 2 54406.00 I  0.042313 0.000038  0.191080 0.000043  I -0.2138860 0.0000042  0.7657 0.0043  I -0.083    0.130 -0.112    0.340		.042290   .191110 -.2139220     0.272    0.026
	//711 3 54407.00 I  0.040320 0.000038  0.191465 0.000042  I -0.2147359 0.0000078  0.9320 0.0032  I -0.093    0.188 -0.058    0.340		.040330   .191490 -.2147400     0.230    0.016
	double mjd = jd_utc - 2400000.5 + 3.0/24.0;
	double jd1 = 54406;
	double jd2 = 54407;

	double x1 = 0.042313;
	double x2 = 0.040333;

	double y1 = 0.191057;
	double y2 = 0.191442;

	double ut1 = -0.2138878;
	double ut2 = -0.2147377;

	//double x1 = 0.042313;
	//double x2 = 0.040320;

	//double y1 = 0.191080;
	//double y2 = 0.191465;

	//double ut1 = -0.2138860;
	//double ut2 = -0.2147359;

	double x_pole = x1 + (mjd - jd1)*(x2 - x1) / (jd2 - jd1);
	double y_pole = y1 + (mjd - jd1)*(y2 - y1) / (jd2 - jd1);
	double ut1_utc = ut1 + (mjd - jd1)*(ut2 - ut1) / (jd2 - jd1);
	printf("mx = %lf, my = %lf, u = %lf, mjd = %lf\n", x_pole, y_pole, ut1_utc, mjd);



	// ∆AT is an integer representing the total count of leap seconds in UTC
	// TT = UTC + ∆AT + 32.184s
	// TT = TAI + 32.184s
	// TAI is International Atomic Time
	jd_tt = jd_utc + ((double)leap_secs + 32.184) / 86400.0;
	// UT1 = UTC + (UT1–UTC)
	jd_ut1 = jd_utc + ut1_utc / 86400.0;
	// UT1 = TT – ∆T
	// ∆T = TT - UT1 = UTC + ∆AT + 32.184s - UTC - (UT1–UTC)
	// ∆T = ∆AT + 32.184s - (UT1–UTC)
	delta_t = 32.184 + leap_secs - ut1_utc;

	short int error = 0;
	short int de_num = 0;
	double jd_beg, jd_end;
	if ((error = ephem_open("../data/eph/lnx1900.405", &jd_beg, &jd_end, &de_num)) != 0)
	{
		if (error == 1)
			printf("JPL ephemeris file not found.\n");
		else
			printf("Error reading JPL ephemeris file header.\n");
		return (error);
	}
	else
	{
		printf("JPL ephemeris DE%d open. Start JD = %10.2f  End JD = %10.2f\n",
			de_num, jd_beg, jd_end);
		printf("\n");
	}

	double vter[3];
	double vcel[3];

	vter[0] = 0;// 6370000.0;
	vter[1] = 6370000.0;
	vter[2] = 0;// 6370000.0;

	if ((error = ter2cel(jd_ut1, 0.0, delta_t, 1, accuracy, 0, x_pole, y_pole, vter, vcel)) != 0)
	{
		printf("Error %d from ter2cel.", error);
		return (error);
	}

	printf("x,y,z: %lf %lf %lf\n", vcel[0], vcel[1], vcel[2]);

	double Arot[9] = { -0.999547001090, 0.030086434974, 0.000773978848,
		-0.030086414120, -0.999547300630, 0.000038576390,
		0.000774789095, 0.000015272667, 0.999999699734 };

	double Arot2[9] = { -9.99547001036216E-01,      3.00864367479260E-02,      7.73978889005599E-04,
		- 3.00864158935205E-02, - 9.99547300576900E-01 ,     3.85760783290874E-05,
		7.74789125949687E-04,      1.52723526581853E-05,      9.99999699734238E-01 };

	double vcel2[3];
	for (int j = 0; j < 3; j++)
	{
		vcel2[j] = 0;
		for (int i = 0; i < 3; i++)
		{
			vcel2[j] += Arot2[j * 3 + i] * vter[i];
		}
	}

	printf("x,y,z: %lf %lf %lf\n", vcel2[0], vcel2[1], vcel2[2]);
	printf("dx,dy,dz: %lf %lf %lf\n", vcel[0] - vcel2[0], vcel[1] - vcel2[1], vcel[2] - vcel2[2]);
}

int main2 (void)
{

/*
  NOVAS 3.1 Example Calculations
 
  See Chapter 3 of User's Guide for explanation.
 
  Written for use with solsys version 1.
 
  To adapt for use with solsys version 2, see comments throughout file. 
  Assumes JPL ephemeris file "JPLEPH" located in same directory as
  application.
 */

   const short int year = 2008;
   const short int month = 4;
   const short int day = 24;
   const short int leap_secs = 33;
   const short int accuracy = 0;
   short int error = 0;
   short int de_num = 0;
   
   const double hour = 10.605;
   const double ut1_utc = -0.387845;
   
   const double latitude = 42.0;
   const double longitude = -70;
   const double height = 0.0;
   const double temperature = 10.0;
   const double pressure = 1010.0;
      
   const double x_pole = -0.002;
   const double y_pole = +0.529;
   
   double jd_beg, jd_end, jd_utc, jd_tt, jd_ut1, jd_tdb, delta_t, ra, 
      dec, dis, rat, dect, dist, zd, az, rar, decr, gast, last, theta, 
      jd[2], pos[3], vel[3], pose[3], elon, elat, r, lon_rad, lat_rad, 
      sin_lon, cos_lon, sin_lat, cos_lat, vter[3], vcel[3];
   
   on_surface geo_loc;
   
   observer obs_loc;
   
   cat_entry star, dummy_star;
   
   object moon, mars;
   
   sky_pos t_place;
   
/*
   Make structures of type 'on_surface' and 'observer-on-surface' containing 
   the observer's position and weather (latitude, longitude, height, 
   temperature, and atmospheric pressure).
*/
   
   make_on_surface (latitude,longitude,height,temperature,pressure, &geo_loc);
   make_observer_on_surface (latitude,longitude,height,temperature,pressure,
      &obs_loc);

/*
   Make a structure of type 'cat_entry' containing the ICRS position 
   and motion of star FK6 1307.
*/

   make_cat_entry ("GMB 1830","FK6",1307,11.88299133,37.71867646, 
      4003.27,-5815.07,109.21,-98.8, &star);
      
/*
   Make structures of type 'object' for the Moon and Mars.
*/

   make_cat_entry ("DUMMY","xxx",0,0.0,0.0,0.0,0.0,0.0,0.0, 
      &dummy_star);
  
   if ((error = make_object (0,11,"Moon",&dummy_star, &moon)) != 0)
   {
      printf ("Error %d from make_object (Moon)\n", error);
      return (error);
   }

   if ((error = make_object (0,4,"Mars",&dummy_star, &mars)) != 0)
   {
      printf ("Error %d from make_object (Mars)\n", error);
      return (error);
   }
   
/*
   Open the JPL binary ephemeris file, here named "JPLEPH".
   Remove this block for use with solsys version 2.
*/

   //if ((error = ephem_open("JPLEPH", &jd_beg, &jd_end, &de_num)) != 0) 
   if ((error = ephem_open("../data/eph/lnx1900.405", &jd_beg, &jd_end, &de_num)) != 0)
   {
	   if (error == 1)
		   printf("JPL ephemeris file not found.\n");
	   else
		   printf("Error reading JPL ephemeris file header.\n");
	   return (error);
   }
   else
   {
	   printf("JPL ephemeris DE%d open. Start JD = %10.2f  End JD = %10.2f\n",
		   de_num, jd_beg, jd_end);
	   printf("\n");
   }

/*
  Uncomment block below for use with solsys version 2
  Prints alternate header
*/
/*
   printf ("Using solsys version 2, no description of JPL ephemeris available\n");
   printf ("\n");
*/

/*
   Write banner.
*/

   printf ("NOVAS Sample Calculations\n");
   printf ("-------------------------\n");
   printf ("\n");
   
/*
   Write assumed longitude, latitude, height (ITRS = WGS-84).
*/

   printf ("Geodetic location:\n");
   printf ("%15.10f        %15.10f        %15.10f\n\n", geo_loc.longitude,
      geo_loc.latitude, geo_loc.height);

/*
   Establish time arguments.
*/

   jd_utc = julian_date (year,month,day,hour);
   jd_tt = jd_utc + ((double) leap_secs + 32.184) / 86400.0;
   jd_ut1 = jd_utc + ut1_utc / 86400.0;
   delta_t = 32.184 + leap_secs - ut1_utc;
   
   jd_tdb = jd_tt;          /* Approximation good to 0.0017 seconds. */
   
   printf ("TT and UT1 Julian Dates and Delta-T:\n");
   printf ("%15.6f        %15.6f        %16.11f\n", jd_tt, jd_ut1, delta_t);
   printf ("\n");
      
/*
   Apparent and topocentric place of star FK6 1307 = GMB 1830.
*/

   if ((error = app_star (jd_tt,&star,accuracy, &ra,&dec)) != 0)
   {
      printf ("Error %d from app_star.\n", error);
      return (error);
   }
      
   if ((error = topo_star (jd_tt,delta_t,&star,&geo_loc, accuracy, 
      &rat,&dect)) != 0)
   {
      printf ("Error %d from topo_star.\n", error);
       return (error);
   }
    
   printf ("FK6 1307 geocentric and topocentric positions:\n");
   printf ("%15.10f        %15.10f\n", ra, dec);
   printf ("%15.10f        %15.10f\n", rat, dect);
   printf ("\n");
     
/*
   Apparent and topocentric place of the Moon.
*/

   if ((error = app_planet (jd_tt,&moon,accuracy, &ra,&dec,&dis)) != 0)
   {
      printf ("Error %d from app_planet.", error);
      return (error);
   }

   if ((error = topo_planet (jd_tt,&moon,delta_t,&geo_loc,accuracy,
      &rat,&dect,&dist)) != 0)
   {
      printf ("Error %d from topo_planet.", error);
      return (error);
   }

   printf ("Moon geocentric and topocentric positions:\n");
   printf ("%15.10f        %15.10f        %15.12f\n", ra, dec, dis);
   printf ("%15.10f        %15.10f        %15.12f\n", rat, dect, dist);

/*
   Topocentric place of the Moon using function 'place'-- should be 
   same as above.
*/

   if ((error = place (jd_tt,&moon,&obs_loc,delta_t,1,accuracy, 
      &t_place)) != 0)
   {
      printf ("Error %d from place.", error);
      return (error);
   }
   
   printf ("%15.10f        %15.10f        %15.12f\n", t_place.ra,t_place.dec,
      t_place.dis);
   printf ("\n");
    
/*
   Position of the Moon in local horizon coordinates.  (Polar motion 
   ignored here.)
*/
   
   equ2hor (jd_ut1,delta_t,accuracy,0.0,0.0,&geo_loc,rat,dect,1,
      &zd,&az,&rar,&decr);

   printf ("Moon zenith distance and azimuth:\n");
   printf ("%15.10f        %15.10f\n", zd, az);
   printf ("\n");

/*
   Greenwich and local apparent sidereal time and Earth Rotation Angle.
*/

   if ((error = sidereal_time (jd_ut1,0.0,delta_t,1,1,accuracy, &gast)) != 0)
   {
      printf ("Error %d from sidereal_time.", error);
      return (error);
   }
   
   last = gast + geo_loc.longitude / 15.0;
   if (last >= 24.0)
      last -= 24.0;
   if (last < 0.0)
      last += 24.0;
      
   theta = era (jd_ut1,0.0);

   printf ("Greenwich and local sidereal time and Earth Rotation Angle:\n");
   printf ("%16.11f        %16.11f        %15.10f\n", gast, last, theta);   
   printf ("\n");

/*      
   Heliocentric position of Mars in BCRS.
   
   TDB ~ TT approximation could lead to error of ~50 m in position of Mars.
*/

   jd[0] = jd_tdb;
   jd[1] = 0.0;
   if ((error = ephemeris (jd,&mars,1,accuracy, pos,vel)) != 0)
   {
      printf ("Error %d from ephemeris (Mars).", error);
      return (error);
   }

   if ((error = equ2ecl_vec (T0,2,accuracy,pos, pose)) != 0)  
   {
      printf ("Error %d from equ2ecl_vec.", error);
      return (error);
   }

   if ((error = vector2radec (pose, &elon,&elat)) != 0)
   {
      printf ("Error %d from vector2radec.", error);
      return (error);
   }
   elon *= 15.0;
   
   r = sqrt (pose[0] * pose[0] + pose[1] * pose[1] + pose[2] * pose[2]);
   
   printf ("Mars heliocentric ecliptic longitude and latitude and "
           "radius vector:\n");
   printf ("%15.10f        %15.10f        %15.12f\n", elon, elat, r);   
   printf ("\n");

/*
   Terrestrial to celestial transformation.
*/

   lon_rad = geo_loc.longitude * DEG2RAD;
   lat_rad = geo_loc.latitude * DEG2RAD;
   sin_lon = sin (lon_rad);
   cos_lon = cos (lon_rad);
   sin_lat = sin (lat_rad);
   cos_lat = cos (lat_rad);

/*      
   Form vector toward local zenith (orthogonal to ellipsoid) in ITRS.
*/

   vter[0] = cos_lat * cos_lon;
   vter[1] = cos_lat * sin_lon;
   vter[2] = sin_lat;

/*      
   Transform vector to GCRS.
*/

   if ((error = ter2cel (jd_ut1,0.0,delta_t,1,accuracy,0,x_pole,y_pole,vter,
      vcel)) != 0)
   {
      printf ("Error %d from ter2cel.", error);
      return (error);
   }
   
   if ((error = vector2radec (vcel, &ra,&dec)) != 0)
   {
      printf ("Error %d from vector2radec.", error);
      return (error);
   }

   printf ("Direction of zenith vector (RA & Dec) in GCRS:\n");
   printf ("%15.10f        %15.10f\n", ra, dec);
   printf ("\n");
   

   ephem_close();  /* remove this line for use with solsys version 2 */
      
   return (0);
}
