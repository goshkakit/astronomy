//
// main.cpp 
// This sample code demonstrates how to use the C++ classes in order
// to determine satellite position and look angles.
//
// Copyright (c) 2003-2013 Michael F. Henry
//
// 01/2013
//
#include "Norad_Test.h"

#include "std_add.h"

#include <stdio.h>

#include "coreLib.h"
#include "cOrbit.h"

// This namespace contains all the OrbitTools classes; be sure
// to use this namespace.
using namespace Zeptomoby::OrbitTools;

/////////////////////////////////////////////////////////////////////////////
// Test routine to output position and velocity information
void PrintPosVel(const cTle &tle)
{
	cOrbit       orbit(tle);
	vector<cEci> Pos;

	// Calculate position, velocity
	// mpe = "minutes past epoch"
	for (int mpe = 0; mpe <= (360 * 4); mpe += 360)
	{
		// Get the position of the satellite at time "mpe"
		cEciTime eci = orbit.GetPosition(mpe);

		// Push the coordinates object onto the end of the vector.
		Pos.push_back(eci);
	}

	// Print TLE data
	printf("%s\n",   tle.Name().c_str());
	printf("%s\n",   tle.Line1().c_str());
	printf("%s\n\n", tle.Line2().c_str());

	// Header
	printf("  TSINCE            X                Y                Z\n\n");

	// Iterate over each of the ECI position objects pushed onto the
	// position vector, above, printing the ECI position information
	// as we go.
	for (unsigned int i = 0; i < Pos.size(); i++)
	{
		printf("%8d.00  %16.8f %16.8f %16.8f\n",
			i * 360,
			Pos[i].Position().m_x,
			Pos[i].Position().m_y,
			Pos[i].Position().m_z);
	}

	printf("\n                    XDOT             YDOT             ZDOT\n\n");

	// Iterate over each of the ECI position objects in the position
	// vector again, but this time print the velocity information.
	for (unsigned int i = 0; i < Pos.size(); i++)
	{
		printf("             %16.8f %16.8f %16.8f\n",
			Pos[i].Velocity().m_x,
			Pos[i].Velocity().m_y,
			Pos[i].Velocity().m_z);
	}
}
//==============================================================================//
// test out
//SGP4 Test
//1 88888U          80275.98708465  .00073094  13844-3  66816-4 0     8
//2 88888  72.8435 115.9689 0086731  52.6988 110.5714 16.05824518   105
//
//  TSINCE            X                Y                Z
//
//       0.00     2328.97070437   -5995.22083333    1719.97065611
//     360.00     2456.10787343   -6071.93868176    1222.89554078
//     720.00     2567.56296224   -6112.50380175     713.96182588
//    1080.00     2663.09017513   -6115.48274470     196.39907463
//    1440.00     2742.55440194   -6079.66984186    -326.39149750
//
//                    XDOT             YDOT             ZDOT
//
//                   2.91207225      -0.98341531      -7.09081697
//                   2.67938906      -0.44828838      -7.22879265
//                   2.44024485       0.09811117      -7.31995951
//                   2.19612076       0.65241695      -7.36282421
//                   1.94849696       1.21107421      -7.35619305
//
//SDP4 Test
//1 11801U          80230.29629788  .01431103  00000-0  14311-1       8
//2 11801  46.7916 230.4354 7318036  47.4722  10.4117  2.28537848     6
//
//  TSINCE            X                Y                Z
//
//       0.00     7473.37213351     428.95462549    5828.74786677
//     360.00    -3305.22417985   32410.86360001  -24697.17732308
//     720.00    14271.28695394   24110.46300337   -4725.76750899
//    1080.00    -9990.05752318   22717.36123643  -23616.89356981
//    1440.00     9787.87233694   33753.34427732  -15030.80628319
//
//                    XDOT             YDOT             ZDOT
//
//                   5.10715289       6.44468289      -0.18613182
//                  -1.30113547      -1.15131513      -0.28333528
//                  -0.32050442       2.67984097      -2.08405301
//                  -1.01667268      -2.29026701       0.72892308
//                  -1.09425038       0.92358954      -1.52230979
//
//Example output:
//AZ: 100.192  EL: 12.852
//==============================================================================//
int TestNorad( )
{
	// Test SGP4
	string str1 = "SGP4 Test";
	string str2 = "1 88888U          80275.98708465  .00073094  13844-3  66816-4 0     8";
	string str3 = "2 88888  72.8435 115.9689 0086731  52.6988 110.5714 16.05824518   105";

	cTle tleSGP4(str1, str2, str3);

	PrintPosVel(tleSGP4);

	printf("\n");

	// Test SDP4
	str1 = "SDP4 Test";
	str2 = "1 11801U          80230.29629788  .01431103  00000-0  14311-1       8";
	str3 = "2 11801  46.7916 230.4354 7318036  47.4722  10.4117  2.28537848     6";

	cTle tleSDP4(str1, str2, str3);

	PrintPosVel(tleSDP4);

	printf("\nExample output:\n");

	// Example: Define a location on the earth, then determine the look-angle
	// to the SDP4 satellite defined above.

	// Create an Orbit object using the SDP4 TLE object.
	cOrbit orbitSDP4(tleSDP4);

	// Get the location of the satellite from the Orbit object. The 
	// earth-centered inertial information is placed into eciSDP4.
	// Here we ask for the location of the satellite 90 minutes after
	// the TLE epoch.
	cEciTime eciSDP4 = orbitSDP4.GetPosition(90.0);

	// Now create a site object. Site objects represent a location on the 
	// surface of the earth. Here we arbitrarily select a point on the
	// equator.
	cSite siteEquator(0.0, -100.0, 0); // 0.00 N, 100.00 W, 0 km altitude

	// Now get the "look angle" from the site to the satellite. 
	// Note that the ECI object "eciSDP4" contains a time associated
	// with the coordinates it contains; this is the time at which
	// the look angle is valid.
	cTopo topoLook = siteEquator.GetLookAngle(eciSDP4);

	// Print out the results.
	printf("AZ: %.3f  EL: %.3f\n", 
		topoLook.AzimuthDeg(),
		topoLook.ElevationDeg());

	return 0;
}

//==============================================================================//
// прогноз строчки NORAD на необходимый момент времени
//==============================================================================//
void NoradPredict()
{
	//GLONASS-102	
	//R15
	//0606201
	//29670
	//19140
	//65
	//COSMOS 2425 (716)       
	//1 29670U 06062A   13045.46180505 -.00000066  00000-0  10000-3 0  7232
	//2 29670  65.8927 357.9765 0023418 342.6485 161.3106  2.13105197 47788

	string str1 = "COSMOS 2425";
	string str2 = "1 29670U 06062A   13045.46180505 -.00000066  00000-0  10000-3 0  7232";
	string str3 = "2 29670  65.8927 357.9765 0023418 342.6485 161.3106  2.13105197 47788";

	cTle tle(str1, str2, str3);
	cOrbit orbit(tle);
	cJulian jd = orbit.Epoch();

	printf( "Date = %f\n", jd.Date() );

	cTle tleSGP4(str1, str2, str3);
	PrintPosVel(tleSGP4);

	double DayFull = 45.46180505;
	int Day = (int)DayFull;
	int hour = (( DayFull - (double)Day )*86400.0)/3600;
	int minute = ( (DayFull - (double)Day)*86400.0 - ((double)hour)*3600.0 )/60.0;
	double sec = ( (DayFull - (double)Day)*86400.0 - ((double)hour)*3600.0  - (double)minute*60.0);
	printf("day = %d\nhour = %d\nminute = %d\nsec = %f\n", Day, hour, minute, sec );

	// sec 39899.95632
	double fullsec = (double)sec + 60.0*minute + 3600.0*hour;
	printf("Day SEC = %f (~39899.95632)\n", fullsec );

	// 2456337.961805
	double secjd = (jd.Date() - 2456337.5)*86400;

	printf("jd.Date() = %f\n", secjd );

	// out
	//day = 45
	//hour = 11
	//minute = 4
	//sec = 59.956320
	//Day SEC = 39899.956320 (~39899.95632)
	//jd.Date() = 39899.956302

	// Time
	// 2013 02 14 11 04 59.956320

	//POSITION
	// time start 2013  2 14  0  0  0
	//40500.00000  0 -19402942.045 -11660867.925  11875536.143
	// time pos 2012 02 14 11 15 00

	double timepos = 40500.00000;

	double H = (int)( timepos/3600.0 );
	double M = (int)(( timepos - H*3600.0 )/60.0);
	double S = timepos - H*3600.0 - M*60.0;
	printf( "H = %f,  D = %f, S = %f\n", H, M, S );
	// H = 11.000000,  D = 15.000000, S = 0.000000

	// время предсказания вперед в минутах
	double minpredict = (timepos - secjd)/60.0;

	printf( "minpredict = %f\n", minpredict );

	
	printf("  TSINCE            X                Y                Z\n");
	// Iterate over each of the ECI position objects pushed onto the
	// position vector, above, printing the ECI position information
	// as we go.
	double it = -35.0/60.0;
	//for( double it = -10; it < 10; it += 0.5 )
	{
		cEciTime eci = orbit.GetPosition( minpredict+it );
		printf("%f  %16.8f %16.8f %16.8f\n", minpredict+it, eci.Position().m_x, eci.Position().m_y, eci.Position().m_z );
	}

};
//==============================================================================//
