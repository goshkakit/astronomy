//==============================================================================//
// Andrianov N.G.
// opbit predict 
// module find Influence Force
// Nutation Earth
//==============================================================================//
#include <math.h>
#include <stdio.h>

#include "InfluenceForce.h"

extern "C"
{
#include "eph_manager.h" /* remove this line for use with solsys version 2 */
#include "novas.h"
}

namespace Force
{
	//==============================================================================//
	// функция тестирования правильности вычисления поворота
	// http://hpiers.obspm.fr/eop-pc/index.php?index=matrice_php&lang=fr
	// out 1
	//tab[i4+1] = 35.000000
	//tdb_utc return 67.184000
	//set_time: Data = 20120922.000000, Time = 202215.000000
	//t = 73.335000,   jd = 2456192.500000 delta -10.732816, CH = 20.370833
	//ITRF to ICRF
	//-1.327546022239e-001     9.911481274985e-001     1.266864067341e-003
	//-9.911489262903e-001     -1.327546814281e-001    -2.173887142407e-005
	//1.466356939667e-004      -1.258536895330e-003    9.999991972911e-001
	//IERS res:
	//-1.327545872600e-001     9.911481297030e-001     1.266707241000e-003
	//-9.911489283060e-001     -1.327546663710e-001    -2.179443100000e-005
	//1.465597880000e-004      -1.258388835000e-003    9.999991974890e-001
	//Delta IERS
	//-1.496387300182e-008     -2.204483262602e-009    1.568263414890e-007
	//2.015730249383e-009      -1.505706484495e-008    5.555957592839e-008
	//7.590596667645e-008      -1.480603299404e-007    -1.978939234704e-01110
	//
	// out 2
	//tab[i4+1] = 35.000000
	//tdb_utc return 67.184000
	//set_time: Data = 20130215.000000, Time = 160853.000000
	//t = 58.133000,   jd = 2456338.500000 delta -10.732816, CH = 16.148056
	//ITRF to ICRF
	//9.550444502895e-001      2.964597735295e-001     1.304089754509e-003
	//-2.964599918975e-001     9.550452724745e-001     -2.698714095229e-005
	//-1.253465356621e-003     -3.608365188597e-004    9.999991493104e-001
	//IERS res:
	//9.550441611240e-001      2.964607054200e-001     1.304011320000e-003
	//-2.964609236470e-001     9.550449832430e-001     -2.707748000000e-005
	//-1.253416878000e-003     -3.607282110000e-004    9.999991494100e-001
	//Delta IERS
	//2.891654781179e-007      -9.318904692912e-007    7.843450902634e-008
	//9.317494958361e-007      2.892314655556e-007     9.033904770511e-008
	//-4.847862116559e-008     -1.083078596726e-007    -9.955891666635e-011
	//==============================================================================//
	int InfluenceForce::TestRotation()
	{
		printf("Test Rotation\n");
		// 1
		//2012 9 22     17 h 22 min 15 s   UTC 
		//double date1 = 20120922.0;
		//double time1 = 202215.0;//21594247; //(UTC+3:00)
		// 2
		// 2013 2 15     13 h 08 min 53 s   UTC 
		//double date1 = 20130215.0;
		//double time1 = 160853.0;//(UTC+3:00)
		// 3 - new
		// 2013 2 15     13 h 08 min 53 s   UTC 
		//double date1 = 20090102.0;
		//double time1 = 10000.0;//(UTC+3:00)


		double date1 = 20071102.0;
		double time1 = 122212.0;//(UTC+3:00)

		double t1, ajd1, delt1;
		set_time(date1, time1, &ajd1, &delt1, &t1 );

		// матрица поворота из ITRF в ICRF
		// из солнечной системы в земную систему
		double A_rot_to[9];
		iers_update_matrix( t1, A_rot_to, ajd1, delt1 );
		//double A_rot[9];
		//transpose( A_rot_to, A_rot ); 

		printf( "ITRF to ICRF\n" );
		for( int j = 0; j < 9; j += 3 )
		{
			//printf("%.12e\t %.12e\t %.12e\n", A_rot[j+0], A_rot[j+1], A_rot[j+2] );
		}
		printf("\n");

		/*printf("\n DET \n");
		double det = A_rot[0] * A_rot[4] * A_rot[8] + A_rot[1] * A_rot[5] * A_rot[6]
					+ A_rot[2] * A_rot[3] * A_rot[7] - A_rot[2] * A_rot[4] * A_rot[6]
					- A_rot[1] * A_rot[3] * A_rot[8] - A_rot[0] * A_rot[5] * A_rot[7];
		printf("DET = %.12e\n", det);*/

		// Transformation coordinate M from the international terrestrial reference system (ITRF)
		// to the international celestial reference frame (ICRF) : 
		// Celestial coordinates ( X Y Z) = M x Terrestrial coordinates ( x y z )

		// 1
		//2012 9 22     17 h 22 min 15 s   UTC 
        //double A_iers[9] = { -0.132754587260,     0.991148129703,     0.001266707241, 
		//					  -0.991148928306,    -0.132754666371,    -0.000021794431, 
		//					  0.000146559788,    -0.001258388835,     0.999999197489}; 

		 // 2
		 //2013 2 15     13 h 08 min 53 s   UTC 
       	//double A_iers[9] = {   0.955044161124,     0.296460705420 ,    0.001304011320, 
		//					  -0.296460923647,     0.955044983243,    -0.000027077480, 
		//					  -0.001253416878,    -0.000360728211,     0.999999149410   };

		// 3 - new
		//2013 2 15     13 h 08 min 53 s   UTC 
		//double A_iers[9] = { 0.316208771190, -0.948689202655,  0.000899882427,
			//					0.948689578879,  0.316208921793,  0.000026570039,
				//				-0.000309757561,  0.000845307402,  0.999999594753};

		double A_iers[9] = { -0.999547001090,     0.030086434974,     0.000773978848,
							-0.030086414120, -0.999547300630,     0.000038576390,
							0.000774789095,     0.000015272667,     0.999999699734 };

		double A_rot[9] = { -9.99547001036216E-01,      3.00864367479260E-02,      7.73978889005599E-04,
			- 3.00864158935205E-02, - 9.99547300576900E-01,      3.85760783290874E-05,
			7.74789125949687E-04,      1.52723526581853E-05,      9.99999699734238E-01 };

		printf( "IERS res:\n" );
		for( int j = 0; j < 9; j += 3 )
		{
			printf("%.12e\t %.12e\t %.12e\n", A_iers[j+0], A_iers[j+1], A_iers[j+2] );
		}
		printf("\n");

		printf("MyICRF to IERS:\n");
		double R[9];
		//R = A_iers * A-rot ^ (-1)
		for (int j = 0; j < 9; j += 3)
		{
			for (int i = 0; i < 9; i += 3)
			{
				R[j + (i / 3)] = A_iers[j] * A_rot[i] + A_iers[j + 1] * A_rot[i + 1] + A_iers[j + 2] * A_rot[i + 2];
			}
		}

		for (int j = 0; j < 9; j += 3)
		{
			printf("%.12e\t %.12e\t %.12e\n", R[j + 0], R[j + 1], R[j + 2]);
		}
		printf("\n");

		//delta(x,y,z) = (R - E)*A_rot*(x,y,z)
		R[0] = R[0] - 1.0;
		R[4] = R[4] - 1.0;
		R[8] = R[8] - 1.0;

		double DR[9];	//DR = (R-E)*A_rot
		for (int j = 0; j < 9; j += 3)
		{
			for (int i = 0; i < 3; i++)
			{
				DR[j + i] = (R[j] * A_rot[i] + R[j + 1] * A_rot[i + 3] + R[j + 2] * A_rot[i + 6]);
			}
		}

		double xyz[3] = { 0 };	//Test vector in ICRF. If in ITRF, then (x,y,z) = A_rot(x_itrf, y_itrf, z_itrf)

		xyz[0] = 0; xyz[1] = 0; xyz[2] = 6371000;

		double delta_x = xyz[0]*DR[0] + xyz[1]*DR[1] + xyz[2]*DR[2];
		double delta_y = xyz[0]*DR[3] + xyz[1]*DR[4] + xyz[2]*DR[5];
		double delta_z = xyz[0]*DR[6] + xyz[1]*DR[7] + xyz[2]*DR[8];
		printf("delta_x = %.12e\t delta_y = %.12e\t delta_z = %.12e\n", delta_x, delta_y, delta_z);

		printf("IERS: x = %.12e\t y = %.12e\t z = %.12e\n", (xyz[0] * A_iers[0] + xyz[1] * A_iers[1] + xyz[2] * A_iers[2]), ((xyz[0] * A_iers[3] + xyz[1] * A_iers[4] + xyz[2] * A_iers[5])), ((xyz[0] * A_iers[6] + xyz[1] * A_iers[7] + xyz[2] * A_iers[8])));
		printf("MyICRF: x = %.12e\t y = %.12e\t z = %.12e\n", (xyz[0] * A_rot[0] + xyz[1] * A_rot[1] + xyz[2] * A_rot[2]), ((xyz[0] * A_rot[3] + xyz[1] * A_rot[4] + xyz[2] * A_rot[5])), ((xyz[0] * A_rot[6] + xyz[1] * A_rot[7] + xyz[2] * A_rot[8])));

		return 0;
	}
	//==============================================================================//
};
/*
// поправки полюса
//13 1 1 56293.00 P  0.082004 0.005865  0.298188 0.007111  P 0.2857338 0.0055637                 P   -75.505     .600    -8.527     .600
//13 1 2 56294.00 P  0.080977 0.005927  0.298532 0.007205  P 0.2847340 0.0056502                 P   -75.505     .600    -8.480     .600
//13 1 3 56295.00 P  0.079960 0.005988  0.298894 0.007299  P 0.2836315 0.0057362                 P   -75.519     .600    -8.293     .600

//7 11 2 54406.00 I  0.042326 0.000038  0.191057 0.000043  I-0.2138878 0.0000042  0.7657 0.0043  I   -64.490     .328    -5.038     .340   .042290   .191110  -.2139220   -63.600    -4.900

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
*/
int main()
{
	Force::InfluenceForce *IF = new Force::InfluenceForce();
	IF->Init_CPU();

	IF->TestRotation();

	short int error = 0;
	double jd_beg, jd_end;
	short int de_num = 0;
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

	const short int accuracy = 0;
	double vter[3], vcel[3];
	vter[0] = 0; vter[1] = 0; vter[2] = 6371000;

	double x_pole = 0.042326; //(sec.of arc)
	double y_pole = 0.191057;
	double ut1_utc = 0.2857338;

	double leap_secs = 35.0;
	
	double delta_t = 32.184 + leap_secs - ut1_utc;

	double jd_ut1 = 54406.390416667 + 2400000.5;
		
	if ((error = ter2cel(jd_ut1, 0.0, delta_t, 1, accuracy, 0, x_pole, y_pole, vter,
		vcel)) != 0)
	{
		printf("Error %d from ter2cel.", error);
		return (error);
	}
	printf("Ter2cel: x = %.12e\t y = %.12e\t z = %.12e\n", vcel[0], vcel[1], vcel[2]);

	return 0;
}