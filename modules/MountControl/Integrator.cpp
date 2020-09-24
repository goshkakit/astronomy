#include "integrator.h"

extern "C" {
#include "novac/novas.h"
#include "novac/eph_manager.h"
}

Integrator::Integrator() {
	IForce1 = std::make_unique<Force::InfluenceForce>();
	IForce1->Init_CPU("Data/TAI-UTC.DAT",
					  "Data/DE403.bin",
					  "Data/finals.all",
					  "Data/atm.config");
	POSat1.Init_CPU();

	//—читывание эфемерид движени€ планет
	short int error = 0;
	short int de_num = 0;
	double jd_beg, jd_end;
	if ((error = ephem_open("Data/eph/lnx1900.405", &jd_beg, &jd_end, &de_num)) != 0)
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
}

Integrator::~Integrator() {
}