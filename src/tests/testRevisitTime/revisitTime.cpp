#include <stdio.h>

#include "Norad\coreLib.h"
#include "Norad\cOrbit.h"

#include "common\DataConverter.h"
#include "common\mytypes.h"


#include "InfluenceForce\InfluenceForce.h"

void main()
{
	// 1 90732U 00000A   18206.00000000  .00000000  00000-0  58106+0 0    11
	// 2 90732  19.0914   8.4920 2172097 149.9634 343.3148  1.01981561000013
	std::string str1 = "90732";
	std::string str2 = "1 90732U 00000A   18206.00000000  .00000000  00000-0  58106+0 0    11";
	std::string str3 = "2 90732  19.0914   8.4920 2172097 149.9634 343.3148  1.01981561000013";

	// init rotation
	Force::InfluenceForce *FI = new Force::InfluenceForce();
	FI->Init_CPU();

	// init norad
	cTle tleSGP4(str1, str2, str3);

	// vector
	cOrbit *orbit = new cOrbit(tleSGP4);

	// Now create a site object. Site objects represent a location on the 
	// surface of the earth. Here we arbitrarily select a point on the
	// equator.
	cSite siteEquator(0.0, -100.0, 0); // 0.00 N, 100.00 W, 0 km altitude

	// Calculate position, velocity
	// mpe = "minutes past epoch"
	for (int mpe = 0; mpe <= (360 * 4); mpe += 360)
	{
		// Get the position of the satellite at time "mpe"
		cEciTime eci = orbit->GetPosition(mpe);

		// Now get the "look angle" from the site to the satellite. 
		// Note that the ECI object "eciSDP4" contains a time associated
		// with the coordinates it contains; this is the time at which
		// the look angle is valid.
		cTopo topoLook = siteEquator.GetLookAngle(eci);

		// Print out the results.
		printf("AZ: %.3f  EL: %.3f\n",	topoLook.AzimuthDeg(),	topoLook.ElevationDeg());

	}
	
}