#include <stdio.h>
#include "TrackCalculation.h"
#include "AstroDataManager.h"
#include "TelescopeAling.h"
#include "NucleoConnection.h"

int main_TelescopeAling()
{
	TelescopeAling TA;
	TA.Init();
	//TA.calculatePoints();
	TA.aline3DPoints();

	return 0;
}

int main()
{
	NucleoConnection nucleoConnection;
	nucleoConnection.loadTask();
	nucleoConnection.runTest();
}

int main_task()
{
	printf("Start\n");

	string DBRoot = "data/telescope";
	std::string pinpf_tels = DBRoot + "/tels_simulation.json";
	std::string tle_input = DBRoot + "/TLE20180724.txt";

	// init track calculator
	TrackCalculation *trackCalculation = new TrackCalculation();
	trackCalculation->Init();

	// load data
	AstroDataManager *astroDataManager = new AstroDataManager();
	astroDataManager->loadTelescope(pinpf_tels);
	astroDataManager->loadTLE(tle_input);

	if (astroDataManager->tels.size() > 0)
	{
		TelescopObject tel = astroDataManager->tels[0];

		int idSat = 28646;
		cOrbit* orb = astroDataManager->tleLoader.GetFromID(idSat);
		if (orb != NULL)
		{
			trackCalculation->PrepareTracks(tel, orb);
		}
	}
	return 0;
}

