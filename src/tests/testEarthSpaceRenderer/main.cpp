
#include "EarthSpaceRenderer/es_globals.h"
#include "EarthSpaceRenderer/es_scene_norad.h"
#include "EarthSpaceRenderer/es_scene_orbit.h"
#include "EarthSpaceRenderer/corbit_norad.h"
#include "EarthSpaceRenderer/corbit_predict.h"

#include "common/tleloader.h"
#include "common/DataConverter.h"
#include "OrbitIntegration/IPredictOrbitMod.h"

#include <math.h>
#include <Windows.h>

int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR sCmdLine, int iShow)
{
	char *AppName = "Planet Earth from space";

	IPredictOrbitMod *Predict = CreatePredictOrbitMod();
	Predict->Init();

	double jd_start = DataConverter().YYYYMMDDtoJD(20180724.) - 0.5;
	double atm = 0.3E-2;
	double sun = 0.5E-05;
	double statevec[6];
	{
		const double Rz = EARTH_RADIUS*1e3; // km
		const double mu = 398600.44188; // km3/s2
		const double T = 3600.0 * 2; // s
		const double Ro = cbrt((mu*T*T) / (4.0*M_PI*M_PI)); // km
		const double Vo = sqrt(mu / Ro); // km/s

		statevec[0] = 0.0;
		statevec[1] = 0.0;
		statevec[2] = Ro * 1e-3; // 10^3 km
		statevec[3] = Vo * 0.5 * sqrt(3.); // km/s
		statevec[4] = Vo * 0.5; // km/s
		statevec[5] = 0.0;
	}

	COrbit_Predict Orbit_Predict(Predict, jd_start, statevec, sun, atm);

	if (OpenGLView.Init(hInstance, AppName, 800, 600, 4))
	{
		Scene::CScene scene;
		scene.layers.push_back(dynamic_cast<Scene::CLayer *>(
			new Scene::CNoradLayer("data/test/TLE20180724.txt", 2)
			)
		);
		scene.layers.push_back(dynamic_cast<Scene::CLayer *>(
			new Scene::COrbitLayer(dynamic_cast<Scene::IOrbit *>(&Orbit_Predict))
			)
		);
		{
			dynamic_cast<Scene::COrbitLayer*>(*scene.layers.rbegin())->lineColor = vec3(1.0f, 0.75f, 0.25f);
		}
		OpenGLView.setScene(&scene);

		OpenGLView.Show();
		OpenGLView.MessageLoop();

		scene.FreeInstances();
	}
	else
	{
		MessageBox(NULL, ErrorLog, AppName, MB_OK | MB_ICONERROR);
	}

	OpenGLView.Destroy();
	Predict->DeInit();
	FreePredictOrbitMod(Predict);

	return 0;
}
