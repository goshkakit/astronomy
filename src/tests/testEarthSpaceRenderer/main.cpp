
#include "EarthSpaceRenderer/es_globals.h"
#include "EarthSpaceRenderer/es_scene_norad.h"
#include "EarthSpaceRenderer/es_scene_orbit.h"
#include "EarthSpaceRenderer/corbit_norad.h"
#include "EarthSpaceRenderer/corbit_predict.h"
#include "EarthSpaceRenderer/space_defines.h"
#include "EarthSpaceRenderer/cearth.h"

#include "common/tleloader.h"
#include "common/DataConverter.h"
#include "OrbitIntegration/IPredictOrbitMod.h"

#include <stdio.h>
#include <math.h>
#include <Windows.h>

volatile BOOL bConsole = false;
volatile HWND hConsole = (HWND)NULL;

void FreeConsole_atexit()
{
	if (bConsole > 0) FreeConsole();
}

class CStaticPointOnEarth : public Scene::IOrbit
{
public:
	Space::CEarth *const Earth;
	const double jd_start, lat, lng, alt;

	Space::vec3d v0, v1;
	Space::m3x3d m0t, m1, m;

	CStaticPointOnEarth(Space::CEarth * pEarth, double jd, double _lat, double _lng, double _alt) : Earth(pEarth), jd_start(jd), lat(_lat), lng(_lng), alt(_alt)
	{
		Earth->setTime(jd);
		m0t = Earth->orientation();

		v0 = (Space::mrotZ(lng * M_PI / 180.) * (Space::mrotY(-lat * M_PI / 180.) * Space::vec3d(1., 0., 0.))) * ((EARTH_RADIUS_M + alt) * 1e-6);
	}

	virtual struct Scene::vec3d_t GetNewPosition(double time)
	{
		const double MIN = 1.0 / 1440.0;
		double jd = time * MIN + jd_start;

		Earth->setTime(jd);
		m1 = Earth->rotator();
		m = m1*m0t;

		v1 = m*v0;

		return Scene::vec3d_t(v1.v);
	}
};

int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR sCmdLine, int iShow)
{
	char *AppName = "Planet Earth from space";

	atexit(FreeConsole_atexit);

	if (bConsole = AllocConsole())
	{
		SetConsoleTitle(AppName);

		if (hConsole = GetConsoleWindow())
		{
			freopen("CONIN$", "r", stdin);
			freopen("CONOUT$", "w", stdout);
			freopen("CONOUT$", "w", stderr);
		}
	}

	{
		Space::CEarth Earth(2456339.04783564825);
		Space::m3x3d m = Earth.rotator();

		printf("\nCalculated:\n");
		printf("% .12f % .12f % .12f\n", m.m[0][0], m.m[0][1], m.m[0][2]);
		printf("% .12f % .12f % .12f\n", m.m[1][0], m.m[1][1], m.m[1][2]);
		printf("% .12f % .12f % .12f\n", m.m[2][0], m.m[2][1], m.m[2][2]);
		printf("norma = %g\n", m.norme() / sqrt(3.));

		Space::m3x3d r(
			0.955044161124, 0.296460705420, 0.001304011320,
			-0.296460923647, 0.955044983243, -0.000027077480,
			-0.001253416878, -0.000360728211, 0.999999149410
		);

		printf("\nRight:\n");
		printf("% .12f % .12f % .12f\n", r.m[0][0], r.m[0][1], r.m[0][2]);
		printf("% .12f % .12f % .12f\n", r.m[1][0], r.m[1][1], r.m[1][2]);
		printf("% .12f % .12f % .12f\n", r.m[2][0], r.m[2][1], r.m[2][2]);
		printf("norma = %g\n", r.norme() / sqrt(3.));

		Space::m3x3d d = m*r.transposed() - Space::munit();
		printf("\nDelta:\n");
		printf("% .12f % .12f % .12f\n", d.m[0][0], d.m[0][1], d.m[0][2]);
		printf("% .12f % .12f % .12f\n", d.m[1][0], d.m[1][1], d.m[1][2]);
		printf("% .12f % .12f % .12f\n", d.m[2][0], d.m[2][1], d.m[2][2]);
		printf("norma = %g\n", d.norme() / sqrt(3.));

		Space::vec3d v0 = Space::vec3d(1., 1., 0.).normalized() * EARTH_RADIUS_M;
		Space::vec3d vc = m*v0;
		Space::vec3d vr = r*v0;
		Space::vec3d dr = vc - vr;

		printf("\nError on equator = %g m\n", dr.abs());

		v0 = Space::vec3d(0., 0., 1.) * EARTH_RADIUS_M;
		vc = m*v0;
		vr = r*v0;
		dr = vc - vr;

		printf("\nError on north pole = %g m\n", dr.abs());

		#ifdef _DEBUG
		printf("\nPress ENTER to continue...\n");
		getc(stdin);
		#endif
	}

	IPredictOrbitMod *Predict = CreatePredictOrbitMod();
	Predict->Init();

	double jd_start = DataConverter().YYYYMMDDtoJD(20180724.) - 0.5;
	double atm = 0.3E-2;
	double sun = 0.5E-05;
	double statevec[6];
	{
		const double Rz = EARTH_RADIUS_KM; // km
		const double mu = EARTH_GRAVITATIONAL_PARAMETER; // km3/s2
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

	Space::CEarth Earth;

	CStaticPointOnEarth kislovodsk25(&Earth, jd_start, 43.740331, 42.653458, 2085.0013);

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
		scene.layers.push_back(dynamic_cast<Scene::CLayer *>(
			new Scene::COrbitLayer(dynamic_cast<Scene::IOrbit *>(&kislovodsk25))
			)
		);
		{
			Scene::COrbitLayer *point = dynamic_cast<Scene::COrbitLayer*>(*scene.layers.rbegin());
			point->lineSize = 0.125f;
			point->pointSize = 15.f;
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
