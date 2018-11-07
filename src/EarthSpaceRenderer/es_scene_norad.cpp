#include "es_scene_norad.h"
#include "es_globals.h"

Scene::CNoradLayer::CNoradLayer()
	: CPointArrayLayer()
	, tleloader()
{

}

Scene::CNoradLayer::CNoradLayer(const char *tlepath, int numbLine)
	: CPointArrayLayer()
	, tleloader()
{
	tleloader.LoadData(tlepath, numbLine);
}

void Scene::CNoradLayer::Load(const char *tlepath, int numbLine)
{
	tleloader.clear();
	tleloader.LoadData(tlepath, numbLine);
}

Scene::CNoradLayer::~CNoradLayer()
{
	tleloader.clear();
}

void Scene::CNoradLayer::Update(float time)
{
	const double C = 4.0e-6 * ((double)Camera.PlanetRadius) / EARTH_RADIUS;

	points.clear();
	points.reserve(tleloader.NORADList.size());

	for (Zeptomoby::OrbitTools::cOrbit * elem : tleloader.NORADList)
	{
		Zeptomoby::OrbitTools::cEciTime eci1 = elem->GetPosition(time);

		const Zeptomoby::OrbitTools::cVector & pos = eci1.Position();

		points.push_back(
			vec3( (float)(pos.m_y * C)
				, (float)(pos.m_z * C)
				, (float)(pos.m_x * C)
			)
		);
	}
}
