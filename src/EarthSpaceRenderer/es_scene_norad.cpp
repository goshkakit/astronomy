#include "es_scene_norad.h"
#include "es_globals.h"

Scene::CNoradLayer::CNoradLayer()
	: CPointArrayLayer()
	, tleloader()
{

}

Scene::CNoradLayer::CNoradLayer(const char *tlepath)
	: CPointArrayLayer()
	, tleloader()
{
	tleloader.LoadData(tlepath, 2);
}

void Scene::CNoradLayer::Load(const char *tlepath)
{
	tleloader.clear();
	tleloader.LoadData(tlepath, 2);
}

Scene::CNoradLayer::~CNoradLayer()
{
	tleloader.clear();
}

void Scene::CNoradLayer::Update(float time)
{
	const double C = 4.0 * ((double)Camera.PlanetRadius) / 6356000.0;

	points.clear();
	points.reserve(tleloader.NORADList.size());

	for (Zeptomoby::OrbitTools::cOrbit * elem : tleloader.NORADList)
	{
		Zeptomoby::OrbitTools::cEciTime eci1 = elem->GetPosition(time*3.0f);

		const Zeptomoby::OrbitTools::cVector & pos = eci1.Position();

		points.push_back(
			vec3( (float)(pos.m_x * C)
				, (float)(pos.m_y * C)
				, (float)(pos.m_z * C)
			)
		);
	}
}
