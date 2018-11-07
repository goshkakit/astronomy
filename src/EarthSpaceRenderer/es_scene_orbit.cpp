#include "es_scene_orbit.h"

Scene::COrbitLayer::COrbitLayer(IOrbit * _pOrbit)
	: CPolylineLayer()
	, pointSize(3.f)
	, pOrbit(_pOrbit)
{

}

Scene::COrbitLayer::~COrbitLayer()
{

}

void Scene::COrbitLayer::Update(float time)
{
	const double C = 1. / EARTH_RADIUS;

	vec3d_t pos = pOrbit->GetNewPosition(time);

	vertices.push_back(
		vec3( (float)(pos.y * C)
			, (float)(pos.z * C)
			, (float)(pos.x * C)
		)
	);
}

void Scene::COrbitLayer::Render()
{
	CPolylineLayer::Render();

	if (vertices.size())
	{
		glColor3fv(&lineColor);
		glPointSize(pointSize);
		glBegin(GL_POINTS);
		glVertex3fv(&*vertices.rbegin());
		glEnd();
	}
}
