#include "es_scene.h"

Scene::CLayer::CLayer()
{

}

Scene::CLayer::~CLayer()
{

}

void Scene::CLayer::Update(float time)
{

}

void Scene::CLayer::Render()
{

}

void Scene::CLayer::Clear()
{

}

Scene::CLayerGroup::CLayerGroup()
	: CLayer()
	, layers()
{

}

Scene::CLayerGroup::~CLayerGroup()
{
	for (CLayer *layer : layers)
	{
		delete layer;
	}
}

void Scene::CLayerGroup::Update(float time)
{
	for (CLayer *layer : layers)
	{
		layer->Update(time);
	}
}

void Scene::CLayerGroup::Render()
{
	for (CLayer *layer : layers)
	{
		layer->Render();
	}
}

void Scene::CLayerGroup::Clear()
{
	for (CLayer *layer : layers)
	{
		layer->Clear();
	}
}

Scene::CScene::CScene()
	: CLayerGroup()
{

}

Scene::CScene::~CScene()
{

}

void Scene::CScene::FreeInstances()
{
	for (CLayer *layer : layers)
	{
		delete layer;
	}

	layers.clear();
}

Scene::CPointArrayLayer::CPointArrayLayer()
	: CLayer()
	, points()
	, pointColor(1.f, 0.f, 0.f)
	, pointSize(3.f)
{

}

Scene::CPointArrayLayer::~CPointArrayLayer()
{

}

void Scene::CPointArrayLayer::Render()
{
	if (!points.size()) return;

	glColor3fv(&pointColor);
	glPointSize(pointSize);
	glBegin(GL_POINTS);

	for (vec3 &point : points)
	{
		glVertex3fv(&point);
	}

	glEnd();
}

void Scene::CPointArrayLayer::Clear()
{
	points.clear();
}

Scene::CLineLayer::CLineLayer()
	: CLayer()
	, vertices{ vec3(), vec3() }
	, lineColor(1.f, 0.f, 0.f)
	, lineSize(1.f)
{

}

Scene::CLineLayer::~CLineLayer()
{

}

void Scene::CLineLayer::Render()
{
	glColor3fv(&lineColor);
	glLineWidth(lineSize);
	glBegin(GL_LINES);
	glVertex3fv(&vertices[0]);
	glVertex3fv(&vertices[1]);
	glEnd();
}

Scene::CPolylineLayer::CPolylineLayer()
	: CLayer()
	, vertices()
	, lineColor(1.f, 0.f, 0.f)
	, lineSize(1.f)
{

}

Scene::CPolylineLayer::~CPolylineLayer()
{

}

void Scene::CPolylineLayer::Render()
{
	if (vertices.size() > 1)
	{
		glColor3fv(&lineColor);
		glLineWidth(lineSize);
		glBegin(GL_LINES);

		glVertex3fv(&*vertices.begin());

		if (vertices.size() > 2)
		{
			for (vec3 &vertex : vertices)
			{
				glVertex3fv(&vertex);
				glVertex3fv(&vertex);
			}
		}

		glVertex3fv(&*vertices.rbegin());

		glEnd();
	}
}

void Scene::CPolylineLayer::Clear()
{
	vertices.clear();
}
