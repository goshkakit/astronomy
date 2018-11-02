#include "es_scene.h"

Scene::CScene::CScene()
	: layers()
{

}

Scene::CScene::~CScene()
{

}

void Scene::CScene::Update(float time)
{
	for (CLayer *layer : layers)
	{
		layer->Update(time);
	}
}

void Scene::CScene::Render()
{
	for (CLayer *layer : layers)
	{
		layer->Render();
	}
}

void Scene::CScene::Clear()
{
	for (CLayer *layer : layers)
	{
		layer->Clear();
	}
}

void Scene::CScene::FreeInstances()
{
	for (CLayer *layer : layers)
	{
		delete layer;
	}

	layers.clear();
}

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
