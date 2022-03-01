#pragma once

#include "es_defines.h"

#include <vector>

namespace Scene
{
	class CLayer
	{
	public:
		CLayer();
		virtual ~CLayer();
		virtual void Update(float time);
		virtual void Render();
		virtual void Clear();

	private:
		CLayer(const CLayer &) = delete;
		CLayer & operator=(const CLayer &) = delete;
	};

	class CLayerGroup : public CLayer
	{
	public:
		std::vector<CLayer *> layers;

		CLayerGroup();
		virtual ~CLayerGroup();
		virtual void Update(float time);
		virtual void Render();
		virtual void Clear();
	};

	class CScene : public CLayerGroup
	{
	public:
		CScene();
		virtual ~CScene();
		void FreeInstances();
	};

	class CPointArrayLayer : public CLayer
	{
	public:
		std::vector<vec3> points;
		vec3 pointColor;
		GLfloat pointSize;

		CPointArrayLayer();
		virtual ~CPointArrayLayer();
		virtual void Render();
		virtual void Clear();
	};

	class CLineLayer : public CLayer
	{
	public:
		vec3 vertices[2];
		vec3 lineColor;
		GLfloat lineSize;

		CLineLayer();
		virtual ~CLineLayer();
		virtual void Render();
	};

	class CPolylineLayer : public CLayer
	{
	public:
		std::vector<vec3> vertices;
		vec3 lineColor;
		GLfloat lineSize;

		CPolylineLayer();
		virtual ~CPolylineLayer();
		virtual void Render();
		virtual void Clear();
	};
}
