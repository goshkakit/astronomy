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
	};

	class CScene
	{
	public:
		std::vector<CLayer *> layers;

		CScene();
		virtual ~CScene();
		virtual void Update(float time);
		virtual void Render();
		virtual void Clear();
		
		void FreeInstances();

	private:
		CScene(const CScene &) = delete;
		CScene & operator=(const CScene &) = delete;
	};

	class CPointArrayLayer : public CLayer
	{
	public:
		std::vector<vec3> points;
		vec3 pointColor;
		GLfloat pointSize;

		CPointArrayLayer();
		virtual ~CPointArrayLayer();
		//virtual void Update(float time);
		virtual void Render();
		virtual void Clear();
	};
}
