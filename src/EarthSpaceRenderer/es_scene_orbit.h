#pragma once

#include "es_scene.h"
#include "iorbit.h"

namespace Scene
{
	class COrbitLayer : public CPolylineLayer
	{
	public:
		GLfloat pointSize;
		IOrbit * const pOrbit;

		COrbitLayer(IOrbit * _pOrbit);
		virtual ~COrbitLayer();
		virtual void Update(float time);
		virtual void Render();
	};
}
