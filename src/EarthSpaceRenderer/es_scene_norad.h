#pragma once

#include "es_scene.h"

#include "common/tleloader.h"

namespace Scene
{
	class CNoradLayer : public CPointArrayLayer
	{
	public:
		TLELoader tleloader;

		void Load(const char *tlepath, int numbLine);
		virtual void Update(float time);

		CNoradLayer();
		CNoradLayer(const char *tlepath, int numbLine);
		virtual ~CNoradLayer();
	};
}
