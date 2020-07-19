#pragma once

#include "iorbit.h"

#include "Norad/std_add.h"
#include "Norad/cOrbit.h"

class COrbit_Norad : public Scene::IOrbit
{
public:
	COrbit_Norad(Zeptomoby::OrbitTools::cOrbit * _pOrbit);

	virtual struct Scene::vec3d_t GetNewPosition(double time);

protected:
	Zeptomoby::OrbitTools::cOrbit * const pOrbit;
};
