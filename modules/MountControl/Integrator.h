#pragma once
#include <utility>
#include <memory>

#include "InfluenceForce\InfluenceForce.h"
#include "OrbitIntegration\IPredictOrbitMod.h"
#include "OrbitIntegration\PredictOrbitMod.h"
#include "common\DataConverter.h"
#include "common\TLELoader.h"

class Integrator {
public:
	Integrator();
	~Integrator();

	std::unique_ptr<Force::InfluenceForce> IForce1;
private:
	Orbit::PredictOrbitSat POSat1;
	DataConverter DC1;
};