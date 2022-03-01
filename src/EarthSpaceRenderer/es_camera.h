#pragma once

#include "es_defines.h"

class CCamera
{
protected:
	mat4x4 *ViewMatrix, *ViewMatrixInverse;

public:
	vec3 X, Y, Z, Position;
	float PlanetRadius;

	CCamera();
	~CCamera();

	void CalculateViewMatrix();
	void Look(vec3 Position, vec3 Reference, float PlanetRadius);
	void OnKeys(SHORT Keys, float FrameTime);
	void OnMouseMove(int dx, int dy);
	void SetViewMatrixPointer(float *ViewMatrix, float *ViewMatrixInverse = NULL);
};
