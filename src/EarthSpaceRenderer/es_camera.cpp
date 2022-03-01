#include "es_camera.h"

CCamera::CCamera()
{
	ViewMatrix = NULL;
	ViewMatrixInverse = NULL;

	X = vec3(1.0f, 0.0f, 0.0f);
	Y = vec3(0.0f, 1.0f, 0.0f);
	Z = vec3(0.0f, 0.0f, 1.0f);

	Position = vec3(0.0f, 0.0f, 0.0f);
}

CCamera::~CCamera()
{
}

void CCamera::CalculateViewMatrix()
{
	if (ViewMatrix != NULL)
	{
		*ViewMatrix = mat4x4(X.x, Y.x, Z.x, 0.0f, X.y, Y.y, Z.y, 0.0f, X.z, Y.z, Z.z, 0.0f, -dot(X, Position), -dot(Y, Position), -dot(Z, Position), 1.0f);

		if (ViewMatrixInverse != NULL)
		{
			*ViewMatrixInverse = inverse(*ViewMatrix);
		}
	}
}

void CCamera::Look(vec3 Position, vec3 Reference, float PlanetRadius)
{
	this->Position = Position;
	this->PlanetRadius = PlanetRadius;

	Z = normalize(Position - Reference);

	vec3 T = cross(vec3(0.0f, 1.0f, 0.0f), Z);
	if (length2(T) < 1.e-6f)
		 T = cross(vec3(0.0f, 0.0f, ((Z.y < 0.f) ? (1.f) : (-1.f))), Z);

	X = normalize(T);
	Y = cross(Z, X);

	CalculateViewMatrix();
}

void CCamera::OnKeys(SHORT Keys, float FrameTime)
{
	float Speed = max(1.0f, min(PlanetRadius / 2.0f, length(Position) - PlanetRadius));

	if (Keys & 0x100) Speed *= 2.0f;
	if (Keys & 0x200) Speed *= 0.5f;

	float Distance = Speed * FrameTime;

	vec3 Up = Y * Distance;
	vec3 Right = X * Distance;
	vec3 Forward = -Z * Distance;

	vec3 Movement;

	if (Keys & 0x001) Movement += Forward;
	if (Keys & 0x002) Movement -= Forward;
	if (Keys & 0x004) Movement -= Right;
	if (Keys & 0x008) Movement += Right;
	if (Keys & 0x010) Movement += Up;
	if (Keys & 0x020) Movement -= Up;

	if (Keys & 0x040)
	{
		X = rotate(X, 45.0f * FrameTime, Z);
		Y = rotate(Y, 45.0f * FrameTime, Z);
	}

	if (Keys & 0x080)
	{
		X = rotate(X, -45.0f * FrameTime, Z);
		Y = rotate(Y, -45.0f * FrameTime, Z);
	}

	Position += Movement;

	if (length(Position) < PlanetRadius + 0.125f)
	{
		Position = normalize(Position) * (PlanetRadius + 0.125f);
	}

	CalculateViewMatrix();
}

void CCamera::OnMouseMove(int dx, int dy)
{
	float Sensitivity = 0.25f;

	float DeltaX = (float)dx * Sensitivity;

	X = rotate(X, DeltaX, Y);
	Z = rotate(Z, DeltaX, Y);

	float DeltaY = (float)dy * Sensitivity;

	Y = rotate(Y, DeltaY, X);
	Z = rotate(Z, DeltaY, X);

	if ((GetKeyState(VK_LBUTTON) & 0x80))
	{
		Position = rotate(Position, DeltaX, Y);
		Position = rotate(Position, DeltaY, X);
	}

	CalculateViewMatrix();
}

void CCamera::SetViewMatrixPointer(float *ViewMatrix, float *ViewMatrixInverse)
{
	this->ViewMatrix = (mat4x4*)ViewMatrix;
	this->ViewMatrixInverse = (mat4x4*)ViewMatrixInverse;

	CalculateViewMatrix();
}
