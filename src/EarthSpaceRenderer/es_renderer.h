#pragma once

#include "es_defines.h"
#include "es_texture.h"

class COpenGLRenderer
{
protected:
	int Width, Height;
	mat3x3 NormalMatrix;
	mat4x4 ModelMatrix, ViewMatrix, ProjectionMatrix;

protected:
	CTexture EarthMap, CloudsMap, LightsMap;
	CShaderProgram gfs, gfa, sfs, sfa;
	GLuint VerticesVBO, TexCoordsVBO, VerticesCount;
	float InnerRadius, OuterRadius;

public:
	CString Text;

public:
	COpenGLRenderer();
	~COpenGLRenderer();

	bool Init();
	void Render(float FrameTime);
	void RenderFinish();
	void Resize(int Width, int Height);
	void Destroy();
};
