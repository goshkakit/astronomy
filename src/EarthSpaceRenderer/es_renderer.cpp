#include "es_renderer.h"
#include "es_camera.h"
#include "es_globals.h"

COpenGLRenderer::COpenGLRenderer()
{
	Camera.SetViewMatrixPointer(&ViewMatrix);
}

COpenGLRenderer::~COpenGLRenderer()
{
}

bool COpenGLRenderer::Init()
{
	bool Error = false;

	Error |= !EarthMap.LoadTexture2D("earthmap.jpg");
	Error |= !CloudsMap.LoadTexture2D("cloudsmap.jpg");
	Error |= !LightsMap.LoadTexture2D("lightsmap.jpg");

	Error |= !gfs.Load("groundfromspace.vs", "groundfromspace.fs");
	Error |= !gfa.Load("groundfromatmosphere.vs", "groundfromatmosphere.fs");
	Error |= !sfs.Load("skyfromspace.vs", "skyfromspace.fs");
	Error |= !sfa.Load("skyfromatmosphere.vs", "skyfromatmosphere.fs");

	if (Error)
	{
		return false;
	}

	float Kr = 0.0030f;
	float Km = 0.0015f;
	float ESun = 16.0f;
	float g = -0.75f;
	InnerRadius = 10.0f * 25.0f;
	OuterRadius = 10.25f * 25.0f;
	float Scale = 1.0f / (OuterRadius - InnerRadius);
	float ScaleDepth = 0.25f;
	float ScaleOverScaleDepth = Scale / ScaleDepth;

	GLuint programs[] = { gfs, gfa, sfs, sfa };

	for (int i = 0; i < 4; i++)
	{
		glUseProgram(programs[i]);
		glUniform3f(glGetUniformLocation(programs[i], "v3LightPos"), 0.0f, 0.0f, 1.0f);
		glUniform3f(glGetUniformLocation(programs[i], "v3InvWavelength"), 1.0f / powf(0.650f, 4.0f), 1.0f / powf(0.570f, 4.0f), 1.0f / powf(0.475f, 4.0f));
		glUniform1f(glGetUniformLocation(programs[i], "fInnerRadius"), InnerRadius);
		glUniform1f(glGetUniformLocation(programs[i], "fInnerRadius2"), InnerRadius * InnerRadius);
		glUniform1f(glGetUniformLocation(programs[i], "fOuterRadius"), OuterRadius);
		glUniform1f(glGetUniformLocation(programs[i], "fOuterRadius2"), OuterRadius * OuterRadius);
		glUniform1f(glGetUniformLocation(programs[i], "fKrESun"), Kr * ESun);
		glUniform1f(glGetUniformLocation(programs[i], "fKmESun"), Km * ESun);
		glUniform1f(glGetUniformLocation(programs[i], "fKr4PI"), Kr * 4.0f * (float)M_PI);
		glUniform1f(glGetUniformLocation(programs[i], "fKm4PI"), Km * 4.0f * (float)M_PI);
		glUniform1f(glGetUniformLocation(programs[i], "fScale"), Scale);
		glUniform1f(glGetUniformLocation(programs[i], "fScaleDepth"), ScaleDepth);
		glUniform1f(glGetUniformLocation(programs[i], "fScaleOverScaleDepth"), ScaleOverScaleDepth);
		glUniform1f(glGetUniformLocation(programs[i], "g"), g);
		glUniform1f(glGetUniformLocation(programs[i], "g2"), g * g);
		glUniform1i(glGetUniformLocation(programs[i], "Samples"), 4);
		glUniform1i(glGetUniformLocation(programs[i], "s2Tex1"), 0);
		glUniform1i(glGetUniformLocation(programs[i], "s2Tex2"), 1);
		glUniform1i(glGetUniformLocation(programs[i], "s2Tex3"), 2);
		glUseProgram(0);
	}

	int X = 128, Y = X / 2, vpos = 0, tpos = 0;
	float a, stepa = (float)M_PI * 2.0f / (float)X, stepb = (float)M_PI / (float)Y, b = -(float)M_PI_2 + stepb;

	vec3 *vertices = new vec3[X * (Y - 1)];

	for (int y = 0; y < (Y - 1); y++)
	{
		a = -(float)M_PI;

		for (int x = 0; x < X; x++)
		{
			vertices[y * X + x] = normalize(vec3(sin(a) * cos(b), sin(b), cos(a) * cos(b)));
			a += stepa;
		}

		b += stepb;
	}

	VerticesCount = (X * (Y - 2) * 2 + X * 2) * 3;

	vec3 *Vertices = new vec3[VerticesCount];
	vec2 *TexCoords = new vec2[VerticesCount];

	for (int x = 0; x < X; x++)
	{
		Vertices[vpos++] = vec3(0.0f, -1.0f, 0.0f);
		Vertices[vpos++] = vertices[(0 + 0) * X + ((x + 1) % X)];
		Vertices[vpos++] = vertices[(0 + 0) * X + ((x + 0) % X)];

		TexCoords[tpos++] = vec2((float)(x + 0.5f) / (float)X, 0.0f);
		TexCoords[tpos++] = vec2((float)(x + 1) / (float)X, (float)(0 + 1) / (float)Y);
		TexCoords[tpos++] = vec2((float)(x + 0) / (float)X, (float)(0 + 1) / (float)Y);
	}

	for (int y = 0; y < Y - 2; y++)
	{
		for (int x = 0; x < X; x++)
		{
			Vertices[vpos++] = vertices[(y + 0) * X + ((x + 0) % X)];
			Vertices[vpos++] = vertices[(y + 0) * X + ((x + 1) % X)];
			Vertices[vpos++] = vertices[(y + 1) * X + ((x + 1) % X)];

			TexCoords[tpos++] = vec2((float)(x + 0) / (float)X, (float)(1 + y + 0) / (float)Y);
			TexCoords[tpos++] = vec2((float)(x + 1) / (float)X, (float)(1 + y + 0) / (float)Y);
			TexCoords[tpos++] = vec2((float)(x + 1) / (float)X, (float)(1 + y + 1) / (float)Y);

			Vertices[vpos++] = vertices[(y + 1) * X + ((x + 1) % X)];
			Vertices[vpos++] = vertices[(y + 1) * X + ((x + 0) % X)];
			Vertices[vpos++] = vertices[(y + 0) * X + ((x + 0) % X)];

			TexCoords[tpos++] = vec2((float)(x + 1) / (float)X, float(1 + y + 1) / (float)Y);
			TexCoords[tpos++] = vec2((float)(x + 0) / (float)X, float(1 + y + 1) / (float)Y);
			TexCoords[tpos++] = vec2((float)(x + 0) / (float)X, float(1 + y + 0) / (float)Y);
		}
	}

	for (int x = 0; x < X; x++)
	{
		Vertices[vpos++] = vertices[(Y - 2) * X + ((x + 0) % X)];
		Vertices[vpos++] = vertices[(Y - 2) * X + ((x + 1) % X)];
		Vertices[vpos++] = vec3(0.0f, 1.0f, 0.0f);

		TexCoords[tpos++] = vec2((float)(x + 0) / (float)X, (float)(Y - 1) / (float)Y);
		TexCoords[tpos++] = vec2((float)(x + 1) / (float)X, (float)(Y - 1) / (float)Y);
		TexCoords[tpos++] = vec2((float)(x + 0.5f) / (float)X, 1.0f);
	}

	glGenBuffers(1, &VerticesVBO);
	glBindBuffer(GL_ARRAY_BUFFER, VerticesVBO);
	glBufferData(GL_ARRAY_BUFFER, VerticesCount * 3 * 4, Vertices, GL_STATIC_DRAW);

	glGenBuffers(1, &TexCoordsVBO);
	glBindBuffer(GL_ARRAY_BUFFER, TexCoordsVBO);
	glBufferData(GL_ARRAY_BUFFER, VerticesCount * 2 * 4, TexCoords, GL_STATIC_DRAW);

	glBindBuffer(GL_ARRAY_BUFFER, 0);

	delete[] vertices;
	delete[] Vertices;
	delete[] TexCoords;

	Camera.Look(vec3(0.0f, 0.0f, 768.0f), vec3(0.0f, 0.0f, 0.0f), InnerRadius);

	return true;
}

void COpenGLRenderer::Render(float FrameTime)
{
	static float a = 0.0f;

	//float ainc = 0.25f * FrameTime;

	//a += ainc;

	Space::m3x3d earthModelRotation = EarthModel.Update(TimeModel.now());

	//if (length(Camera.Position) < OuterRadius * 1.1f)
	if (bCameraRotationEnabled)
	{
		//vec3 Y = vec3(0.0f, 1.0f, 0.0f);

		//Camera.Position = rotate(Camera.Position, ainc, Y);

		//Camera.X = rotate(Camera.X, ainc, Y);
		//Camera.Y = rotate(Camera.Y, ainc, Y);
		//Camera.Z = rotate(Camera.Z, ainc, Y);

		Camera.Position = EarthModel.rotateVec3(Camera.Position, earthModelRotation);

		Camera.X = EarthModel.rotateVec3(Camera.X, earthModelRotation);
		Camera.Y = EarthModel.rotateVec3(Camera.Y, earthModelRotation);
		Camera.Z = EarthModel.rotateVec3(Camera.Z, earthModelRotation);

		Camera.CalculateViewMatrix();
	}

	mat4x4 sceneRotation = EarthModel.getMat4x4(EarthModel.EarthAtNow().orientation());

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);

	GLuint program;

	glMatrixMode(GL_MODELVIEW);
	glLoadMatrixf(&ViewMatrix);
	glScalef(InnerRadius, InnerRadius, InnerRadius);
	//glRotatef(a, 0.0f, 1.0f, 0.0f);
	glMultMatrixf(&sceneRotation);

	glMatrixMode(GL_TEXTURE);
	glLoadIdentity();
	glScalef(InnerRadius, InnerRadius, InnerRadius);
	//glRotatef(a, 0.0f, 1.0f, 0.0f);
	glMultMatrixf(&sceneRotation);

	glActiveTexture(GL_TEXTURE0); glBindTexture(GL_TEXTURE_2D, EarthMap);
	glActiveTexture(GL_TEXTURE1); glBindTexture(GL_TEXTURE_2D, CloudsMap);
	glActiveTexture(GL_TEXTURE2); glBindTexture(GL_TEXTURE_2D, LightsMap);

	if (length(Camera.Position) > OuterRadius) program = gfs;
	else program = gfa;

	glUseProgram(program);
	glUniform3fv(glGetUniformLocation(program, "v3CameraPos"), 1, &Camera.Position);
	glUniform1f(glGetUniformLocation(program, "fCameraHeight"), length(Camera.Position));
	glUniform1f(glGetUniformLocation(program, "fCameraHeight2"), length(Camera.Position) * length(Camera.Position));

	glEnableClientState(GL_VERTEX_ARRAY);
	glBindBuffer(GL_ARRAY_BUFFER, VerticesVBO);
	glVertexPointer(3, GL_FLOAT, 0, NULL);

	glEnableClientState(GL_TEXTURE_COORD_ARRAY);
	glBindBuffer(GL_ARRAY_BUFFER, TexCoordsVBO);
	glTexCoordPointer(2, GL_FLOAT, 0, NULL);

	glDrawArrays(GL_TRIANGLES, 0, VerticesCount);

	glDisableClientState(GL_TEXTURE_COORD_ARRAY);
	glDisableClientState(GL_VERTEX_ARRAY);

	glUseProgram(0);

	glActiveTexture(GL_TEXTURE2); glBindTexture(GL_TEXTURE_2D, 0);
	glActiveTexture(GL_TEXTURE1); glBindTexture(GL_TEXTURE_2D, 0);
	glActiveTexture(GL_TEXTURE0); glBindTexture(GL_TEXTURE_2D, 0);

	glMatrixMode(GL_MODELVIEW);
	glLoadMatrixf(&ViewMatrix);
	glScalef(OuterRadius, OuterRadius, OuterRadius);

	glMatrixMode(GL_TEXTURE);
	glLoadIdentity();
	glScalef(OuterRadius, OuterRadius, OuterRadius);

	if (length(Camera.Position) > OuterRadius) program = sfs;
	else program = sfa;

	glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_BLEND);

	glCullFace(GL_FRONT);

	glUseProgram(program);
	glUniform3fv(glGetUniformLocation(program, "v3CameraPos"), 1, &Camera.Position);
	glUniform1f(glGetUniformLocation(program, "fCameraHeight"), length(Camera.Position));
	glUniform1f(glGetUniformLocation(program, "fCameraHeight2"), length(Camera.Position) * length(Camera.Position));

	glEnableClientState(GL_VERTEX_ARRAY);
	glBindBuffer(GL_ARRAY_BUFFER, VerticesVBO);
	glVertexPointer(3, GL_FLOAT, 0, NULL);

	glEnableClientState(GL_TEXTURE_COORD_ARRAY);
	glBindBuffer(GL_ARRAY_BUFFER, TexCoordsVBO);
	glTexCoordPointer(2, GL_FLOAT, 0, NULL);

	glDrawArrays(GL_TRIANGLES, 0, VerticesCount);

	glDisableClientState(GL_TEXTURE_COORD_ARRAY);
	glDisableClientState(GL_VERTEX_ARRAY);

	glUseProgram(0);

	glCullFace(GL_BACK);

}

void COpenGLRenderer::RenderFinish()
{
	glDisable(GL_CULL_FACE);
	glDisable(GL_DEPTH_TEST);
}

void COpenGLRenderer::Resize(int Width, int Height)
{
	this->Width = Width;
	this->Height = Height;

	glViewport(0, 0, Width, Height);

	ProjectionMatrix = perspective(45.0f, (float)Width / (float)Height, 0.125f, 4096.0f);

	glMatrixMode(GL_PROJECTION);
	glLoadMatrixf(&ProjectionMatrix);
}

void COpenGLRenderer::Destroy()
{
	EarthMap.Destroy();
	CloudsMap.Destroy();
	LightsMap.Destroy();

	glDeleteBuffers(1, &VerticesVBO);
	glDeleteBuffers(1, &TexCoordsVBO);

	gfs.Destroy();
	gfa.Destroy();
	sfs.Destroy();
	sfa.Destroy();
}
