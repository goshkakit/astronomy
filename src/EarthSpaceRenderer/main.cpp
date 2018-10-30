#include "main.h"
#include <stddef.h> 

#include "astronomy/tleloader.h"

TwType TW_TYPE_OGLDEV_VECTOR3F;
TwType TW_TYPE_OGLDEV_VECTOR2I;
TwType TW_TYPE_OGLDEV_ATTENUATION;

int gl_max_texture_size = 0, gl_max_texture_max_anisotropy_ext = 0;

CCamera Camera;

COpenGLRenderer OpenGLRenderer;

CString ModuleDirectory, ErrorLog;

COpenGLView OpenGLView;

void GetModuleDirectory()
{
	char *moduledirectory = new char[256];
	GetModuleFileName(GetModuleHandle(NULL), moduledirectory, 256);
	*(strrchr(moduledirectory, '\\') + 1) = 0;
	ModuleDirectory = moduledirectory;
	delete[] moduledirectory;
}

TLELoader tleload;

// ----------------------------------------------------------------------------------------------------------------------------

CTexture::CTexture()
{
	Texture = 0;
}

CTexture::~CTexture()
{
}

CTexture::operator GLuint ()
{
	return Texture;
}

bool CTexture::LoadTexture2D(char *FileName)
{
	CString DirectoryFileName = ModuleDirectory + FileName;

	int Width, Height, BPP;

	FIBITMAP *dib = GetBitmap(DirectoryFileName, Width, Height, BPP);

	if (dib == NULL)
	{
		ErrorLog.Append("Error loading texture " + DirectoryFileName + "!\r\n");
		return false;
	}

	GLenum Format = 0;

	if (BPP == 32) Format = GL_BGRA;
	if (BPP == 24) Format = GL_BGR;
	if (BPP == 8) Format = GL_LUMINANCE;

	if (Format == 0)
	{
		ErrorLog.Append("Unsupported texture format (%s)!\r\n", FileName);
		FreeImage_Unload(dib);
		return false;
	}

	Destroy();

	glGenTextures(1, &Texture);

	glBindTexture(GL_TEXTURE_2D, Texture);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

	if (GLEW_EXT_texture_filter_anisotropic)
	{
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_ANISOTROPY_EXT, gl_max_texture_max_anisotropy_ext);
	}

	glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE);

	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, Width, Height, 0, Format, GL_UNSIGNED_BYTE, FreeImage_GetBits(dib));

	glBindTexture(GL_TEXTURE_2D, 0);

	FreeImage_Unload(dib);

	return true;
}

bool CTexture::LoadTextureCubeMap(char **FileNames)
{
	int Width, Height, BPP;

	FIBITMAP *dib[6];

	bool Error = false;

	for (int i = 0; i < 6; i++)
	{
		CString DirectoryFileName = ModuleDirectory + FileNames[i];

		dib[i] = GetBitmap(DirectoryFileName, Width, Height, BPP);

		if (dib[i] == NULL)
		{
			ErrorLog.Append("Error loading texture " + DirectoryFileName + "!\r\n");
			Error = true;
		}
	}

	if (Error)
	{
		for (int i = 0; i < 6; i++)
		{
			FreeImage_Unload(dib[i]);
		}

		return false;
	}

	GLenum Format = 0;

	if (BPP == 32) Format = GL_BGRA;
	if (BPP == 24) Format = GL_BGR;

	if (Format == 0)
	{
		ErrorLog.Append("Unsupported texture format (%s)!\r\n", FileNames[5]);

		for (int i = 0; i < 6; i++)
		{
			FreeImage_Unload(dib[i]);
		}

		return false;
	}

	Destroy();

	glGenTextures(1, &Texture);

	glBindTexture(GL_TEXTURE_CUBE_MAP, Texture);

	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

	if (GLEW_EXT_texture_filter_anisotropic)
	{
		glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAX_ANISOTROPY_EXT, gl_max_texture_max_anisotropy_ext);
	}

	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_GENERATE_MIPMAP, GL_TRUE);

	for (int i = 0; i < 6; i++)
	{
		glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X + i, 0, GL_RGBA8, Width, Height, 0, Format, GL_UNSIGNED_BYTE, FreeImage_GetBits(dib[i]));
	}

	glBindTexture(GL_TEXTURE_CUBE_MAP, 0);

	for (int i = 0; i < 6; i++)
	{
		FreeImage_Unload(dib[i]);
	}

	return true;
}

void CTexture::Destroy()
{
	glDeleteTextures(1, &Texture);
	Texture = 0;
}

FIBITMAP *CTexture::GetBitmap(char *FileName, int &Width, int &Height, int &BPP)
{
	FREE_IMAGE_FORMAT fif = FreeImage_GetFileType(FileName);

	if (fif == FIF_UNKNOWN)
	{
		fif = FreeImage_GetFIFFromFilename(FileName);
	}

	if (fif == FIF_UNKNOWN)
	{
		return NULL;
	}

	FIBITMAP *dib = NULL;

	if (FreeImage_FIFSupportsReading(fif))
	{
		dib = FreeImage_Load(fif, FileName);
	}

	if (dib != NULL)
	{
		int OriginalWidth = FreeImage_GetWidth(dib);
		int OriginalHeight = FreeImage_GetHeight(dib);

		Width = OriginalWidth;
		Height = OriginalHeight;

		if (Width == 0 || Height == 0)
		{
			FreeImage_Unload(dib);
			return NULL;
		}

		BPP = FreeImage_GetBPP(dib);

		if (Width > gl_max_texture_size) Width = gl_max_texture_size;
		if (Height > gl_max_texture_size) Height = gl_max_texture_size;

		if (!GLEW_ARB_texture_non_power_of_two)
		{
			Width = 1 << (int)floor((log((float)Width) / log(2.0f)) + 0.5f);
			Height = 1 << (int)floor((log((float)Height) / log(2.0f)) + 0.5f);
		}

		if (Width != OriginalWidth || Height != OriginalHeight)
		{
			FIBITMAP *rdib = FreeImage_Rescale(dib, Width, Height, FILTER_BICUBIC);
			FreeImage_Unload(dib);
			dib = rdib;
		}
	}

	return dib;
}

// ----------------------------------------------------------------------------------------------------------------------------

CShaderProgram::CShaderProgram()
{
	SetDefaults();
}

CShaderProgram::~CShaderProgram()
{
}

CShaderProgram::operator GLuint ()
{
	return Program;
}

bool CShaderProgram::Load(char *VertexShaderFileName, char *FragmentShaderFileName)
{
	bool Error = false;

	Destroy();

	Error |= ((VertexShader = LoadShader(VertexShaderFileName, GL_VERTEX_SHADER)) == 0);
	Error |= ((FragmentShader = LoadShader(FragmentShaderFileName, GL_FRAGMENT_SHADER)) == 0);

	if (Error)
	{
		Destroy();
		return false;
	}

	Program = glCreateProgram();
	glAttachShader(Program, VertexShader);
	glAttachShader(Program, FragmentShader);
	glLinkProgram(Program);

	int LinkStatus;
	glGetProgramiv(Program, GL_LINK_STATUS, &LinkStatus);

	if (LinkStatus == GL_FALSE)
	{
		ErrorLog.Append("Error linking program (%s, %s)!\r\n", VertexShaderFileName, FragmentShaderFileName);

		int InfoLogLength = 0;
		glGetProgramiv(Program, GL_INFO_LOG_LENGTH, &InfoLogLength);

		if (InfoLogLength > 0)
		{
			char *InfoLog = new char[InfoLogLength];
			int CharsWritten = 0;
			glGetProgramInfoLog(Program, InfoLogLength, &CharsWritten, InfoLog);
			ErrorLog.Append(InfoLog);
			delete[] InfoLog;
		}

		Destroy();

		return false;
	}

	return true;
}

void CShaderProgram::Destroy()
{
	glDetachShader(Program, VertexShader);
	glDetachShader(Program, FragmentShader);

	glDeleteShader(VertexShader);
	glDeleteShader(FragmentShader);

	glDeleteProgram(Program);

	delete[] UniformLocations;
	delete[] AttribLocations;

	SetDefaults();
}

GLuint CShaderProgram::LoadShader(char *FileName, GLenum Type)
{
	CString DirectoryFileName = ModuleDirectory + FileName;

	FILE *File;

	if (fopen_s(&File, DirectoryFileName, "rb") != 0)
	{
		ErrorLog.Append("Error loading file " + DirectoryFileName + "!\r\n");
		return 0;
	}

	fseek(File, 0, SEEK_END);
	long Size = ftell(File);
	fseek(File, 0, SEEK_SET);
	char *Source = new char[Size + 1];
	fread(Source, 1, Size, File);
	fclose(File);
	Source[Size] = 0;

	GLuint Shader = glCreateShader(Type);

	glShaderSource(Shader, 1, (const char**)&Source, NULL);
	delete[] Source;
	glCompileShader(Shader);

	int CompileStatus;
	glGetShaderiv(Shader, GL_COMPILE_STATUS, &CompileStatus);

	if (CompileStatus == GL_FALSE)
	{
		ErrorLog.Append("Error compiling shader %s!\r\n", FileName);

		int InfoLogLength = 0;
		glGetShaderiv(Shader, GL_INFO_LOG_LENGTH, &InfoLogLength);

		if (InfoLogLength > 0)
		{
			char *InfoLog = new char[InfoLogLength];
			int CharsWritten = 0;
			glGetShaderInfoLog(Shader, InfoLogLength, &CharsWritten, InfoLog);
			ErrorLog.Append(InfoLog);
			delete[] InfoLog;
		}

		glDeleteShader(Shader);

		return 0;
	}

	return Shader;
}

void CShaderProgram::SetDefaults()
{
	VertexShader = 0;
	FragmentShader = 0;

	Program = 0;

	UniformLocations = NULL;
	AttribLocations = NULL;
}

// ----------------------------------------------------------------------------------------------------------------------------

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
	X = normalize(cross(vec3(0.0f, 1.0f, 0.0f), Z));
	Y = cross(Z, X);

	CalculateViewMatrix();
}

void CCamera::OnKeys(SHORT Keys, float FrameTime)
{
	float Speed = max(1.0f, min(PlanetRadius / 2.0f, length(Position) - PlanetRadius));

	if (Keys & 0x100) Speed *= 5.0f;
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

// ----------------------------------------------------------------------------------------------------------------------------

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
	//Error |= !LightsMap.LoadTexture2D("lightsmap.jpg");

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
		//glUniform1i(glGetUniformLocation(programs[i], "s2Tex3"), 2);
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

	float ainc = 0.25f * FrameTime;

	a += ainc;

	if (length(Camera.Position) < OuterRadius * 1.1f)
	{
		vec3 Y = vec3(0.0f, 1.0f, 0.0f);

		Camera.Position = rotate(Camera.Position, ainc, Y);

		Camera.X = rotate(Camera.X, ainc, Y);
		Camera.Y = rotate(Camera.Y, ainc, Y);
		Camera.Z = rotate(Camera.Z, ainc, Y);

		Camera.CalculateViewMatrix();
	}

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);

	GLuint program;

	glMatrixMode(GL_MODELVIEW);
	glLoadMatrixf(&ViewMatrix);
	glScalef(InnerRadius, InnerRadius, InnerRadius);
	glRotatef(a, 0.0f, 1.0f, 0.0f);

	glMatrixMode(GL_TEXTURE);
	glLoadIdentity();
	glScalef(InnerRadius, InnerRadius, InnerRadius);
	glRotatef(a, 0.0f, 1.0f, 0.0f);

	glActiveTexture(GL_TEXTURE0); glBindTexture(GL_TEXTURE_2D, EarthMap);
	glActiveTexture(GL_TEXTURE1); glBindTexture(GL_TEXTURE_2D, CloudsMap);
	//glActiveTexture(GL_TEXTURE2); glBindTexture(GL_TEXTURE_2D, LightsMap);

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

	//glActiveTexture(GL_TEXTURE2); glBindTexture(GL_TEXTURE_2D, 0);
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

	//glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
	//glEnable(GL_BLEND);

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
	//LightsMap.Destroy();

	glDeleteBuffers(1, &VerticesVBO);
	glDeleteBuffers(1, &TexCoordsVBO);

	gfs.Destroy();
	gfa.Destroy();
	sfs.Destroy();
	sfa.Destroy();
}

// ----------------------------------------------------------------------------------------------------------------------------

COpenGLView::COpenGLView()
{
}

COpenGLView::~COpenGLView()
{
}

bool COpenGLView::Init(HINSTANCE hInstance, char *Title, int Width, int Height, int Samples)
{
	this->Title = Title;
	this->Width = Width;
	this->Height = Height;

	WNDCLASSEX WndClassEx;

	memset(&WndClassEx, 0, sizeof(WNDCLASSEX));

	WndClassEx.cbSize = sizeof(WNDCLASSEX);
	WndClassEx.style = CS_OWNDC | CS_HREDRAW | CS_VREDRAW;
	WndClassEx.lpfnWndProc = WndProc;
	WndClassEx.hInstance = hInstance;
	WndClassEx.hIcon = LoadIcon(NULL, IDI_APPLICATION);
	WndClassEx.hIconSm = LoadIcon(NULL, IDI_APPLICATION);
	WndClassEx.hCursor = LoadCursor(NULL, IDC_ARROW);
	WndClassEx.lpszClassName = "Win32OpenGLWindowClass";

	if (RegisterClassEx(&WndClassEx) == 0)
	{
		ErrorLog.Set("RegisterClassEx failed!");
		return false;
	}

	DWORD Style = WS_OVERLAPPEDWINDOW | WS_CLIPSIBLINGS | WS_CLIPCHILDREN;

	hWnd = CreateWindowEx(WS_EX_APPWINDOW, WndClassEx.lpszClassName, Title, Style, 0, 0, Width, Height, NULL, NULL, hInstance, NULL);

	if (hWnd == NULL)
	{
		ErrorLog.Set("CreateWindowEx failed!");
		return false;
	}

	HDC hDC = GetDC(hWnd);

	if (hDC == NULL)
	{
		ErrorLog.Set("GetDC failed!");
		return false;
	}

	PIXELFORMATDESCRIPTOR pfd;

	memset(&pfd, 0, sizeof(PIXELFORMATDESCRIPTOR));

	pfd.nSize = sizeof(PIXELFORMATDESCRIPTOR);
	pfd.nVersion = 1;
	pfd.dwFlags = PFD_DRAW_TO_WINDOW | PFD_SUPPORT_OPENGL | PFD_DOUBLEBUFFER;
	pfd.iPixelType = PFD_TYPE_RGBA;
	pfd.cColorBits = 32;
	pfd.cDepthBits = 24;
	pfd.iLayerType = PFD_MAIN_PLANE;

	int PixelFormat = ChoosePixelFormat(hDC, &pfd);

	if (PixelFormat == 0)
	{
		ErrorLog.Set("ChoosePixelFormat failed!");
		return false;
	}

	static int MSAAPixelFormat = 0;

	if (SetPixelFormat(hDC, MSAAPixelFormat == 0 ? PixelFormat : MSAAPixelFormat, &pfd) == FALSE)
	{
		ErrorLog.Set("SetPixelFormat failed!");
		return false;
	}

	hGLRC = wglCreateContext(hDC);

	if (hGLRC == NULL)
	{
		ErrorLog.Set("wglCreateContext failed!");
		return false;
	}

	if (wglMakeCurrent(hDC, hGLRC) == FALSE)
	{
		ErrorLog.Set("wglMakeCurrent failed!");
		return false;
	}

	if (glewInit() != GLEW_OK)
	{
		ErrorLog.Set("glewInit failed!");
		return false;
	}

	if (!GLEW_VERSION_2_1)
	{
		ErrorLog.Set("OpenGL 2.1 not supported!");
		return false;
	}

	if (MSAAPixelFormat == 0 && Samples > 0)
	{
		if (GLEW_ARB_multisample && WGLEW_ARB_pixel_format)
		{
			while (Samples > 0)
			{
				UINT NumFormats = 0;

				int PFAttribs[] =
				{
					WGL_DRAW_TO_WINDOW_ARB, GL_TRUE,
					WGL_SUPPORT_OPENGL_ARB, GL_TRUE,
					WGL_DOUBLE_BUFFER_ARB, GL_TRUE,
					WGL_PIXEL_TYPE_ARB, WGL_TYPE_RGBA_ARB,
					WGL_COLOR_BITS_ARB, 32,
					WGL_DEPTH_BITS_ARB, 24,
					WGL_ACCELERATION_ARB, WGL_FULL_ACCELERATION_ARB,
					WGL_SAMPLE_BUFFERS_ARB, GL_TRUE,
					WGL_SAMPLES_ARB, Samples,
					0
				};

				if (wglChoosePixelFormatARB(hDC, PFAttribs, NULL, 1, &MSAAPixelFormat, &NumFormats) == TRUE && NumFormats > 0) break;

				Samples--;
			}

			wglDeleteContext(hGLRC);
			DestroyWindow(hWnd);
			UnregisterClass(WndClassEx.lpszClassName, hInstance);

			return Init(hInstance, Title, Width, Height, Samples);
		}
		else
		{
			Samples = 0;
		}
	}

	this->Samples = Samples;

	GetModuleDirectory();

	glGetIntegerv(GL_MAX_TEXTURE_SIZE, &gl_max_texture_size);

	if (GLEW_EXT_texture_filter_anisotropic)
	{
		glGetIntegerv(GL_MAX_TEXTURE_MAX_ANISOTROPY_EXT, &gl_max_texture_max_anisotropy_ext);
	}

	if (WGLEW_EXT_swap_control)
	{
		wglSwapIntervalEXT(0);
	}

	// Initialize AntTweakBar
	TwInit(TW_OPENGL, NULL);

	TwStructMember Vector3fMembers[] = {
		{ "x", TW_TYPE_FLOAT, offsetof(vec3, x), "" },
		{ "y", TW_TYPE_FLOAT, offsetof(vec3, y), "" },
		{ "z", TW_TYPE_FLOAT, offsetof(vec3, z), "" }
	};

	TW_TYPE_OGLDEV_VECTOR3F = TwDefineStruct("Vector3f", Vector3fMembers, 3, sizeof(vec3), NULL, NULL);

	
	TwStructMember Vector2iMembers[] = {
		{ "Norad id", TW_TYPE_FLOAT, offsetof(vec2i, x), "" },
		{ "Status", TW_TYPE_FLOAT, offsetof(vec2i, y), "" },
	};

	TW_TYPE_OGLDEV_VECTOR2I = TwDefineStruct("Vector2i", Vector2iMembers, 2, sizeof(vec2i), NULL, NULL);


	//TwStructMember AttenuationMembers[] = {
	//	{ "Const", TW_TYPE_FLOAT, offsetof(LightAttenuation, Constant), "" },
	//	{ "Linear", TW_TYPE_FLOAT, offsetof(LightAttenuation, Linear), "" },
	//	{ "Exp", TW_TYPE_FLOAT, offsetof(LightAttenuation, Exp), "" }
	//};
	//TW_TYPE_OGLDEV_ATTENUATION = TwDefineStruct("Attenuation", AttenuationMembers, 3, sizeof(LightAttenuation), NULL, NULL);

	tw_Zoom = 1.0f;
	tw_Fps = 0.0f;
	tw_Dist = 0.0f;
	tw_R = 0.0f;
	tw_pos = vec3();
	tw_dir = vec3();
	globalTime = 0;

	bar = TwNewBar("OGLDEV");
	TwDefine(" GLOBAL help='This example shows how to integrate AntTweakBar with GLUT and OpenGL.' "); // Message added to the help bar.
	TwDefine(" TweakBar size='200 400' color='96 216 224' "); // change default tweak bar size and color

	// Add 'g_Zoom' to 'bar': this is a modifable (RW) variable of type TW_TYPE_FLOAT. Its key shortcuts are [z] and [Z].
	TwAddVarRW(bar, "Zoom", TW_TYPE_FLOAT, &tw_Zoom, " min=0.01 max=2.5 step=0.01 keyIncr=z keyDecr=Z help='Scale the object (1=original size).' ");

	TwAddVarRW(bar, "FPS", TW_TYPE_FLOAT, &tw_Fps, " min=0.01 max=2.5 step=0.01 keyIncr=z keyDecr=Z help='Scale the object (1=original size).' ");

	TwAddVarRW(bar, "Dist", TW_TYPE_FLOAT, &tw_Dist, " min=0.01 max=2.5 step=0.01 keyIncr=z keyDecr=Z help='Scale the object (1=original size).' ");

	TwAddVarRW(bar, "Earth R", TW_TYPE_FLOAT, &tw_R, " min=0.01 max=2.5 step=0.01 keyIncr=z keyDecr=Z help='Scale the object (1=original size).' ");

	TwAddButton(bar, "Camera", NULL, NULL, "");
	TwAddVarRW(bar, "Position", TW_TYPE_OGLDEV_VECTOR3F, (void*)&tw_pos, NULL);
	TwAddVarRO(bar, "Direction", TW_TYPE_DIR3F, &tw_dir, " axisz=-z ");


	// load tle
	tleload.LoadData("TLE20180724.txt", 2);
	//std::string content = requestHttps.getResponse();
	//bool success = parseJsonResponse.ParseTelescopeStatus(status, content);
	//if (success){
    //
	//	TwAddButton(bar, "Telescope status", NULL, NULL, "");
	//	for (int i = 0; i < status.size(); i++)
	//	{
	//		statusFlag[i].x = status[i].norad_id;
	//		statusFlag[i].y = status[i].activity_code;
	//		TwAddVarRW(bar, status[i].id.c_str(), TW_TYPE_OGLDEV_VECTOR2I, &statusFlag[i], " min=0.01 max=2.5 step=0.01 keyIncr=z keyDecr=Z help='Scale the object (1=original size).' ");
	//	}
	//}

	return OpenGLRenderer.Init();
}

void COpenGLView::Show(bool Maximized)
{
	RECT dRect, wRect, cRect;

	GetWindowRect(GetDesktopWindow(), &dRect);
	GetWindowRect(hWnd, &wRect);
	GetClientRect(hWnd, &cRect);

	wRect.right += Width - cRect.right;
	wRect.bottom += Height - cRect.bottom;
	wRect.right -= wRect.left;
	wRect.bottom -= wRect.top;
	wRect.left = dRect.right / 2 - wRect.right / 2;
	wRect.top = dRect.bottom / 2 - wRect.bottom / 2;

	MoveWindow(hWnd, wRect.left, wRect.top, wRect.right, wRect.bottom, FALSE);

	ShowWindow(hWnd, Maximized ? SW_SHOWMAXIMIZED : SW_SHOWNORMAL);
}

void COpenGLView::MessageLoop()
{
	MSG Msg;

	while (GetMessage(&Msg, NULL, 0, 0) > 0)
	{
		TranslateMessage(&Msg);
		DispatchMessage(&Msg);
	}
}

void COpenGLView::Destroy()
{
	if (GLEW_VERSION_2_1)
	{
		OpenGLRenderer.Destroy();
	}

	TwTerminate();

	wglDeleteContext(hGLRC);
	DestroyWindow(hWnd);
}

void COpenGLView::OnKeyDown(UINT Key)
{
	/*switch(Key)
	{
	case VK_F1:
	break;

	case VK_SPACE:
	break;
	}*/
}

void COpenGLView::OnLButtonDown(int X, int Y)
{
	LastX = X;
	LastY = Y;

	TwMouseButton(TW_MOUSE_PRESSED, TW_MOUSE_LEFT);
	TwMouseMotion(X, Y);
}

void COpenGLView::OnMouseMove(int X, int Y)
{
	if ((GetKeyState(VK_LBUTTON) & 0x80) || GetKeyState(VK_RBUTTON) & 0x80)
	{
		Camera.OnMouseMove(LastX - X, LastY - Y);

		LastX = X;
		LastY = Y;
	}
}

void COpenGLView::OnMouseWheel(short zDelta)
{
}

void COpenGLView::OnPaint()
{
	static DWORD LastFPSTime = GetTickCount(), LastFrameTime = LastFPSTime, FPS = 0;

	PAINTSTRUCT ps;

	HDC hDC = BeginPaint(hWnd, &ps);

	DWORD Time = GetTickCount();

	float FrameTime = (Time - LastFrameTime) * 0.001f;

	LastFrameTime = Time;

	if (Time - LastFPSTime > 1000)
	{
		CString Text = Title;

		if (OpenGLRenderer.Text[0] != 0)
		{
			Text.Append(" - " + OpenGLRenderer.Text);
		}

		Text.Append(" - %dx%d", Width, Height);
		Text.Append(", ATF %dx", gl_max_texture_max_anisotropy_ext);
		Text.Append(", MSAA %dx", Samples);
		Text.Append(", FPS: %d", FPS);
		Text.Append(" - %s", glGetString(GL_RENDERER));

		SetWindowText(hWnd, Text);

		tw_Fps = FPS;

		LastFPSTime = Time;
		FPS = 0;
	}
	else
	{
		FPS++;
	}

	SHORT Keys = 0x000;

	if (GetKeyState('W') & 0x80) Keys |= 0x001;
	if (GetKeyState('S') & 0x80) Keys |= 0x002;
	if (GetKeyState('A') & 0x80) Keys |= 0x004;
	if (GetKeyState('D') & 0x80) Keys |= 0x008;
	if (GetKeyState('R') & 0x80) Keys |= 0x010;
	if (GetKeyState('F') & 0x80) Keys |= 0x020;
	if (GetKeyState('Q') & 0x80) Keys |= 0x040;
	if (GetKeyState('E') & 0x80) Keys |= 0x080;

	if (GetKeyState(VK_SHIFT) & 0x80) Keys |= 0x100;
	if (GetKeyState(VK_CONTROL) & 0x80) Keys |= 0x200;

	if (Keys & 0xFF)
	{
		Camera.OnKeys(Keys, FrameTime);
	}

	OpenGLRenderer.Render(FrameTime);

	//----------------------------------------------------------------------------------------------------------------------
	glColor3f(1.0f, 0.0f, 0.0f);
	glBegin(GL_LINES);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(2.0f*Camera.PlanetRadius, 0.0f, 0.0f);
	glEnd();

	glColor3f(0.0f, 1.0f, 0.0f);
	glBegin(GL_LINES);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(0.0f, 2.0f*Camera.PlanetRadius, 0.0f);
	glEnd();

	glColor3f(0.0f, 0.0f, 1.0f);
	glBegin(GL_LINES);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(0.0f, 0.0f, 2.0f*Camera.PlanetRadius);
	glEnd();

	glColor3f(1.0f, 0.0f, 0.0f);
	glPointSize(3.0f);

	int bandValue = 0;
	vec3 pos;

	for (int i = 0; i < tleload.NORADList.size(); i++)
	{
		cEciTime eci1 = tleload.NORADList[i]->GetPosition(globalTime*3.0f);

		
		vec3 pos2;
		// прогноз положения в километрах
		pos.x = eci1.Position().m_x;
		pos.y = eci1.Position().m_y;
		pos.z = eci1.Position().m_z;

		pos2.x = 4.0*pos.x / 6356000.0f * Camera.PlanetRadius;
		pos2.y = 4.0*pos.y / 6356000.0f * Camera.PlanetRadius;
		pos2.z = 4.0*pos.z / 6356000.0f * Camera.PlanetRadius;

		if (length(pos) < 6371.0f + 100.0f)
			bandValue++;

		glBegin(GL_POINTS);
			//glVertex3f(0.0f, 0.0f, 0.0f);
			glVertex3f(pos2.x, pos2.y, pos2.z );
		glEnd();
	}

	//----------------------------------------------------------------------------------------------------------------------

	tw_Zoom = bandValue;

	globalTime += FrameTime;

	OpenGLRenderer.RenderFinish();

	tw_R = Camera.PlanetRadius;
	tw_Dist = length(Camera.Position);

	tw_pos = pos;// Camera.Position;
	tw_dir = Camera.Z;

	TwDraw();

	SwapBuffers(hDC);

	EndPaint(hWnd, &ps);

	InvalidateRect(hWnd, NULL, FALSE);
}

void COpenGLView::OnRButtonDown(int X, int Y)
{
	LastX = X;
	LastY = Y;
}

void COpenGLView::OnSize(int Width, int Height)
{
	this->Width = Width;
	this->Height = Height;

	OpenGLRenderer.Resize(Width, Height);

	// Send the new window size to AntTweakBar
	TwWindowSize(Width, Height);
}

// ----------------------------------------------------------------------------------------------------------------------------

LRESULT CALLBACK WndProc(HWND hWnd, UINT uiMsg, WPARAM wParam, LPARAM lParam)
{
	switch (uiMsg)
	{
	case WM_CLOSE:
		PostQuitMessage(0);
		break;

	case WM_KEYDOWN:
		OpenGLView.OnKeyDown((UINT)wParam);
		break;

	case WM_LBUTTONDOWN:
		OpenGLView.OnLButtonDown(LOWORD(lParam), HIWORD(lParam));
		break;

	case WM_MOUSEMOVE:
		OpenGLView.OnMouseMove(LOWORD(lParam), HIWORD(lParam));
		break;

	case 0x020A: // WM_MOUSWHEEL
		OpenGLView.OnMouseWheel(HIWORD(wParam));
		break;

	case WM_PAINT:
		OpenGLView.OnPaint();
		break;

	case WM_RBUTTONDOWN:
		OpenGLView.OnRButtonDown(LOWORD(lParam), HIWORD(lParam));
		break;

	case WM_SIZE:
		OpenGLView.OnSize(LOWORD(lParam), HIWORD(lParam));
		break;

	default:
		return DefWindowProc(hWnd, uiMsg, wParam, lParam);
	}

	return 0;
}

// ----------------------------------------------------------------------------------------------------------------------------

int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR sCmdLine, int iShow)
{
	char *AppName = "Planet Earth from space";

	if (OpenGLView.Init(hInstance, AppName, 800, 600, 4))
	{
		OpenGLView.Show();
		OpenGLView.MessageLoop();
	}
	else
	{
		MessageBox(NULL, ErrorLog, AppName, MB_OK | MB_ICONERROR);
	}

	OpenGLView.Destroy();

	return 0;
}
