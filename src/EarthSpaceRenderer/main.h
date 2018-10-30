// ----------------------------------------------------------------------------------------------------------------------------
//
// Version 2.02
//
// ----------------------------------------------------------------------------------------------------------------------------

#include <windows.h>

#include "glmath.h"
#include "string.h"

#include <gl/glew.h> // http://glew.sourceforge.net/
#include <gl/wglew.h>

#include <FreeImage.h> // http://freeimage.sourceforge.net/

#include <AntTweakBar.h> //https://triplepointfive.github.io/ogltutor/tutorials/tutorial48.html

//#include "request/requestData.h"
// ----------------------------------------------------------------------------------------------------------------------------

#pragma comment(lib, "opengl32.lib")
#pragma comment(lib, "glu32.lib")
#pragma comment(lib, "glew32.lib")
#pragma comment(lib, "FreeImage.lib")
#pragma comment(lib, "AntTweakBar.lib")
//#pragma comment(lib, "libcurl.lib")

// ----------------------------------------------------------------------------------------------------------------------------

extern CString ModuleDirectory, ErrorLog;

// ----------------------------------------------------------------------------------------------------------------------------

#define BUFFER_SIZE_INCREMENT 1048576

// ----------------------------------------------------------------------------------------------------------------------------

class CBuffer
{
private:
	BYTE *Buffer;
	int BufferSize, Position;

public:
	CBuffer();
	~CBuffer();

	void AddData(void *Data, int DataSize);
	void Empty();
	void *GetData();
	int GetDataSize();

private:
	void SetDefaults();
};

// ----------------------------------------------------------------------------------------------------------------------------

extern int gl_max_texture_size, gl_max_texture_max_anisotropy_ext;

// ----------------------------------------------------------------------------------------------------------------------------

class CTexture
{
protected:
	GLuint Texture;

public:
	CTexture();
	~CTexture();

	operator GLuint ();

	bool LoadTexture2D(char *FileName);
	bool LoadTextureCubeMap(char **FileNames);
	void Destroy();

protected:
	FIBITMAP *CTexture::GetBitmap(char *FileName, int &Width, int &Height, int &BPP);
};

// ----------------------------------------------------------------------------------------------------------------------------

class CShaderProgram
{
protected:
	GLuint VertexShader, FragmentShader, Program;

public:
	GLuint *UniformLocations, *AttribLocations;

public:
	CShaderProgram();
	~CShaderProgram();

	operator GLuint ();

	bool Load(char *VertexShaderFileName, char *FragmentShaderFileName);
	void Destroy();

protected:
	GLuint LoadShader(char *FileName, GLenum Type);
	void SetDefaults();
};

// ----------------------------------------------------------------------------------------------------------------------------

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

// ----------------------------------------------------------------------------------------------------------------------------

class COpenGLRenderer
{
protected:
	int Width, Height;
	mat3x3 NormalMatrix;
	mat4x4 ModelMatrix, ViewMatrix, ProjectionMatrix;

protected:
	CTexture EarthMap, CloudsMap;//, LightsMap;
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

// ----------------------------------------------------------------------------------------------------------------------------
#define MAX_TELESCOPES 32

class COpenGLView
{
protected:
	char *Title;
	int Width, Height, Samples;
	HWND hWnd;
	HGLRC hGLRC;

	TwBar *bar;
	float tw_Zoom;
	float tw_Fps;
	float tw_Dist;
	float tw_R;
	vec3 tw_pos, tw_dir;
	float globalTime;

	//RequestHttps requestHttps;
	//ParseJsonResponse parseJsonResponse;
	//std::vector<telescopeStatus> status;
	vec2i statusFlag[MAX_TELESCOPES];

	// ������ �����������
	//Force::InfluenceForce *IForce;

protected:
	int LastX, LastY;

public:
	COpenGLView();
	~COpenGLView();

	bool Init(HINSTANCE hInstance, char *Title, int Width, int Height, int Samples);
	void Show(bool Maximized = false);
	void MessageLoop();
	void Destroy();

	void OnKeyDown(UINT Key);
	void OnLButtonDown(int X, int Y);
	void OnMouseMove(int X, int Y);
	void OnMouseWheel(short zDelta);
	void OnPaint();
	void OnRButtonDown(int X, int Y);
	void OnSize(int Width, int Height);
};

// ----------------------------------------------------------------------------------------------------------------------------

LRESULT CALLBACK WndProc(HWND hWnd, UINT uiMsg, WPARAM wParam, LPARAM lParam);
int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR sCmdLine, int iShow);
