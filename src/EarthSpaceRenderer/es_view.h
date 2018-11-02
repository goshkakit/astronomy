#pragma once

#include "es_defines.h"
#include "es_scene.h"

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

	Scene::CScene *pScene;

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

	void setScene(Scene::CScene * scene);
};
