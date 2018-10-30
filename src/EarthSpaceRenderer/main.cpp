#include "main.h"
#include "es_globals.h"

#pragma comment(lib, "opengl32.lib")
#pragma comment(lib, "glew32.lib")
#pragma comment(lib, "FreeImage.lib")
#pragma comment(lib, "AntTweakBar.lib")

void GetModuleDirectory()
{
	char *moduledirectory = new char[256];
	GetModuleFileName(GetModuleHandle(NULL), moduledirectory, 256);
	*(strrchr(moduledirectory, '\\') + 1) = 0;
	ModuleDirectory = moduledirectory;
	delete[] moduledirectory;
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
