#include "es_view.h"
#include "es_globals.h"

LRESULT CALLBACK WndProc(HWND hWnd, UINT uiMsg, WPARAM wParam, LPARAM lParam);

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
		ErrorLog.Append("RegisterClassEx failed!\r\n");
		return false;
	}

	DWORD Style = WS_OVERLAPPEDWINDOW | WS_CLIPSIBLINGS | WS_CLIPCHILDREN;

	hWnd = CreateWindowEx(WS_EX_APPWINDOW, WndClassEx.lpszClassName, Title, Style, 0, 0, Width, Height, NULL, NULL, hInstance, NULL);

	if (hWnd == NULL)
	{
		ErrorLog.Append("CreateWindowEx failed!\r\n");
		return false;
	}

	HDC hDC = GetDC(hWnd);

	if (hDC == NULL)
	{
		ErrorLog.Append("GetDC failed!\r\n");
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
		ErrorLog.Append("ChoosePixelFormat failed!\r\n");
		return false;
	}

	static int MSAAPixelFormat = 0;

	if (SetPixelFormat(hDC, MSAAPixelFormat == 0 ? PixelFormat : MSAAPixelFormat, &pfd) == FALSE)
	{
		ErrorLog.Append("SetPixelFormat failed!\r\n");
		return false;
	}

	hGLRC = wglCreateContext(hDC);

	if (hGLRC == NULL)
	{
		ErrorLog.Append("wglCreateContext failed!\r\n");
		return false;
	}

	if (wglMakeCurrent(hDC, hGLRC) == FALSE)
	{
		ErrorLog.Append("wglMakeCurrent failed!\r\n");
		return false;
	}

	if (glewInit() != GLEW_OK)
	{
		ErrorLog.Append("glewInit failed!\r\n");
		return false;
	}

	if (!GLEW_VERSION_2_1)
	{
		ErrorLog.Append("OpenGL 2.1 not supported!\r\n");
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
	tleload.LoadData("data/TLE20180724.txt", 2);

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
		glVertex3f(pos2.x, pos2.y, pos2.z);
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
