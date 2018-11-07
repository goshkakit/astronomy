
#include "EarthSpaceRenderer/es_globals.h"
#include "EarthSpaceRenderer/es_scene_norad.h"

#include <Windows.h>

int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR sCmdLine, int iShow)
{
	char *AppName = "Planet Earth from space";

	if (OpenGLView.Init(hInstance, AppName, 800, 600, 4))
	{
		Scene::CScene scene;
		scene.layers.push_back(dynamic_cast<Scene::CLayer *>(
			new Scene::CNoradLayer("data/test/TLE20180724.txt", 2)
			)
		);
		OpenGLView.setScene(&scene);

		OpenGLView.Show();
		OpenGLView.MessageLoop();

		scene.FreeInstances();
	}
	else
	{
		MessageBox(NULL, ErrorLog, AppName, MB_OK | MB_ICONERROR);
	}

	OpenGLView.Destroy();

	return 0;
}
