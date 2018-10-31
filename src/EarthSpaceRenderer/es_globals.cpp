#include "es_globals.h"

TwType TW_TYPE_OGLDEV_VECTOR3F;
TwType TW_TYPE_OGLDEV_VECTOR2I;
TwType TW_TYPE_OGLDEV_ATTENUATION;

int gl_max_texture_size = 0, gl_max_texture_max_anisotropy_ext = 0;

CCamera Camera;

COpenGLRenderer OpenGLRenderer;

CString ModuleDirectory = "./", DataDirectory = "../data/render/", ErrorLog;

COpenGLView OpenGLView;

TLELoader tleload;

#ifdef _DEBUG
#include <stdio.h>
#include <time.h>

class CErrorLogSaver
{
public:
	CErrorLogSaver();
	~CErrorLogSaver();
};

static CErrorLogSaver ErrorLogSaver;

CErrorLogSaver::CErrorLogSaver()
{

}

CErrorLogSaver::~CErrorLogSaver()
{
	FILE *f = fopen((char*)(ModuleDirectory + "EarthSpaceRenderer.log"), "ab");

	if (f)
	{
		fprintf(f, "\n");
		fprintf(f, "--------------------------------------------------------------------------------\n");

		time_t rawtime;
		struct tm * timeinfo;
		time(&rawtime);
		timeinfo = localtime(&rawtime);
		fprintf(f, "%s", asctime(timeinfo));
		fprintf(f, "--------------------------------------------------------------------------------\n");

		char * s = (char*)ErrorLog;
		fwrite(s, strlen(s), 1, f);

		fclose(f);
	}
}
#endif//_DEBUG
