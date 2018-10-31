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
