#pragma once

#include "es_defines.h"
#include "es_camera.h"
#include "es_renderer.h"
#include "es_view.h"
#include "es_time_model.h"
#include "es_earth_model.h"

extern TwType TW_TYPE_OGLDEV_VECTOR3F;
extern TwType TW_TYPE_OGLDEV_VECTOR2I;
extern TwType TW_TYPE_OGLDEV_ATTENUATION;

extern int gl_max_texture_size, gl_max_texture_max_anisotropy_ext;

extern CCamera Camera;

extern COpenGLRenderer OpenGLRenderer;

extern CString ModuleDirectory, DataDirectory, ErrorLog;

extern COpenGLView OpenGLView;

extern CTimeModel TimeModel;

extern CEarthModel EarthModel;

extern volatile int bCameraRotationEnabled;
