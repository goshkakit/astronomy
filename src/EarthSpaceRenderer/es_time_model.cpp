#include "es_time_model.h"

#include <Windows.h>

double CTimeModel::getTime()
{
#ifdef WIN64
    ULONGLONG v_ticks_at_now = GetTickCount64();
    static ULONGLONG v_ticks_at_start = v_ticks_at_now;
#else
    DWORD v_ticks_at_now = GetTickCount();
    static DWORD v_ticks_at_start = v_ticks_at_now;
#endif

    return (double)(v_ticks_at_now - v_ticks_at_start) / 1000.;
}
