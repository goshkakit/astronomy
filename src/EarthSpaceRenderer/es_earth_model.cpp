#include "es_earth_model.h"

CEarthModel::CEarthModel(double epoch_start /*= 2400000.*/)
    : Earth_at_start(epoch_start)
    , Earth(epoch_start)
{

}

CEarthModel::CEarthModel(const Time::long_jd &epoch_start)
    : Earth_at_start(epoch_start)
    , Earth(epoch_start)
{

}

CEarthModel::~CEarthModel()
{

}

Space::m3x3d CEarthModel::Update(float time)
{
    Space::m3x3d m1 = Earth.orientation();

    Earth.setTime(Time::long_jd(epochStart()).addSeconds(60.*time));
    Space::m3x3d m2 = Earth.rotator();

    return (m2*m1);
}

vec3 CEarthModel::rotateVec3(const vec3 &v, const Space::m3x3d &m)
{
    Space::vec3d r = m * Space::vec3d(v.z, v.x, v.y);

    return vec3((float)r.y, (float)r.z, (float)r.x);
}

mat4x4 CEarthModel::getMat4x4(const Space::m3x3d &m)
{
    mat4x4 R;

    R[0] = (float)m.yy; R[1] = (float)m.yz; R[2] = (float)m.yx;
    R[4] = (float)m.zy; R[5] = (float)m.zz; R[6] = (float)m.zx;
    R[8] = (float)m.xy; R[9] = (float)m.xz; R[10] = (float)m.xx;
    R[15] = 1.f;

    return R;
}
