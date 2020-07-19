#pragma once

#include "cearth.h"
#include "glmath.h"

class CEarthModel
{
public:
    /*! Creates new instance of CEarthModel
        \param[in] epoch_start JD of the model time start
    */
    CEarthModel(double epoch_start = 2400000.);
    /*! Creates new instance of CEarthModel
        \param[in] epoch_start long JD of the model time start
    */
    CEarthModel(const Time::long_jd &epoch_start);
    //! Frees the instance
    ~CEarthModel();

    //! Sets another JD of the model time start
    void setEpochStart(double epoch_start) {
        Earth_at_start.setTime(epoch_start);
    }
    //! Sets another long JD of the model time start
    void setEpochStart(const Time::long_jd &epoch_start) {
        Earth_at_start.setTime(epoch_start);
    }
    //! Model time start
    const Time::long_jd & epochStart() const {
        return Earth_at_start.ljd();
    }

    //! Earth state at the model time start
    const Space::CEarth & EarthAtStart() const {
        return Earth_at_start;
    }
    //! Earth state at the current model time
    const Space::CEarth & EarthAtNow() const {
        return Earth;
    }
    //! Earth state at the chosen long JD
    static Space::CEarth getEarthAtTime(const Time::long_jd &ljd) {
        return Space::CEarth(ljd);
    }
    //! Earth state at the chosen model time
    Space::CEarth EarthAtTime(float time) const {
        return Space::CEarth(Time::long_jd(epochStart()).addSeconds(60.*time));
    }

    /*! Updates Earth state to the chosen model time
        \returns Rotation matrix from the previous to the chosen model time
    */
    Space::m3x3d Update(float time);

    /*! Rotates GL vector by the rotation matrix
        \param[in] v GL vector at the previous model time in the Scene coordinate system (YZX)
        \param[in] m Rotation matrix from the previous to the chosen model time in the Space coordinate system (XYZ)
        \returns new value of the given GL vector at the chosen model time in the Scene coordinate system (YZX)
    */
    static vec3 rotateVec3(const vec3 &v, const Space::m3x3d &m);

    /*! Transforms Orientation matrix
        \param[in] m Orientation matrix in the Space coordinate system (XYZ)
        \returns GL 4x4 Transformation matrix in the Scene coordinate system (YZX)
    */
    static mat4x4 getMat4x4(const Space::m3x3d &m);

private:
    Space::CEarth Earth_at_start;
    Space::CEarth Earth;
};
