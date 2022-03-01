#pragma once

class CTimeModel
{
public:
    //! \returns the amount of seconds of the real time passed from the first call of this function
    static double getTime();

    /*! Updates the current value of the model time
        \returns the amount of the real time passed from the start in seconds
    */
    double Update() {
        return (t = getTime());
    }

    //! \returns the amount of the real time passed from the start in seconds
    double time() const {
        return t;
    }

    //! \returns the current value of the model time in minutes
    float now() const {
        return (float)(m*t);
    }

    /*! Creates new instance of CTimeModel
        \param[in] mult The ratio between the model time and real time
    */
    CTimeModel(double mult = 60.) : m(mult / 60.), t(0.) {}

    //! Sets the ratio between the model time and real time
    void setTimeMultiplier(double mult) {
        m = mult / 60.;
    }

private:
    double m, t;
};
