#pragma once

#include <vector>
#include <mutex>
/* Each sensor is determined by 3 integers for initial point. Structure of all sensors is the same. */
using namespace std;
const int SENSORS_AMOUNT = 270;

class Sensor {
public:
    inline Sensor(int a1, int a2, int a3);

    // Generate random number by this sensor.
    inline double MakeRandom();
    int a1_;
    int a2_;
    int a3_;
};

//std::vector<Sensor> InitSensors();

inline Sensor::Sensor(int a1, int a2, int a3) :
    a1_(a1), a2_(a2), a3_(a3)
{
}

inline double Sensor::MakeRandom() {
    int ic15 = 32768, ic10 = 1024;
    int mz = 710, my = 17784, mx = 11973;
    double xi = 9.0949470177292824E-13, c = 1.073741824E9;
    double b;
    int i13, i12, i11, ii;
    i13 = mz * a1_ + my * a2_ + mx * a3_;
    i12 = my * a1_ + mx * a2_;
    i11 = mx * a1_;
    ii = i11 / ic15;
    i12 = i12 + ii;
    a1_ = i11 - ic15 * ii;
    ii = i12 / ic15;
    i13 = i13 + ii;
    a2_ = i12 - ic15 * ii;
    a3_ = i13 % ic10;
    b = xi * (c * a3_ + ic15 * a2_ + a1_);
    return b;
}

