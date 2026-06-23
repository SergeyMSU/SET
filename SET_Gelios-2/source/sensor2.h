#pragma once

using namespace std;
class sensor2
{
public:
    sensor2(unsigned int a1, unsigned int a2, unsigned int a3, unsigned int a4);

    // Generate random number by this sensor.
    double MakeRandom();
    unsigned int taus_rand_step(unsigned int& state, int S1, int S2, int S3, int M);
    void Restart();
    unsigned int a1_;
    unsigned int a2_;
    unsigned int a3_;
    unsigned int a4_;
    long unsigned int num;
};

