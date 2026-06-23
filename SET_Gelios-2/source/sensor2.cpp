#include "sensor2.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <math.h>
#include <fstream>

using namespace std;

sensor2::sensor2(unsigned int a1, unsigned int a2, unsigned int a3, unsigned int a4) :
    a1_(a1), a2_(a2), a3_(a3), a4_(a4)
{
    this->num = 1;
}

unsigned int sensor2::taus_rand_step(unsigned int& state, int S1, int S2, int S3, int M)
{
    unsigned int b = (((state << S1) ^ state) >> S2);
    state = (((state & M) << S3) ^ b);
    return state;
}

double sensor2::MakeRandom() {
    this->num++;
    auto l = (taus_rand_step(a1_, 6, 13, 18, 4294967294UL) ^
        taus_rand_step(a2_, 2, 27, 2, 4294967288UL) ^
        taus_rand_step(a3_, 13, 21, 7, 4294967280UL) ^
        taus_rand_step(a4_, 3, 12, 13, 4294967268UL));
    /*if (l >= 4294967297UL)
    {
        cout << l << " " << 4294967297UL << endl;
        cout << " 1 !!!!" << endl;
        exit(-1);
    }*/
    //double k = (1.0 * (l + 1)) / 4294967297.0;
    double k = (1.0 * l + 1.0) / 4294967297.0;

    //cout << l << " " << 4294967295UL << endl;
    /*if (fpclassify(log(k)) != FP_NORMAL && fpclassify(log(k)) != FP_ZERO)
    {
        cout << "ERROR " << endl;
        cout << l << " " << k << " " << log(k) << endl;
        cout << "log(k)" << endl;
        exit(-1);
    }

    if (fpclassify(log(1.0 - k)) != FP_NORMAL && fpclassify(log(1.0 - k)) != FP_ZERO)
    {
        cout << "ERROR " << endl;
        cout << l << " " << k << " " << log(k) << endl;
        cout << "log(1 - k)" << endl;
        exit(-1);
    }*/

    return k;
}