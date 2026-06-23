#pragma once

#include "Help.h"
#include <vector>

using namespace std;

enum Rail_type
{
    // Ниже находятся перечислители - все возможные значения этого типа данных
    // Каждый перечислитель отделяется запятой (НЕ точкой с запятой)
    A,
    B,
    C,
    D,
    No,
};

class Point;

class Rail
{
public:
    Rail_type type;
    vector <Point*> All_point;
    vector <Point*> Key_point;
    int M1;
    int M2;
    int M3;
    int M4;

    double s; // Угол луча или что-то другое для другого rail_type

	Rail(const double& s);

    void Init_start(Rail* T);
};
