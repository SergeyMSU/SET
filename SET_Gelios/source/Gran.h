#pragma once
#include "Help.h"
#include <vector>

using namespace std;

enum Gran_type  // Тип грани нужен для граничных условий
{
    Usualy,                // Обычные
    Inner_sphere,          // Внутренняя сфера
    Extern,                // Выходные условия
    Input,                 // Входной поток
    Axis,                  // Ось симметрии
    Upper_wall,            // Верхняя стенка
    G_no,
};

class Point;
class Cell;
struct Parametr;

class Gran   // Нормаль к грани смотрит направо, если идти от узла A к B. Нормаль должна указывать на соседнюю ячейку.
{
public:
    int number;
    Gran_type type;
    Point* A;
    Point* B;
    Cell* Master;        // Ячейка, которой принадлежит грань.
    Cell* Sosed;         // Ячейка - сосед по этой грани
    Cell* Sosed_down;    // Это для TVD
    Cell* Sosed_up;      // Сосед для ТВД по направлению внешней нормали к грани
    Gran* Gran_copy;     // Грань - копия, но с другой нормалью. Введена для удобства
    double a;
    double b;
    double aa;           // aa * x + bb * y + cc = 0
    double bb;
    double cc;
    bool parallel;       // true если грань параллельна оси У
    bool main_gran;      // true - для основной грани, false для фантомной
    int koef;            // коеффициент для умножения уравнения прямой для определения с какой стороны точка  (y - ax - b = 0)
                         // должен быть меньше нуля
    int metod_HLLC;

    Gran(Point* A, Point* B, Gran_type type = G_no);
    void Get_Center(double& x, double& y);
    void Get_Center_posle(double& x, double& y);
    void Get_normal(double& n1, double& n2);
    void Get_normal_posle(double& n1, double& n2);
    double Get_square(void);
    double Get_lenght(void);
    double Get_square_rotate(const double& angle);

    void Get_par(Parametr& par, int i); // Здесь задаются граничные условия на грани
    void Get_par_TVD(Parametr& par, int i);
    void Get_par_TVD_radial(Parametr& par, int i);

    bool belong(const double& x, const double& y);         // принадлежит ли точка полупространству грани

    bool belong_gran(const double& x, const double& y);         // принадлежит ли точка грани

    void renew(void); // Обновить значения a и b
 

};

