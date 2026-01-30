#pragma once
#include <vector>
#include "Help.h"
#include <mutex>

using namespace std;

enum Cell_type 
{
	C_1,
	C_2,
	C_3,
	C_4,
	C_5,
	C_centr,
	C_no,
};

struct Parametr
{                    // ДОБАВЛЯТЬ ПЕРЕМЕННЫЕ СТРОГО В КОНЕЦ!!!!
	double ro = 0.0;
	double p = 0.0;
	double u = 0.0;
	double v = 0.0;
	double Q = 0.0;
	double ro_H1 = 0.0;
	double p_H1 = 0.0;
	double u_H1 = 0.0;
	double v_H1 = 0.0;
	double ro_H2 = 0.0;
	double p_H2 = 0.0;
	double u_H2 = 0.0;
	double v_H2 = 0.0;
	double ro_H3 = 0.0;
	double p_H3 = 0.0;
	double u_H3 = 0.0;
	double v_H3 = 0.0;
	double ro_H4 = 0.0;
	double p_H4 = 0.0;
	double u_H4 = 0.0;
	double v_H4 = 0.0;
	double F_n = 0.0;
	double F_u = 0.0;
	double F_v = 0.0;
	double F_T = 0.0;
	double I_u = 0.0;
	double I_v = 0.0;
	double I_T = 0.0;
	double II_u = 0.0;
	double II_v = 0.0;
	double II_T = 0.0;
	double M_u = 0.0;   // Мультифлюид источники
	double M_v = 0.0;
	double M_T = 0.0;
	double H_n[pop_atom];
	double H_u[pop_atom];
	double H_v[pop_atom];
	double H_T[pop_atom];
	double H_n2[4];
	double H_u2[4];
	double H_v2[4];
	double H_T2[4];
	double H_n3[4];
	double H_u3[4];
	double H_v3[4];
	double H_T3[4];
	double H_uu[pop_atom];
	double H_vv[pop_atom];
	double H_uv[pop_atom];
	double H_uuu[pop_atom];
	double k_u = 0.0;
	double k_v = 0.0;
	double k_T = 0.0;
	double k_u2 = 0.0;
	double k_v2 = 0.0;
	double k_T2 = 0.0;
	double k_u3 = 0.0;
	double k_v3 = 0.0;
	double k_T3 = 0.0;
	int num_atoms[4];   // Число перезарядок в ячейке
	double w_m[7];       // Средние веса по сортам
	double I1_mf[4];
	double I2_mf[4];
	double I3_mf[4];
	double I1_mc[4];
	double I2_mc[4];
	double I3_mc[4];

	/*double npui = 0.0;
	double Tpui = 0.0;
	double pp = 0.0;
	double ppui = 0.0;
	double divV = 0.0;*/
};

class Point;
class Gran;
class Sensor;
class Dist_func;

class Cell
{
public:
	double Potok[5];
	double Potok_H1[4];
	double Potok_H2[4];
	double Potok_H3[4];
	double Potok_H4[4];
	double L;                   // Характерный размер ячейки
	Parametr par[2];
	vector <Point*> contour;    // Гарантируется расположение точек по кругу
	vector <Gran*> Grans;


	double pogloshenie[4][pogl_rad_];

	int i_pui = 0;             // Количество сортов пуи
	vector <double> fpui;
	vector <double> Wmin_sort;
	vector <double> Wmax_sort;
	vector <double> fpui_max_sort;
	double fpui_max = 1.0;      // Максимальное значение fpui  
	double Wmin = 0.0;
	double Wmax = 1000.0;         // Максимальная скорость при которой fpui не нулевая
	vector <vector <double>> nu_pui;     // Частоты перезарядки от L, где L от 0 до 20 с шагом 100
	vector <vector <double>> nu2_pui;     // Источник импульса перезарядки от L, где L от 0 до 20 с шагом 100
	vector <vector <double>> nu3_pui;     // Источник импульса перезарядки от L, где L от 0 до 20 с шагом 100
	bool pui_ = false;                  // Нужно ли вообще считать эти PUI?

	double S_p[n_S];              // Для  S+  и  S-
	double S_m[n_S];

	Dist_func* df_s4;
	bool df_s4_bool = false;

	Dist_func* df_s3;
	bool df_s3_bool = false;

	Dist_func* df_s2;
	bool df_s2_bool = false;

	Dist_func* df_s1;
	bool df_s1_bool = false;

	int number;
	int zona;
	int zona_alpha;
	Cell_type type;
	mutex mut;                 // Мьютекс для записи в ячейку
	double x_min;
	double x_max;
	double y_min;
	double y_max;

	double r_istoch;           // Радиус источников (в какой точке они расположены, нужно для переинтерполяции)
	Cell* Back;                // ячейка до переинтерполяции
	Cell* Next;

    // Блок переменных чисто для монте-карло
	bool axis_; // Является ли ячейка граничной с осью (для сноса скорости в монте-карло
	double y_ax; // Высота центра граничной ячейки над осью

	double x_center;
	double y_center;
	double alf_center;

	Cell(Point* A, Point* B, Point* C, Point* D); 
	Cell(Point* A, Point* B, Point* C);
	Cell(void);

	void Initial(void);

	void Get_Center(double& x, double& y);
	void Get_Center_posle(double& x, double& y);
	void Get_Center2(double& x, double& y);
	void Get_Center_posle2(double& x, double& y);
	double Get_Volume(void);
	double Get_Volume_rotate(const double& angle);
	double Get_Volume_posle(void);
	double Get_Volume_posle_rotate(const double& angle);

	double get_nu_pui(const double& L, int i);
	double get_nu_pui2(const double& L, int i);
	double get_nu_pui3(const double& L, int i);
	double get_fpui(const double& W, const double& Wmin, const double& Wmax);
	bool Change_Velosity_PUI(Sensor* sens, const double& Vh1, const double& Vh2, const double& Vh3, //
		const double& Vp1, const double& Vp2, const double& Vp3, double& W1, double& W2, double& W3, int Nw, //
		const double& wmin, const double& wmax, int ifg);
	bool Change_Velosity_PUI2(Sensor* sens, const double& Vh1, const double& Vh2, const double& Vh3, //
		const double& Vp1, const double& Vp2, const double& Vp3, double& W1, double& W2, double& W3, int Nw, //
		const double& wmin, const double& wmax, int ifg);

	// Посчитать источники Монте-Карло
	void Get_Sourse_MK1(double& q1, double& q2, double& q3, const double& u, const double& v, const double& ro, const double& p, bool interpol = true);
	void Get_Sourse_MK2(double& q1, double& q2, double& q3, const double& u, const double& v, const double& ro, const double& p);
	void Get_Sourse_MK3(double& q1, double& q2, double& q3, const double& u, const double& v, const double& ro, const double& p);


	bool belong(const double& x, const double& y);

	void renew(void); // Обновить значения L

	void Calc_Sourse(void);  // Функция вычисляющая источники

};



