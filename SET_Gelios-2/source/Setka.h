#pragma once

#include "Help.h"
#include <vector>
#include <string>
#include <mutex>
#include "sensor2.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>

class Rail;
class Point;
class Cell;
class Gran;
class Sensor;
class MKmethod;
class Dist_func;
using namespace std;

class Setka
{
public:
	vector <Rail*> A_Rails;
	vector <Rail*> B_Rails;
	vector <Rail*> C_Rails;
	vector <Rail*> D_Rails;

	vector <Point*> All_Points;
	vector <Gran*> All_Gran;           // Грани основные
	vector <Gran*> All_Gran_copy;      // Грани фантомные (нормаль в другую сторону)

	vector <Cell*> All_Cells;          // Все ячейки
	vector <Cell*> All_Cells_Inner;    // Ячейки внутри маленького радиуса (там отдельно считает)
	vector <Cell*> All_Cells_zero;    // Ячейки в нуле

	vector <Gran*> Line_Contact;     // Контакт
	vector <Gran*> Line_Inner;       // Внутренняя волна
	vector <Gran*> Line_Outer;		 // Внешняя волна

	vector <Point*> Contact;	     // Контакт
	vector <Point*> Inner;		     // Внутренняя волна
	vector <Point*> Outer;		     // Внешняя волна

	vector <Cell*> Cell_sphere;      // Ячейки внешнии 
	vector <Cell*> Cell_side;      // Ячейки боковые
	vector <Cell*> Cell_back;      // Ячейки задние
	vector <Cell*> Cell_disk;      // Ячейки диска
	vector <Cell*> Cell_other;      // Ячейки , которые не считаются
	Cell* Cell_m;                  // Для 4-го типа вылета частиц

	vector <Dist_func*> Dist_func_all;

	vector<Sensor*> Sensors;
	vector<sensor2*> Sensors2;
	vector <double> Ri;               // Геометрические зоны для расщепления Монте-Карло
	double Mu[4][9];                  // Это и есть веса для каждой зоны (по радиусу) и каждого сорта атома (из четырёх)
	double Mu_stat[4][9];      // Считаем статистику для весов
	int I_stat[4][9];      // Считаем статистику для весов


	
	//double pogloshenie[4][pogl_alf_][pogl_rad_];
	//const double pogVmin = -10.0;
	//const double pogVmax = 5.0;
	//const double pogRmax = 2.0;

	double Mu_statistic[4][I_][J_];      // Считаем статистику для весов
	double Mu_[4][I_][J_];      // Применяем статистику для весов
	double SINKR[J_];            // Критические синусы для зон по углу по МОДУЛЮ
	mutex m_m;
	ofstream f_way;
	int f_num;

	mutex mut_1;
	int k_1;

	// Сбор статистики весов
	double mmu1;
	double mmu2;
	double mmu3;
	double mmu4;
	double mmu5;
	int mn1;
	int mn2;
	int mn3;
	int mn4;
	int mn5;



	int N1;
	int N2;
	int N3;
	int N4;
	int M1;
	int M2;
	int M3;
	int M4;


	int Number1;
	int Number2;
	int Number3;
	int Number4;
	int AllNumber;
	double sqv_1;
	double sqv_2;
	double sqv_3;
	double sqv_4;
	double sum_s;

	mutex Smut;

	double* V_r_stat;
	double* V_t_stat;
	double* V_p_stat;
	double* mu_stat;
	double* phi_stat;
	int* num_stat;
	int number_stat;
	int number_stat_2;
	mutex mut_stat;

	// для ПУИ
	double wmin = 0.0;
	double wmax = 0.0;
	int Nw = 0;


	// Моменты функции распределения для Насти
	double mu_mom[4][54];
	double Vx_mom[4][54];
	double Vy_mom[4][54];
	double Vxx_mom[4][54];
	double Vyy_mom[4][54];
	double Vxy_mom[4][54];
	double Vxxx_mom[4][54];
	double T_mom[4][54];
	mutex mut_mom;


	// Следующие массивы для BS, HP, TS для движения сетки, они используются (инициализируются и удаляются)
	// в функции Download_surface
	double* BS_al;
	double* BS_r;
	double* BS_x;
	double* BS_y;
	double* HP_al;
	double* HP_r;
	double* HP_x;
	double* HP_y;
	double* TS_al;
	double* TS_r;
	int TS_n;
	int HP_n;
	int BS_n;


	Setka(int N1, int N2, int N3, int N4, int M1, int M2, int M3, int M4);
	Setka();
	void Inizialization(void);

	void Write_file_for_FCMHD(void);
	void Culc_ENA(Setka* SS3);
	void Read_file_for_FCMHD(void);
	void friend Pereinterpol(Setka* S1, Setka* S2);


	// Сохранение сетки
	void Save_surface(string name); // Сохранение положений поверхностей разрыва для дальнейшего считывания
	// Для TS сохраняет в полярных координатах, для других поверхностей сохранение и в полярных и в декартовых

	//

	// Считывание положений основных поверхностей с файла и перестройка сетки под эти поверхности
	void Download_surface(string name);
	double find_TS(const double& alp);
	double find_HP_angle(const double& alp);
	double find_HP_x(const double& x);
	double find_BS_x(const double& x);
	double find_BS_angle(const double& alp);
	void Move_Setka_Calculate_stat();  // Здесь нет сгущения точек нулю по умолчанию, но можно организовать


	void normir(int ii);

	void Save_Setka_ALL_ALPHA(string name); // Большая и сложная функция сохранения полной сетки
	void Save_Source_MK(string name);
	void Download_Source_MK(string name);
	void Download_Setka_ALL_ALPHA(string name);       // Старая версия для фалов до 122 включительно
	void Download_Setka_ALL_ALPHA_2_0(string name);   // Новая функция с какого-то момента.... 
	void Copy(Setka* S); // Копирует сетку S (положения разрывов), при этом сетки разного разрешения
	// Для того, чтобы копировать сетки, они должны находиться в одном масштабе (левая и правая граница должна быть такая-же)
	// потому что это параметр не сетки, а параметр самой программы. #define


	void Print_point();      // Печатает точки в сетке (не ячейки, а узлы)
	void Print_Gran();
	void Print_Gran(string name);
	void Print_Gran_type();
	void Print_cell(void);   // Печатает ячейки (но много лишних узлов выписывает, т.к. просто каждую ячейку рисует и линии накладываются
	void Print_cell_type(void);
	void Print_cell2(string name = "Setka_all_cell.txt");  // В отличие от предыдущей не выводит лишних точек, нужна для дополнительной проверки правильности геометрии
	void Print_connect(void); // Посмотреть, как связаны ячейки
	void Print_point_connect(void); // Посмотреть, как связаны узлы с ячейками
	void Print_Tecplot(void);
	void Print_Tecplot_MK(string name0 = "none");
	void Print_Sourse(void);
	void Print_Tecplot_MK_2d(string name0);

	void Print_for_Igor(void);



	void Proverka(void);   // Проверка правильности построения сетки (геометрии и связей). Сюда можно добавлять все тесты сетки в будующем

	// Движение сетки
	void Move_Setka_Calculate(const double& dt);
	void Move_Setka_Calculate_2(const double& dt);
	void Move_Setka_Calculate_3(const double& dt);
	// здесь внутренняя ударная волны выделяется только до Pi/2, а внешняя до -500 а.е. примерно
	void Move_surface(int ii, const double& dt);  // Вычисление скоростей поверхностей   ii - какие параметры активные par[ii]
	void Move_surface_hand(void);  // Ручное движение сетки
	void Smooth_kvadr3(const double& x1, const double& y1, const double& z1, //
		const double& x2, const double& y2, const double& z2, const double& x3, const double& y3,//
		const double& z3, const double& x4, const double& y4, const double& z4, const double& x5,//
		const double& y5, const double& z5, double& xx, double& yy, double& zz);

	void TVD_prepare(void);
	void M_K_prepare(void);
	void Print_TVD(void);

	// Газовая динамика

	void Save_G_D(void);
	void Download_G_D(void);
	void Save_G_D_5_komponent(void);
	void Download_G_D_5_komponent(void);
	void Go_stationary_inner_infty(int step);
	void Go_5_komponent_infty(int step, bool dvig);
	void Go_stationary_5_komponent_inner(int step);
	void Go_stationary_5_komponent_inner_2(int step);
	void Go_stationary_5_komponent_inner_MK(int step);  // Источники берутся из Монте-Карло
	void Go_stationary_5_komponent_inner_MK2(int step);
	void Go_stationary(int step);
	void Go_stationary_TVD(int step);
	void Go_stationary_5_komponent(int step);
	void Go(int step); // Выделение разрывов  Газовая динамика
	void Go_2(int step);
	void Go_3(int step);
	void Go_5_komponent(int step);
	void Go_5_komponent_2(int step);  // Другой способ записи законов сохранения и, соответственно многое поменялось
	void Go_5_komponent_MK(int step, bool dvig = true);  // Версия с источниками из Монте-Карло
	void Go_5_komponent__MK2(int step);
	double HLLC_2d_Korolkov_b_s(const double& ro_L, const double& Q_L, const double& p_L, const double& v1_L, const double& v2_L,//
		const double& ro_R, const double& Q_R, const double& p_R, const double& v1_R, const double& v2_R, const double& W, //
		double* P, double& PQ, const double& n1, const double& n2, const double& rad, int metod, double& Vl, double& Vc, double& Vp, bool nul_potok = false);
	double HLLC_2d_Korolkov_b_s_2(const double& ro_L, const double& Q_L, const double& p_L, const double& v1_L, const double& v2_L,//
		const double& ro_R, const double& Q_R, const double& p_R, const double& v1_R, const double& v2_R, const double& W, //
		double* P, double& PQ, const double& n1, const double& n2, const double& rad, int metod, double& Vl, double& Vc, double& Vp, bool nul_potok = false);
		// Функция возвращает не потоки, а сами большие величины в ячейке. 
	void Init_conditions(void);


	// Монте-Карло

	void Init_Velosity(Sensor* sens, const double& A2, vector <double>& mu, vector <double>& Wt, vector <double>& Wp, vector <double>& Wr, const double& the);
	void Velosity_initial2(Sensor* s, double& Vx, double& Vy, double& Vz); //   Для вылета сзади

	template<typename Random_type, typename Distribution_type, typename accuracy_type>
	void Velosity_initial2(Random_type& gen, Distribution_type& dis, accuracy_type& Vx, accuracy_type& Vy, accuracy_type& Vz);

	void Velosity_initial(Sensor* s, double& Vx, double& Vy, double& Vz);  //  Для вылета с плоскости спереди (4 тип)

	template<typename Random_type, typename Distribution_type, typename accuracy_type>
	void Velosity_initial(Random_type& gen, Distribution_type& dis, accuracy_type& Vx, accuracy_type& Vy, accuracy_type& Vz);

	void Init_Pozision(Sensor* sens, const double& A1, double& phi, double& the);    // Моделирование начальных положений на полусфере
	double F_mk(const double& gamma, const double& Yr);
	void MK_start(void);
	void MK_start_new(void);
	void MK_start_2_0(void);

	// БЛОК ДЛЯ РАБОТЫ НА MPI 
	void MPI_MK_start(int argc, char** argv);  // Версия Монте-Карло на MPI


	Cell* Belong_point(int b, const double& x, const double& y);
	Cell* Find_cell(int& b, const double& x, const double& y);

	void Fly_exchenge(Sensor* sens, double x_0, double y_0, double z_0, double Vx, double Vy, double Vz, Cell* now, //
		double mu, const double& mu_0, bool ExCh);
	void Fly_exchenge_Split(Sensor* sens, double x_0, double y_0, double z_0, double Vx, double Vy, double Vz, Cell* now, double mu,//
		const double& mu_0, bool ExCh, int zone);
	void Change_Velosity(Sensor* s, const double& Ur, const double& Uthe, const double& Uphi, //
		const double& Vr, const double& Vthe, const double& Vphi, double& X, double& Y, double& Z, const double& cp);
	
	void Change_Velosity_Split(Sensor* s, const double& Ur, const double& Uthe, const double& Uphi, //
		const double& Vr, const double& Vthe, const double& Vphi, vector <double>& Wr, vector <double>& Wthe,//
		vector <double>& Wphi, vector <double>& mu_, const double& cp, const double& r, int I,//
		const double& x_ex = 0.0, const double& y_ex = 0.0, const double& z_ex = 0.0);
	
	double w_c_v_s(const double& r, const double& v, int i);
	
	void Fly_exchenge_Imit(MKmethod& MK, Sensor* sens, double x_0, double y_0, double z_0, double Vx, double Vy, double Vz, Cell* now, double mu,//
		 double KSI, double I_do, int area, const double& mu_start, int to_I , int iii); // Имитационный метод
	
	void Fly_exchenge_Imit_Korol(MKmethod& MK, Sensor* sens, double x_0, double y_0, double z_0, double Vx, double Vy, double Vz, Cell* now, double mu, //
		int area, bool ExCh, const double& mu_start, int to_I = 0, int to_J = 0, bool georaschep = true, int zon_stat = -1); // Смотри описание функции в коде функции
	
	void Fly_exchenge_Imit_Korol(MKmethod& MK, int& s1, int& s2, int& s3, double x_0, double y_0, double z_0, double Vx, double Vy, double Vz, Cell* now, double mu, //
		int area, bool ExCh, const double& mu_start, int to_I = 0, int to_J = 0, bool georaschep = true, int zon_stat = -1);

	void Fly_exchenge_Imit_Korol_2(MKmethod& MK, Sensor* sens, double x_0, double y_0, double z_0, double Vx, double Vy, double Vz, Cell* now, double mu, double KSI, //
		double I_do, int area, const double& mu_start);

	void Fly_exchenge_Imit_Korol(MKmethod& MK, Sensor* sens, bool** AZ, double x_0, double y_0, double z_0, double Vx, double Vy, double Vz, Cell* now, double mu, //
		int area, bool ExCh, const double& mu_start, int to_I = 0, int to_J = 0, bool georaschep = true, int zon_stat = -1);

	void Fly_exchenge_Imit_Korol_auto_weight(MKmethod& MK, int& s1, int& s2, int& s3, double x_0, double y_0, double z_0, double Vx, double Vy, double Vz, Cell* now, double mu, //
		int area, bool ExCh, const double& mu_start, int to_I = 0, int to_J = 0, bool georaschep = true, int zon_stat = -1);

	void Fly_exchenge_Imit_Korol_PUI(MKmethod& MK, Sensor* sens, double x_0, double y_0, double z_0, double Vx, double Vy, double Vz, Cell* now, double mu, //
		int area, bool ExCh, const double& mu_start, int to_I, int to_J, bool georaschep, int zon_stat = -1, bool pui__ = false);
	int geo_zones(const double& r, const double& k = 1.0);
	int alpha_zones(const double& x, const double& y);
	double distination(const double& x0, const double& y0, const double& z0,
		const double& Vx, const double& Vy, const double& Vz, const int& sort, int& to_iii, int& to_i, int& to_j);

	void culc_K_Istok(void);
	void proverka_Istok(int ni);

	void culc_PUI(void);
	void GD_prepare(void);
	double get_w_init(const int& k);

	// Поглощение

	void func_pogloshenie(void);

	// Перезарядка с расщеплением на траектории
	double Velosity_1(const double& u, const double& cp);
	double Velosity_2(const double& u, const double& cp);
	double Velosity_3(const double& u, const double& cp);


	// Тест скорости выполнения программы

	void test_velosity(const double& a, const double& b);
	void test_main(void);

private:

};

