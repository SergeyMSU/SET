#pragma once
#include "Help.h"
#include <string>
#include <vector>
#include <mutex>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>

class Dist_func
{
public:
	// Класс для сохренения трёх-мерной функции распределения 
	unsigned short n1; // Должно делится на 2!
	unsigned short n2;
	unsigned short n3;
	double c1;          // Центр распределения
	double d1;          // Ширина распределения
	double a1;
	double c2;
	double d2;
	double a2;
	double c3;
	double d3;
	double a3;

	double xxx = 0.0;        // В какой точке мерится функция распределения
	double yyy = 0.0;

	string name;

	double*** V;

	Dist_func(const unsigned short& a1, const unsigned short& a2, const unsigned short& a3, //
		const double& c1, const double& d1, //
		const double& c2, const double& d2, //
		const double& c3, const double& d3);
	bool call_name(string name);         // Задать имя (лучше номер ячейки)
	bool Add_point(const double& V1xyz, const double& V2xyz, const double& V3xyz, const double& y, const double& z, const double& mu);
	bool normir(const double& ccc);

	void v_cyl_to_v_xyz(double v_rho, double v_phi, double x, double y, double& vx, double& vy); //z направлен в upwind!!!!!!! x,y - остальные две оси
	void v_cyl(double vx, double vy, double x, double y, double& w_rho, double& w_phi);

	bool print_1d(int koord); // По какой координате срез

	bool print_3d(void); // По какой координате срез

private:
};

