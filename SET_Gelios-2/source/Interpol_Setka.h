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

class Setka;
class Point;
class Gran;
class Cell;

using namespace std;

class Interpol_Setka
{
public:
	vector <Point*> All_Points;   // Все точки (центры ячеек исходной сетки + дополнительные)
	vector <Gran*> All_Gran;           // Грани основные
	vector <Gran*> All_Gran_copy;      // Грани фантомные (нормаль в другую сторону)
	vector <Cell*> All_Cells;          // Все ячейки

	Interpol_Setka(Setka* SS);     // Построение интерполяционной сетки из основной SS

	void Print_Cell(string name);       // Печать сетки
};

