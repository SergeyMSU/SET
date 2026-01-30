#include "Setka.h"
#include <algorithm>
#include <omp.h>
#include <mutex>
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>

#if USEMPI 
#include "mpi.h"
#endif

using namespace std;

template void Setka::Velosity_initial2(std::mt19937& gen, std::uniform_real_distribution<double>& dis, double& Vx, double& Vy, double& Vz);
template void Setka::Velosity_initial2(std::mt19937& gen, std::uniform_real_distribution<float>& dis, float& Vx, float& Vy, float& Vz);
template void Setka::Velosity_initial(std::mt19937& gen, std::uniform_real_distribution<double>& dis, double& Vx, double& Vy, double& Vz);
template void Setka::Velosity_initial(std::mt19937& gen, std::uniform_real_distribution<float>& dis, float& Vx, float& Vy, float& Vz);
	
Setka::Setka()
{
	this->Inizialization();
}

Setka::Setka(int N1, int N2, int N3, int N4, int M1, int M2, int M3, int M4)
// Начальная инициализация сетки,  эту функцию лучше не трогать, она писалась поэтапно и результаты постоянно проверялись
// Изменение чего-либо может привести к ошибке
{
	// Число А - лучей, Число B - лучей, 
	this->Inizialization();
	this->N1 = N1;
	this->N2 = N2;
	this->N3 = N3;
	this->N4 = N4;
	this->M1 = M1;
	this->M2 = M2;
	this->M3 = M3;
	this->M4 = M4;
	N1++;
	N3++;
	M1++;
	double s = 0.0;
	int kk, kkk, kkkk;
	int n = M1 + M2 + M3 + M4 - 1; // Число граней по лучу. 
	kk = n * (N2 + N1 - 2) + (N3 - 1) * (N4 + M1 - 1);
	kkk = kk + (N4 - 1) * (M2 + M3 + M4);

	auto Cell_Centr = new Cell(); // Центральная особая ячейка (создаётся здесь отдельно, чтобы положить в общий массив только в конце связки всего)

	for (int i = 0; i < N1; i++)
	{
		//double x = 1.0 * i / (N1 - 1);
		//s = (0.35 * x * (1.0 - x) + kv(x)) * (pi_ / 2.0);  //i * (pi_ / 2.0) / (N1 - 1.0);
		s = i * (pi_ / 2.0) / (N1 - 1.0);
		auto R = new Rail(s);
		R->type = A;
		this->A_Rails.push_back(R);
		R->M1 = M1;
		R->M2 = M2;
		R->M3 = M3;
		R->M4 = M4;

		R->Init_start(R);
		for (auto& j : R->All_point)
		{
			this->All_Points.push_back(j);
		}
	}

	for (int i = 0; i < N2; i++)
	{
		//s = (pi_ / 2.0) + (i + 1) * (pi_ / 4.0) / (N2);
		s = (pi_ / 2.0) + (i + 1) * (pi_ / 3.0) / (N2);
		auto R = new Rail(s);
		R->type = B;
		this->B_Rails.push_back(R);
		R->M1 = M1;
		R->M2 = M2;
		R->M3 = M3;
		R->M4 = M4;
		R->Init_start(R);
		for (auto& j : R->All_point)
		{
			this->All_Points.push_back(j);
		}
	}

	for (int i = 0; i < N3; i++)
	{
		//s = (3.0 * pi_ / 4.0) + (i) * (pi_ / 4.0) / (N3 - 1);
		s = (5.0 * pi_ / 6.0) + (i) * (pi_ / 6.0) / (N3 - 1);
		auto R = new Rail(s);
		R->type = C;
		this->C_Rails.push_back(R);
		R->M1 = M1;
		R->M2 = N4;
		R->M3 = 0.0;
		R->M4 = 0.0;
		if (i > 0)
		{
			R->Init_start(R);
		}
	}

	for (int i = 0; i < N4; i++)
	{
		s = -R1_ - (i + 1) * (-R1_ - Left_) / (N4);
		auto R = new Rail(s);
		R->type = D;
		this->D_Rails.push_back(R);
		R->M1 = N3;
		R->M2 = M2;
		R->M3 = M3;
		R->M4 = M4;
	}

	// Отдельно нужно обработать С-type линию от тройной точки
	auto C = this->C_Rails[0];
	for (int i = 0; i < M1; i++)
	{
		C->All_point.push_back(this->B_Rails[N2 - 1]->All_point[i]);
		if (i == M1 - 1)
		{
			C->Key_point.push_back(this->B_Rails[N2 - 1]->All_point[i]);
		}
	}
	double l;
	double pr = (5.0 * pi_ / 6.0);
	for (int i = 0; i < N4; i++)
	{
		//double x = pow(Left_ / (R2_ * cos(3.0 * pi_ / 4.0)), 1.0 / (N4));
		//l = R2_ * cos(3.0 * pi_ / 4.0) * pow(x, i + 1); ; // R2_* cos(3.0 * pi_ / 4.0) - (i + 1) * (R2_ * cos(3.0 * pi_ / 4.0) - Left_) / N4;
		//auto K = new Point(l, R2_ * sin(3.0 * pi_ / 4.0));
		double x = pow(Left_ / (R2_ * cos(pr)), 1.0 / (N4));
		l = R2_ * cos(pr) * pow(x, i + 1); ; // R2_* cos(3.0 * pi_ / 4.0) - (i + 1) * (R2_ * cos(3.0 * pi_ / 4.0) - Left_) / N4;
		auto K = new Point(l, R2_ * sin(pr));
		C->All_point.push_back(K);
		this->All_Points.push_back(K);
		K->type = P_U_2;
		if (i == N4 - 1)
		{
			C->Key_point.push_back(K);
			K->type = P_Outer_Boandary;
		}
	}
	// Обработка завершена

	for (int i = 1; i < N3; i++)
	{
		for (auto& j : this->C_Rails[i]->All_point)
		{
			this->All_Points.push_back(j);
		}
	}

	// Отдельно обрабатываем четвёртый тип линий
	for (int i = 0; i < N4; i++)
	{
		auto P = C->All_point[M1 + i];
		this->D_Rails[i]->All_point.push_back(P);
		this->D_Rails[i]->Key_point.push_back(P);
		for (int j = 0; j < M2; j++)
		{
			l = P->y + (j + 1) * (R3_ - P->y) / (M2);
			auto K = new Point(P->x, l);
			this->D_Rails[i]->All_point.push_back(K);
			this->All_Points.push_back(K);
			K->type = P_U_3;
			if (j == M2 - 1)
			{
				this->D_Rails[i]->Key_point.push_back(K);
				K->type = P_Contact;
			}
		}
		for (int j = 0; j < M3; j++)
		{
			l = R3_ + (j + 1) * (R4_ - R3_) / M3;
			auto K = new Point(P->x, l);
			this->D_Rails[i]->All_point.push_back(K);
			this->All_Points.push_back(K);
			K->type = P_U_4;
			if (j == M3 - 1)
			{
				this->D_Rails[i]->Key_point.push_back(K);
				K->type = P_Outer_shock;
			}
		}
		double x = pow(R5_ / R4_, 1.0 / (M4));
		for (int j = 0; j < M4; j++)
		{
			l = R4_ * pow(x, j + 1);  //R4_ + (j + 1) * (R5_ - R4_) / M4;
			auto K = new Point(P->x, l);
			this->D_Rails[i]->All_point.push_back(K);
			this->All_Points.push_back(K);
			K->type = P_U_5;
			if (j == M4 - 1)
			{
				this->D_Rails[i]->Key_point.push_back(K);
				K->type = P_Outer_Boandary;
			}
		}
	}

	// Теперь начинается создание самих ячеек
	for (int i = 0; i < this->A_Rails.size() - 1; i++)
	{
		for (int j = 0; j < this->A_Rails[0]->All_point.size() - 1; j++)
		{
			auto C = new Cell(this->A_Rails[i]->All_point[j], this->A_Rails[i]->All_point[j + 1], //
				this->A_Rails[i + 1]->All_point[j + 1], this->A_Rails[i + 1]->All_point[j]);
			this->A_Rails[i]->All_point[j]->my_cell.push_back(C);
			this->A_Rails[i]->All_point[j + 1]->my_cell.push_back(C);
			this->A_Rails[i + 1]->All_point[j + 1]->my_cell.push_back(C);
			this->A_Rails[i + 1]->All_point[j]->my_cell.push_back(C);
			this->All_Cells.push_back(C);
		}
	}

	for (int i = 0; i < this->B_Rails.size() - 1; i++)
	{
		for (int j = 0; j < this->B_Rails[0]->All_point.size() - 1; j++)
		{
			auto C = new Cell(this->B_Rails[i]->All_point[j], this->B_Rails[i]->All_point[j + 1], //
				this->B_Rails[i + 1]->All_point[j + 1], this->B_Rails[i + 1]->All_point[j]);
			this->B_Rails[i]->All_point[j]->my_cell.push_back(C);
			this->B_Rails[i]->All_point[j + 1]->my_cell.push_back(C);
			this->B_Rails[i + 1]->All_point[j + 1]->my_cell.push_back(C);
			this->B_Rails[i + 1]->All_point[j]->my_cell.push_back(C);
			this->All_Cells.push_back(C);
		}
	}

	for (int i = 0; i < this->C_Rails.size() - 1; i++)
	{
		for (int j = 0; j < this->C_Rails[0]->All_point.size() - 1; j++)
		{
			auto C = new Cell(this->C_Rails[i]->All_point[j], this->C_Rails[i]->All_point[j + 1], //
				this->C_Rails[i + 1]->All_point[j + 1], this->C_Rails[i + 1]->All_point[j]);
			this->C_Rails[i]->All_point[j]->my_cell.push_back(C);
			this->C_Rails[i]->All_point[j + 1]->my_cell.push_back(C);
			this->C_Rails[i + 1]->All_point[j + 1]->my_cell.push_back(C);
			this->C_Rails[i + 1]->All_point[j]->my_cell.push_back(C);
			this->All_Cells.push_back(C);
		}
	}

	for (int i = 0; i < this->D_Rails.size() - 1; i++)
	{
		for (int j = 0; j < this->D_Rails[0]->All_point.size() - 1; j++)
		{
			auto C = new Cell(this->D_Rails[i]->All_point[j], this->D_Rails[i]->All_point[j + 1], //
				this->D_Rails[i + 1]->All_point[j + 1], this->D_Rails[i + 1]->All_point[j]);
			this->D_Rails[i]->All_point[j]->my_cell.push_back(C);
			this->D_Rails[i]->All_point[j + 1]->my_cell.push_back(C);
			this->D_Rails[i + 1]->All_point[j + 1]->my_cell.push_back(C);
			this->D_Rails[i + 1]->All_point[j]->my_cell.push_back(C);
			this->All_Cells.push_back(C);
		}
	}


	// Связка по-середине
	auto P1 = this->A_Rails[this->A_Rails.size() - 1];
	auto P2 = this->B_Rails[0];
	for (int j = 0; j < P1->All_point.size() - 1; j++)
	{
		auto C = new Cell(P1->All_point[j], P1->All_point[j + 1], //
			P2->All_point[j + 1], P2->All_point[j]);
		P1->All_point[j]->my_cell.push_back(C);
		P1->All_point[j + 1]->my_cell.push_back(C);
		P2->All_point[j + 1]->my_cell.push_back(C);
		P2->All_point[j]->my_cell.push_back(C);
		this->All_Cells.push_back(C);
	}

	// Связка от тройной точки
	P1 = this->B_Rails[this->B_Rails.size() - 1];
	P2 = this->D_Rails[0];
	for (int j = 0; j < P2->All_point.size() - 1; j++)
	{
		auto C = new Cell(P1->All_point[j + M1 - 1], P1->All_point[j + M1], //
			P2->All_point[j + 1], P2->All_point[j]);
		P1->All_point[j + M1 - 1]->my_cell.push_back(C);
		P1->All_point[j + M1]->my_cell.push_back(C);
		P2->All_point[j + 1]->my_cell.push_back(C);
		P2->All_point[j]->my_cell.push_back(C);
		this->All_Cells.push_back(C);
	}

	// На этом этапе все узлы связаны, ячейки созданы. Теперь нужно связать ячейки. Создать грани. ---------------------------------------------

	// Пронумеруем узлы и ячейки.
	int num = 0;
	for (auto& i : this->All_Cells)
	{
		i->number = num;
		num++;
	}

	num = 0;  
	// Нумеруем точки
	for (auto& i : this->All_Points)
	{
		i->number = num;
		num++;
	}

	n = M1 + M2 + M3 + M4 - 1; // Число граней по лучу. 

	if (true)
	{

		for (int i = 0; i < n * (N1 - 1); i++) // Связываем ячейки по лучу (продольная связка)
		{
			if ((i + 1) % (n) == 0)
			{
				auto K = this->All_Cells[i];
				auto G = new Gran(K->contour[1], K->contour[2], Input);
				K->Grans.push_back(G);
				G->main_gran = true;
				G->Master = K;
				this->All_Gran.push_back(G);
				continue;
			}
			auto K = this->All_Cells[i];
			auto G = new Gran(K->contour[1], K->contour[2], Usualy);
			auto G2 = new Gran(K->contour[2], K->contour[1], Usualy);
			G->main_gran = true;
			G2->main_gran = false;
			this->All_Gran.push_back(G);
			this->All_Gran_copy.push_back(G2);
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			G->Master = K;
			G->Sosed = this->All_Cells[i + 1];
			G2->Master = this->All_Cells[i + 1];
			G2->Sosed = K;
			K->Grans.push_back(G);
			this->All_Cells[i + 1]->Grans.push_back(G2);
		}

		for (int i = 0; i < n * (N1 - 2); i++)  // Связываем поперёк
		{
			auto K = this->All_Cells[i];
			auto G = new Gran(K->contour[2], K->contour[3], Usualy);
			auto G2 = new Gran(K->contour[3], K->contour[2], Usualy);
			G->main_gran = true;
			G2->main_gran = false;
			this->All_Gran.push_back(G);
			this->All_Gran_copy.push_back(G2);
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			G->Master = K;
			G->Sosed = this->All_Cells[i + n];
			G2->Master = this->All_Cells[i + n];
			G2->Sosed = K;
			K->Grans.push_back(G);
			this->All_Cells[i + n]->Grans.push_back(G2);
		}


		for (int i = n * (N1 - 1); i < n * (N2 + N1 - 2); i++) // Связываем ячейки по лучу (продольная связка)
		{
			if ((i + 1) % (n) == 0)
			{
				auto K = this->All_Cells[i];
				auto G = new Gran(K->contour[1], K->contour[2], Upper_wall);
				G->main_gran = true;
				this->All_Gran.push_back(G);
				K->Grans.push_back(G);
				G->Master = K;
				continue;
			}
			auto K = this->All_Cells[i];
			auto G = new Gran(K->contour[1], K->contour[2], Usualy);
			auto G2 = new Gran(K->contour[2], K->contour[1], Usualy);
			G->main_gran = true;
			G2->main_gran = false;
			this->All_Gran.push_back(G);
			this->All_Gran_copy.push_back(G2);
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			G->Master = K;
			G->Sosed = this->All_Cells[i + 1];
			G2->Master = this->All_Cells[i + 1];
			G2->Sosed = K;
			K->Grans.push_back(G);
			this->All_Cells[i + 1]->Grans.push_back(G2);
		}

		for (int i = n * (N1 - 1); i < n * (N2 + N1 - 3); i++)  // Связываем поперёк
		{
			auto K = this->All_Cells[i];
			auto G = new Gran(K->contour[2], K->contour[3], Usualy);
			auto G2 = new Gran(K->contour[3], K->contour[2], Usualy);
			G->main_gran = true;
			G2->main_gran = false;
			this->All_Gran.push_back(G);
			this->All_Gran_copy.push_back(G2);
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			G->Master = K;
			G->Sosed = this->All_Cells[i + n];
			G2->Master = this->All_Cells[i + n];
			G2->Sosed = K;
			K->Grans.push_back(G);
			this->All_Cells[i + n]->Grans.push_back(G2);
		}


		for (int i = n * (N2 + N1 - 2); i < n * (N2 + N1 - 2) + (N3 - 1) * (N4 + M1 - 1); i++) // Связываем ячейки по лучу (продольная связка)
		{
			if ((i + 1 - n * (N2 + N1 - 2)) % (N4 + M1 - 1) == 0)
			{
				auto K = this->All_Cells[i];
				auto G = new Gran(K->contour[1], K->contour[2], Extern);
				G->main_gran = true;
				this->All_Gran.push_back(G);
				K->Grans.push_back(G);
				G->Master = K;
				continue;
			}
			auto K = this->All_Cells[i];
			auto G = new Gran(K->contour[1], K->contour[2], Usualy);
			auto G2 = new Gran(K->contour[2], K->contour[1], Usualy);
			G->main_gran = true;
			G2->main_gran = false;
			this->All_Gran.push_back(G);
			this->All_Gran_copy.push_back(G2);
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			G->Master = K;
			G->Sosed = this->All_Cells[i + 1];
			G2->Master = this->All_Cells[i + 1];
			G2->Sosed = K;
			K->Grans.push_back(G);
			this->All_Cells[i + 1]->Grans.push_back(G2);
		}

		for (int i = n * (N2 + N1 - 2); i < n * (N2 + N1 - 2) + (N3 - 2) * (N4 + M1 - 1); i++)  // Связываем поперёк
		{
			auto K = this->All_Cells[i];
			auto G = new Gran(K->contour[2], K->contour[3], Usualy);
			auto G2 = new Gran(K->contour[3], K->contour[2], Usualy);
			G->main_gran = true;
			G2->main_gran = false;
			this->All_Gran.push_back(G);
			this->All_Gran_copy.push_back(G2);
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			G->Master = K;
			G->Sosed = this->All_Cells[i + (N4 + M1 - 1)];
			G2->Master = this->All_Cells[i + (N4 + M1 - 1)];
			G2->Sosed = K;
			K->Grans.push_back(G);
			this->All_Cells[i + (N4 + M1 - 1)]->Grans.push_back(G2);
		}


		kk = n * (N2 + N1 - 2) + (N3 - 1) * (N4 + M1 - 1);

		for (int i = kk; i < kk + (N4 - 1) * (M2 + M3 + M4); i++) // Связываем ячейки по лучу (продольная связка)
		{
			if ((i + 1 - kk) % (M2 + M3 + M4) == 0)
			{
				auto K = this->All_Cells[i];
				auto G = new Gran(K->contour[1], K->contour[2], Upper_wall);
				G->main_gran = true;
				K->Grans.push_back(G);
				G->Master = K;
				this->All_Gran.push_back(G);
				continue;
			}
			auto K = this->All_Cells[i];
			auto G = new Gran(K->contour[1], K->contour[2], Usualy);
			auto G2 = new Gran(K->contour[2], K->contour[1], Usualy);
			G->main_gran = true;
			G2->main_gran = false;
			this->All_Gran.push_back(G);
			this->All_Gran_copy.push_back(G2);
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			G->Master = K;
			G->Sosed = this->All_Cells[i + 1];
			G2->Master = this->All_Cells[i + 1];
			G2->Sosed = K;
			K->Grans.push_back(G);
			this->All_Cells[i + 1]->Grans.push_back(G2);
		}

		for (int i = kk; i < kk + (N4 - 2) * (M2 + M3 + M4); i++)  // Связываем поперёк
		{
			auto K = this->All_Cells[i];
			auto G = new Gran(K->contour[2], K->contour[3], Usualy);
			auto G2 = new Gran(K->contour[3], K->contour[2], Usualy);
			G->main_gran = true;
			G2->main_gran = false;
			this->All_Gran.push_back(G);
			this->All_Gran_copy.push_back(G2);
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			G->Master = K;
			G->Sosed = this->All_Cells[i + (M2 + M3 + M4)];
			G2->Master = this->All_Cells[i + (M2 + M3 + M4)];
			G2->Sosed = K;
			K->Grans.push_back(G);
			this->All_Cells[i + (M2 + M3 + M4)]->Grans.push_back(G2);
		}


		// Связываем шов поцентру
		kkk = kk + (N4 - 1) * (M2 + M3 + M4);

		for (int i = kkk; i < kkk + M1 + M2 + M3 + M4 - 2; i++)
		{
			auto K = this->All_Cells[i];
			auto G = new Gran(K->contour[1], K->contour[2], Usualy);
			auto G2 = new Gran(K->contour[2], K->contour[1], Usualy);
			G->main_gran = true;
			G2->main_gran = false;
			this->All_Gran.push_back(G);
			this->All_Gran_copy.push_back(G2);
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			G->Master = K;
			G->Sosed = this->All_Cells[i + 1];
			G2->Master = this->All_Cells[i + 1];
			G2->Sosed = K;
			K->Grans.push_back(G);
			this->All_Cells[i + 1]->Grans.push_back(G2);
		}

		for (int i = kkk; i < kkk + M1 + M2 + M3 + M4 - 1; i++)
		{
			auto K = this->All_Cells[i];
			auto F = this->All_Cells[n * (N1 - 2) + i - kkk];
			auto G = new Gran(K->contour[0], K->contour[1], Usualy);
			auto G2 = new Gran(K->contour[1], K->contour[0], Usualy);
			G->main_gran = true;
			G2->main_gran = false;
			this->All_Gran.push_back(G);
			this->All_Gran_copy.push_back(G2);
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			G->Master = K;
			G->Sosed = F;
			G2->Master = F;
			G2->Sosed = K;
			K->Grans.push_back(G);
			F->Grans.push_back(G2);
		}

		for (int i = kkk; i < kkk + M1 + M2 + M3 + M4 - 1; i++)
		{
			auto K = this->All_Cells[i];
			auto F = this->All_Cells[n * (N1 - 1) + i - kkk];
			auto G = new Gran(K->contour[2], K->contour[3], Usualy);
			auto G2 = new Gran(K->contour[3], K->contour[2], Usualy);
			G->main_gran = true;
			G2->main_gran = false;
			this->All_Gran.push_back(G);
			this->All_Gran_copy.push_back(G2);
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			G->Master = K;
			G->Sosed = F;
			G2->Master = F;
			G2->Sosed = K;
			K->Grans.push_back(G);
			F->Grans.push_back(G2);
		}


		kkkk = n * (N2 + N1 - 2);
		for (int i = kkkk; i < kkkk + M1 - 1; i++)
		{
			auto K = this->All_Cells[i];
			auto F = this->All_Cells[i - n];
			auto G = new Gran(K->contour[0], K->contour[1], Usualy);
			auto G2 = new Gran(K->contour[1], K->contour[0], Usualy);
			G->main_gran = true;
			G2->main_gran = false;
			this->All_Gran.push_back(G);
			this->All_Gran_copy.push_back(G2);
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			G->Master = K;
			G->Sosed = F;
			G2->Master = F;
			G2->Sosed = K;
			K->Grans.push_back(G);
			F->Grans.push_back(G2);
		}


		for (int i = this->All_Cells.size() - M2 - M3 - M4; i < this->All_Cells.size() - 1; i++)
		{
			auto K = this->All_Cells[i];
			auto G = new Gran(K->contour[1], K->contour[2], Usualy);
			auto G2 = new Gran(K->contour[2], K->contour[1], Usualy);
			G->main_gran = true;
			G2->main_gran = false;
			this->All_Gran.push_back(G);
			this->All_Gran_copy.push_back(G2);
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			G->Master = K;
			G->Sosed = this->All_Cells[i + 1];
			G2->Master = this->All_Cells[i + 1];
			G2->Sosed = K;
			K->Grans.push_back(G);
			this->All_Cells[i + 1]->Grans.push_back(G2);
		}

		for (int i = this->All_Cells.size() - M2 - M3 - M4; i < this->All_Cells.size(); i++)
		{
			auto K = this->All_Cells[i];
			auto F = this->All_Cells[kk + i - this->All_Cells.size() + M2 + M3 + M4];
			auto G = new Gran(K->contour[2], K->contour[3], Usualy);
			auto G2 = new Gran(K->contour[3], K->contour[2], Usualy);
			G->main_gran = true;
			G2->main_gran = false;
			this->All_Gran.push_back(G);
			this->All_Gran_copy.push_back(G2);
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			G->Master = K;
			G->Sosed = F;
			G2->Master = F;
			G2->Sosed = K;
			K->Grans.push_back(G);
			F->Grans.push_back(G2);
		}

		for (int i = this->All_Cells.size() - M2 - M3 - M4; i < this->All_Cells.size(); i++)
		{
			auto K = this->All_Cells[i];
			auto F = this->All_Cells[kkkk + M1 - n + i - (this->All_Cells.size() - M2 - M3 - M4) - 1];
			auto G = new Gran(K->contour[0], K->contour[1], Usualy);
			auto G2 = new Gran(K->contour[1], K->contour[0], Usualy);
			G->main_gran = true;
			G2->main_gran = false;
			this->All_Gran.push_back(G);
			this->All_Gran_copy.push_back(G2);
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			G->Master = K;
			G->Sosed = F;
			G2->Master = F;
			G2->Sosed = K;
			K->Grans.push_back(G);
			F->Grans.push_back(G2);
		}

		if (true)  // Ручное связывание одного из узлов
		{
			auto K = this->All_Cells[this->All_Cells.size() - M2 - M3 - M4];
			auto F = this->All_Cells[kkkk + M1 - 1];
			auto G = new Gran(K->contour[3], K->contour[0], Usualy);
			auto G2 = new Gran(K->contour[0], K->contour[3], Usualy);
			G->main_gran = true;
			G2->main_gran = false;
			this->All_Gran.push_back(G);
			this->All_Gran_copy.push_back(G2);
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			G->Master = K;
			G->Sosed = F;
			G2->Master = F;
			G2->Sosed = K;
			K->Grans.push_back(G);
			F->Grans.push_back(G2);
		}

		for (int i = kkkk + M1; i < kkkk + M1 + N4 - 1; i++)  // Последний шов от тройной точки вбок
		{
			auto K = this->All_Cells[i];
			auto F = this->All_Cells[kk + (i - kkkk - M1) * (M2 + M3 + M4)];
			auto G = new Gran(K->contour[0], K->contour[1], Usualy);
			auto G2 = new Gran(K->contour[1], K->contour[0], Usualy);
			G->main_gran = true;
			G2->main_gran = false;
			this->All_Gran.push_back(G);
			this->All_Gran_copy.push_back(G2);
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			G->Master = K;
			G->Sosed = F;
			G2->Master = F;
			G2->Sosed = K;
			K->Grans.push_back(G);
			F->Grans.push_back(G2);
		}

		// Добавляем грани - границы

		// Начнём с оси симметрии

		for (int i = 0; i < M1 + M2 + M3 + M4 - 1; i++)
		{
			auto K = this->All_Cells[i];
			auto G = new Gran(K->contour[0], K->contour[1], Axis);
			auto G2 = new Gran(K->contour[1], K->contour[0], Axis);
			G->main_gran = true;
			G2->main_gran = false;
			K->Grans.push_back(G);
			G->Master = K;
			G2->Master = nullptr;
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			G2->Sosed = K;
			this->All_Gran.push_back(G);
			this->All_Gran_copy.push_back(G2);
		}

		for (int i = n * (N2 + N1 - 2) + (N3 - 2) * (N4 + M1 - 1); i < n * (N2 + N1 - 2) + (N3 - 1) * (N4 + M1 - 1); i++)  // Связываем поперёк
		{
			auto K = this->All_Cells[i];
			auto G = new Gran(K->contour[2], K->contour[3], Axis);
			auto G2 = new Gran(K->contour[3], K->contour[2], Axis);
			G->main_gran = true;
			G2->main_gran = false;
			K->Grans.push_back(G);
			G->Master = K;
			G2->Master = nullptr;
			G2->Sosed = K;
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			this->All_Gran.push_back(G);
			this->All_Gran_copy.push_back(G2);
		}


		for (int i = 0; i < N1 + N2 - 2; i++)
		{
			auto K = this->All_Cells[i * (M1 + M2 + M3 + M4 - 1)];
			auto G = new Gran(K->contour[3], K->contour[0], Inner_sphere);
			auto G2 = new Gran(K->contour[0], K->contour[3], Inner_sphere);  // Usualy);
			G->main_gran = true;
			G2->main_gran = false;
			K->Grans.push_back(G);
			Cell_Centr->Grans.push_back(G2);
			G->Master = K;
			G2->Master = Cell_Centr;
			G->Sosed = Cell_Centr;
			G2->Sosed = K;
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			this->All_Gran.push_back(G);
			this->All_Gran_copy.push_back(G2);
		}


		// Шов по центру с верхней стенкой
		if (true)
		{
			auto K = this->All_Cells[kkk + M1 + M2 + M3 + M4 - 2];
			auto G = new Gran(K->contour[1], K->contour[2], Upper_wall);
			auto G2 = new Gran(K->contour[2], K->contour[1], Upper_wall);
			G->main_gran = true;
			G2->main_gran = false;
			K->Grans.push_back(G);
			G->Master = K;
			G2->Sosed = K;
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			this->All_Gran.push_back(G);
			this->All_Gran_copy.push_back(G2);
		}

		if (true)
		{
			auto K = this->All_Cells[this->All_Cells.size() - 1];
			auto G = new Gran(K->contour[1], K->contour[2], Upper_wall);
			auto G2 = new Gran(K->contour[2], K->contour[1], Upper_wall);
			G->main_gran = true;
			G2->main_gran = false;
			K->Grans.push_back(G);
			G->Master = K;
			G2->Sosed = K;
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			this->All_Gran.push_back(G);
			this->All_Gran_copy.push_back(G2);
		}

	}

	int h = (N1 + N2 - 2) * (M1 + M2 + M3 + M4 - 1);
	for (int i = 0; i < N3 - 1; i++)
	{
		auto K = this->All_Cells[i * (M1 + N4 - 1) + h];
		auto G = new Gran(K->contour[3], K->contour[0], Inner_sphere);
		auto G2 = new Gran(K->contour[0], K->contour[3], Inner_sphere);  //Usualy);
		G->main_gran = true;
		G2->main_gran = false;
		K->Grans.push_back(G);
		Cell_Centr->Grans.push_back(G2);
		G->Master = K;
		G2->Master = Cell_Centr;
		G->Sosed = Cell_Centr;
		G2->Sosed = K;
		G->Gran_copy = G2;
		G2->Gran_copy = G;
		this->All_Gran.push_back(G);
		this->All_Gran_copy.push_back(G2);
	}

	// Шов по центру связываем с внутренней границей
	auto K = this->All_Cells[kkk];
	auto G = new Gran(K->contour[3], K->contour[0], Inner_sphere);
	auto G2 = new Gran(K->contour[0], K->contour[3], Inner_sphere); // Usualy);
	G->main_gran = true;
	G2->main_gran = false;
	K->Grans.push_back(G);
	Cell_Centr->Grans.push_back(G2);
	G->Master = K;
	G2->Master = Cell_Centr;
	G->Sosed = Cell_Centr;
	G2->Sosed = K;
	G->Gran_copy = G2;
	G2->Gran_copy = G;
	this->All_Gran.push_back(G);
	this->All_Gran_copy.push_back(G2);

	for (int i = kk + (N4 - 2) * (M2 + M3 + M4); i < kk + (N4 - 1) * (M2 + M3 + M4); i++) // Связываем ячейки по лучу (продольная связка)
	{
		auto K = this->All_Cells[i];
		auto G = new Gran(K->contour[2], K->contour[3], Extern);
		auto G2 = new Gran(K->contour[3], K->contour[2], Extern);
		G->main_gran = true;
		G2->main_gran = true;
		K->Grans.push_back(G);
		G->Master = K;
		G2->Sosed = K;
		G->Gran_copy = G2;
		G2->Gran_copy = G;
		this->All_Gran.push_back(G);
		this->All_Gran_copy.push_back(G2);
	}

	for (auto& i : this->All_Cells)  // Нужно для выделения контактной поверхности
	{
		if (i->contour[0]->type == P_Inner_Boandary || i->contour[0]->type == P_U_1)
		{
			i->type = C_1;
		}
		else if (i->contour[3]->type == P_U_2 || i->contour[2]->type == P_U_2)
		{
			i->type = C_2;
		}
		
		if (i->contour[0]->type == P_U_3 || i->contour[1]->type == P_U_3 || i->contour[2]->type == P_U_3)
		{
			i->type = C_3;
		}
		else if (i->contour[0]->type == P_U_4 || i->contour[1]->type == P_U_4 || i->contour[2]->type == P_U_4)
		{
			i->type = C_4;
		}
		else if(i->contour[0]->type == P_U_5 || i->contour[1]->type == P_U_5 || i->contour[2]->type == P_U_5)
		{
			i->type = C_5;
		}
	}


	// Хотим добавить особые ячейки.
	vector<Cell*> centr_cell;
	if (true)
	{
		sort(Cell_Centr->Grans.begin(), Cell_Centr->Grans.end(), [&](Gran* i, Gran* j)
			{
				double x1, y1, x2, y2;
				i->Get_Center(x1, y1);
				j->Get_Center(x2, y2);
				return (polar_angle(x1, y1) < polar_angle(x2, y2));
			});

		int kl = 8;
		while (Cell_Centr->Grans.size() % kl != 0)
		{
			kl++;
		}
		int how = Cell_Centr->Grans.size() / kl;

		for (int i = 0; i < kl; i++)
		{
			auto A = new Cell();
			A->type = C_centr;
			for (int j = i * how; j < (i + 1) * how; j++)
			{
				A->Grans.push_back(Cell_Centr->Grans[j]);
				Cell_Centr->Grans[j]->Master = A;
				Cell_Centr->Grans[j]->Gran_copy->Sosed = A;
			}
			centr_cell.push_back(A);
		}

		auto O = new Point(0.0, 0.0);
		O->type = P_Null;
		this->All_Points.push_back(O);

		for (auto& i : centr_cell)
		{
			i->contour.push_back(O);
			O->my_cell.push_back(i);
			for (auto& j : i->Grans)
			{
				i->contour.push_back(j->A);
				j->A->my_cell.push_back(i);
			}
			i->contour.push_back(i->Grans[i->Grans.size() - 1]->B);
			i->Grans[i->Grans.size() - 1]->B->my_cell.push_back(i);
		}

		for (int i = 0; i < centr_cell.size() - 1; i++)
		{
			auto K = centr_cell[i];
			int vm = K->contour.size();
			auto G = new Gran(K->contour[vm - 1], K->contour[0], Usualy);
			auto G2 = new Gran(K->contour[0], K->contour[vm - 1], Usualy);
			G->main_gran = true;
			G2->main_gran = false;
			K->Grans.push_back(G);
			centr_cell[i+1]->Grans.push_back(G2);
			G->Master = K;
			G2->Master = centr_cell[i + 1];
			G->Sosed = centr_cell[i + 1];
			G2->Sosed = K;
			G->Gran_copy = G2;
			G2->Gran_copy = G;
			this->All_Gran.push_back(G);
			this->All_Gran_copy.push_back(G2);
		}

		auto PO = centr_cell[0];
		auto G2 = new Gran(PO->contour[0], PO->contour[1], Axis);
		auto G3 = new Gran(PO->contour[1], PO->contour[0], Axis);
		G2->main_gran = true;
		G3->main_gran = false;
		PO->Grans.push_back(G2);
		G2->Master = PO;
		G3->Sosed = PO;
		G3->Gran_copy = G2;
		G2->Gran_copy = G3;
		this->All_Gran.push_back(G2);
		this->All_Gran_copy.push_back(G3);

		PO = centr_cell[centr_cell.size() - 1];
		G2 = new Gran(PO->contour[PO->contour.size() - 1], PO->contour[0], Axis);
		G3 = new Gran(PO->contour[0], PO->contour[PO->contour.size() - 1], Axis);
		G2->main_gran = true;
		G3->main_gran = false;
		PO->Grans.push_back(G2);
		G2->Master = PO;
		G3->Sosed = PO;
		G3->Gran_copy = G2;
		G2->Gran_copy = G3;
		this->All_Gran.push_back(G2);
		this->All_Gran_copy.push_back(G3);


		/*this->All_Cells.push_back(Cell_Centr);
		Cell_Centr->type = C_centr;
		for (int i = 0; i < this->A_Rails.size(); i++)
		{
			Cell_Centr->contour.push_back(this->A_Rails[i]->All_point[0]);
		}
		for (int i = 0; i < this->B_Rails.size(); i++)
		{
			Cell_Centr->contour.push_back(this->B_Rails[i]->All_point[0]);
		}
		for (int i = 1; i < this->C_Rails.size(); i++)
		{
			Cell_Centr->contour.push_back(this->C_Rails[i]->All_point[0]);
		}

		auto G2 = new Gran(Cell_Centr->contour[Cell_Centr->contour.size() - 1], Cell_Centr->contour[0], Axis);
		G2->main_gran = true;
		Cell_Centr->Grans.push_back(G2);
		G2->Master = Cell_Centr;
		this->All_Gran.push_back(G2);*/
	}

	for (auto& i : centr_cell)
	{
		this->All_Cells.push_back(i);
	}


	delete Cell_Centr;


	// Нумерация всего
	num = 0;
	for (auto& i : this->All_Gran)  // Нужно для выделения контактной поверхности
	{
		i->number = num;
		num++;
		if (i->A->type == P_Contact && i->B->type == P_Contact)
		{
			this->Line_Contact.push_back(i);
		}
		sort(this->Line_Contact.begin(), this->Line_Contact.end(), [&](Gran* i, Gran* j)
			{
				double x1, y1, x2, y2;
				i->Get_Center(x1, y1);
				j->Get_Center(x2, y2);
				return (polar_angle(x1, y1) < polar_angle(x2, y2));
			});

		if (i->A->type == P_Inner_shock && i->B->type == P_Inner_shock)
		{
			this->Line_Inner.push_back(i);
		}

		sort(this->Line_Inner.begin(), this->Line_Inner.end(), [&](Gran* i, Gran* j)
			{
				double x1, y1, x2, y2;
				i->Get_Center(x1, y1);
				j->Get_Center(x2, y2);
				return (polar_angle(x1, y1) < polar_angle(x2, y2));
			});

		if (i->A->type == P_Outer_shock && i->B->type == P_Outer_shock)
		{
			this->Line_Outer.push_back(i);
		}

		sort(this->Line_Outer.begin(), this->Line_Outer.end(), [&](Gran* i, Gran* j)
			{
				double x1, y1, x2, y2;
				i->Get_Center(x1, y1);
				j->Get_Center(x2, y2);
				return (polar_angle(x1, y1) < polar_angle(x2, y2));
			});
	}

	num = 0;
	for (auto& i : this->All_Gran_copy)  // Нужно для выделения контактной поверхности
	{
		i->number = num;
		num++;
	}

	num = 0;
	for (auto& i : this->All_Cells)
	{
		i->number = num;
		num++;
	}

	num = 0;
	// Нумеруем точки
	for (auto& i : this->All_Points)
	{
		i->number = num;
		num++;
	}


	cout << "Setka.cpp    " << "Vsego tochek = " << this->All_Points.size() << endl;
	cout << "Setka.cpp    " << "Vsego Yacheek = " << this->All_Cells.size() << endl;
}

void Setka::Inizialization(void)
{
	mmu1 = 0.0;
	mmu2 = 0.0;
	mmu3 = 0.0;
	mmu4 = 0.0;
	mmu5 = 0.0;
	mn1 = 0;
	mn2 = 0;
	mn3 = 0;
	mn4 = 0;
	mn5 = 0;
	f_num = 0;
	this->Cell_m = nullptr;
	f_way.open("file_work.txt");

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < I_ + 1; j++)
		{
				Mu_stat[i][j] = 0.0;
		}
	}

	// Считаем веса
	int i1, i2, i3;
	/*ifstream fout;
	fout.open("weyght_1.txt");
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < I_ + 1; j++)
		{
			for (int k = 0; k < J_ + 1; k++)
			{
				fout >> i1 >> i2 >> i3 >> Mu_stat[i][j];
			}
		}
	}*/

	Ri.resize(I_ + 1);
	//Ri[0] = R5_ - 2.0;
	Ri[0] = 1.0/RR_;                 // Здесь задаются радиусы
	Ri[1] = 3.0 / RR_;
	Ri[2] = 10.0 / RR_;
	Ri[3] = 30.0 / RR_;
	Ri[4] = 100.0 / RR_;
	Ri[5] = 250.0 / RR_;
	Ri[6] = 500.0 / RR_;                 // Здесь задаются радиусы
	Ri[7] = 1000.0 / RR_;
	Ri[8] = Rmax_;

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < Al_stat; j++)
		{
			this->mu_mom[i][j] = 0.0;
			this->Vx_mom[i][j] = 0.0;
			this->Vy_mom[i][j] = 0.0;
			this->Vxx_mom[i][j] = 0.0;
			this->Vyy_mom[i][j] = 0.0;
			this->Vxy_mom[i][j] = 0.0;
			this->Vxxx_mom[i][j] = 0.0;
			this->T_mom[i][j] = 0.0;
		}
	}

	this->k_1 = 0;


}

void Setka::Write_file_for_FCMHD(void)
{
	cout << "A" << endl;
	// Нумерация всего
	if (true)
	{
		int num = 1;
		for (auto& i : this->All_Gran)  // Нужно для выделения контактной поверхности
		{
			auto A = i->Sosed;
			if (A == nullptr) A = i->Master;
			if (i->Master->type == Cell_type::C_1 || i->Master->type == Cell_type::C_2 || i->Master->type == Cell_type::C_3 ||
				A->type == Cell_type::C_1 || A->type == Cell_type::C_2 || A->type == Cell_type::C_3)
			{
				i->number = num;
				num++;
			}
			else
			{
				i->number = -1;
			}
		}

		num = 1;
		for (auto& i : this->All_Gran_copy)  // Нужно для выделения контактной поверхности
		{
			i->number = -1;
		}

		num = 1;
		for (auto& i : this->All_Cells)
		{
			if (i->type == Cell_type::C_1 || i->type == Cell_type::C_2 || i->type == Cell_type::C_3)
			{
				i->number = num;
				num++;
			}
			else
			{
				i->number = -1;
			}
		}

		num = 1;
		// Нумеруем точки
		for (auto& i : this->All_Points)
		{
			i->number = num;
			num++;
		}
	}

	cout << "B" << endl;
	ofstream file("FCMHD_1.bin", ios::binary);
	if (!file.is_open())
	{
		cout << "ERROR FCMHD_1.bin" << endl;
		exit(-1);
	}
	cout << "C" << endl;
	int n1 = 0;
	for (const auto& i : this->All_Cells)
	{
		if (i->type == Cell_type::C_1 || i->type == Cell_type::C_2 || i->type == Cell_type::C_3)
		{
			n1++;
		}
	}
	cout << "D" << endl;
	int n2 = 0;
	for (const auto& i : this->All_Gran)
	{
		auto A = i->Sosed;
		if (A == nullptr) A = i->Master;
		if (i->Master->type == Cell_type::C_1 || i->Master->type == Cell_type::C_2 || i->Master->type == Cell_type::C_3 ||
			A->type == Cell_type::C_1 || A->type == Cell_type::C_2 || A->type == Cell_type::C_3)
		{
			n2++;
		}
	}
	
	cout << "E" << endl;
	double time = 0.0;
	file.write(reinterpret_cast<const char*>(&time), sizeof(double));
	file.write(reinterpret_cast<const char*>(&n1), sizeof(int));
	file.write(reinterpret_cast<const char*>(&n2), sizeof(int));

	cout << "F" << endl;
	// host_Cell_par
	for (const auto& i : this->All_Cells)
	{
		if (i->type == Cell_type::C_1 || i->type == Cell_type::C_2 || i->type == Cell_type::C_3)
		{
			double value = i->par[0].ro;
			file.write(reinterpret_cast<const char*>(&value), sizeof(double));

			value = i->par[0].u;
			file.write(reinterpret_cast<const char*>(&value), sizeof(double));

			value = i->par[0].v;
			file.write(reinterpret_cast<const char*>(&value), sizeof(double));

			value = 0.0;
			file.write(reinterpret_cast<const char*>(&value), sizeof(double));

			value = i->par[0].p;
			file.write(reinterpret_cast<const char*>(&value), sizeof(double));

			value = 0.0;
			file.write(reinterpret_cast<const char*>(&value), sizeof(double));

			value = 0.0;
			file.write(reinterpret_cast<const char*>(&value), sizeof(double));

			value = 0.0;
			file.write(reinterpret_cast<const char*>(&value), sizeof(double));

			// Атомы
			value = i->par[0].H_n[2];
			if (value < 0.0)
			{
				cout << "Error   3t4g34tv4vcgv456" << endl;
				cout << i->par[0].H_n[2] << " " << i->par[0].H_u[2] << " " << i->par[0].H_v[2] << " " << i->par[0].H_T[2] << endl;
				exit(-1);
			}
			file.write(reinterpret_cast<const char*>(&value), sizeof(double));

			value = i->par[0].H_u[2];
			file.write(reinterpret_cast<const char*>(&value), sizeof(double));

			value = i->par[0].H_v[2];
			file.write(reinterpret_cast<const char*>(&value), sizeof(double));

			value = i->par[0].H_T[2];
			file.write(reinterpret_cast<const char*>(&value), sizeof(double));

			value = i->par[0].H_n[3];
			if (value < 0.0)
			{
				cout << "Error   rtyhjtr6yh5436y4365" << endl;
				cout << i->par[0].H_n[3] << " " << i->par[0].H_u[3] << " " << i->par[0].H_v[3] << " " << i->par[0].H_T[3] << endl;
				exit(-1);
			}
			file.write(reinterpret_cast<const char*>(&value), sizeof(double));

			value = i->par[0].H_u[3];
			file.write(reinterpret_cast<const char*>(&value), sizeof(double));

			value = i->par[0].H_v[3];
			file.write(reinterpret_cast<const char*>(&value), sizeof(double));

			value = i->par[0].H_T[3];
			file.write(reinterpret_cast<const char*>(&value), sizeof(double));
		}
	}

	cout << "G" << endl;
	// host_Cell_center
	for (const auto& i : this->All_Cells)
	{
		if (i->type == Cell_type::C_1 || i->type == Cell_type::C_2 || i->type == Cell_type::C_3)
		{
			double x, y, z;
			i->Get_Center(x, y);
			z = 0.0;

			file.write(reinterpret_cast<const char*>(&x), sizeof(double));
			file.write(reinterpret_cast<const char*>(&y), sizeof(double));
			file.write(reinterpret_cast<const char*>(&z), sizeof(double));
		}
	}

	cout << "H" << endl;
	// host_Cell_Volume
	for (const auto& i : this->All_Cells)
	{
		if (i->type == Cell_type::C_1 || i->type == Cell_type::C_2 || i->type == Cell_type::C_3)
		{
			double VV;
			VV = i->Get_Volume();

			file.write(reinterpret_cast<const char*>(&VV), sizeof(double));
		}
	}

	cout << "K" << endl;
	// host_Cell_gran
	for (const auto& i : this->All_Cells)
	{
		if (i->type == Cell_type::C_1 || i->type == Cell_type::C_2 || i->type == Cell_type::C_3)
		{
			Gran* A = i->Grans[0];
			if (A->main_gran == false) A = A->Gran_copy;
			int value = A->number;
			if (value == -1)
			{
				cout << "erori uheuyrfguoehfri" << endl;
				exit(-1);
			}
			file.write(reinterpret_cast<const char*>(&value), sizeof(int));

			A = i->Grans[1];
			if (A->main_gran == false) A = A->Gran_copy;
			value = A->number;
			if (value == -1)
			{
				cout << "erori uheuyrfguoehfri" << endl;
				exit(-1);
			}
			file.write(reinterpret_cast<const char*>(&value), sizeof(int));

			A = i->Grans[2];
			if (A->main_gran == false) A = A->Gran_copy;
			value = A->number;
			if (value == -1)
			{
				cout << "erori uheuyrfguoehfri" << endl;
				exit(-1);
			}
			file.write(reinterpret_cast<const char*>(&value), sizeof(int));

			A = i->Grans[3];
			if (A->main_gran == false) A = A->Gran_copy;
			value = A->number;
			if (value == -1)
			{
				cout << "erori uheuyrfguoehfri" << endl;
				exit(-1);
			}
			file.write(reinterpret_cast<const char*>(&value), sizeof(int));

			value = 0;
			file.write(reinterpret_cast<const char*>(&value), sizeof(int));

			value = 0;
			file.write(reinterpret_cast<const char*>(&value), sizeof(int));

		}
	}

	cout << "L" << endl;
	// host_Gran_normal
	for (const auto& i : this->All_Gran)
	{
		auto A = i->Sosed;
		if (A == nullptr) A = i->Master;
		if (i->Master->type == Cell_type::C_1 || i->Master->type == Cell_type::C_2 || i->Master->type == Cell_type::C_3 ||
			A->type == Cell_type::C_1 || A->type == Cell_type::C_2 || A->type == Cell_type::C_3)
		{
			double n1, n2, n3;
			i->Get_normal(n1, n2);
			n3 = 0.0;
			if (i->Master->type != Cell_type::C_1 && i->Master->type != Cell_type::C_2 && i->Master->type != Cell_type::C_3)
			{
				n1 = -n1;
				n2 = -n2;
			}

			file.write(reinterpret_cast<const char*>(&n1), sizeof(double));
			file.write(reinterpret_cast<const char*>(&n2), sizeof(double));
			file.write(reinterpret_cast<const char*>(&n3), sizeof(double));
		}
	}

	cout << "M" << endl;
	// host_Gran_square
	for (const auto& i : this->All_Gran)
	{
		auto A = i->Sosed;
		if (A == nullptr) A = i->Master;
		if (i->Master->type == Cell_type::C_1 || i->Master->type == Cell_type::C_2 || i->Master->type == Cell_type::C_3 ||
			A->type == Cell_type::C_1 || A->type == Cell_type::C_2 || A->type == Cell_type::C_3)
		{
			double value = i->Get_square();
			file.write(reinterpret_cast<const char*>(&value), sizeof(double));
		}
	}

	cout << "N" << endl;
	// host_Gran_center
	for (const auto& i : this->All_Gran)
	{
		auto A = i->Sosed;
		if (A == nullptr) A = i->Master;
		if (i->Master->type == Cell_type::C_1 || i->Master->type == Cell_type::C_2 || i->Master->type == Cell_type::C_3 ||
			A->type == Cell_type::C_1 || A->type == Cell_type::C_2 || A->type == Cell_type::C_3)
		{
			double x, y, z;
			i->Get_Center(x, y);
			z = 0.0;

			file.write(reinterpret_cast<const char*>(&x), sizeof(double));
			file.write(reinterpret_cast<const char*>(&y), sizeof(double));
			file.write(reinterpret_cast<const char*>(&z), sizeof(double));
		}
	}

	cout << "O" << endl;
	// host_Gran_neighbour
	for (const auto& i : this->All_Gran)
	{
		auto A = i->Sosed;
		if (A == nullptr) A = i->Master;
		if (i->Master->type == Cell_type::C_1 || i->Master->type == Cell_type::C_2 || i->Master->type == Cell_type::C_3 ||
			A->type == Cell_type::C_1 || A->type == Cell_type::C_2 || A->type == Cell_type::C_3)
		{
			if (i->Master->type != Cell_type::C_1 && i->Master->type != Cell_type::C_2 && i->Master->type != Cell_type::C_3)
			{
				auto F = i->Master;
				i->Master = i->Sosed;
				i->Sosed = F;
			}

			int value = i->Master->number;
			if (value < 1)
			{
				cout << "Error ejrgkohueghjfe " << endl;
				exit(-1);
			}
			file.write(reinterpret_cast<const char*>(&value), sizeof(int));

			value = -1;
			if(i->Sosed != nullptr) value = i->Sosed->number;
			file.write(reinterpret_cast<const char*>(&value), sizeof(int));
		}
	}

	cout << "P" << endl;
	// host_Gran_neighbour_TVD
	for (const auto& i : this->All_Gran)
	{
		auto A = i->Sosed;
		if (A == nullptr) A = i->Master;
		if (i->Master->type == Cell_type::C_1 || i->Master->type == Cell_type::C_2 || i->Master->type == Cell_type::C_3 ||
			A->type == Cell_type::C_1 || A->type == Cell_type::C_2 || A->type == Cell_type::C_3)
		{
			int value = -1;
			file.write(reinterpret_cast<const char*>(&value), sizeof(int));

			value = -1;
			file.write(reinterpret_cast<const char*>(&value), sizeof(int));
		}
	}

	double vv = 121.0;
	file.write(reinterpret_cast<const char*>(&vv), sizeof(double));

	cout << "Q" << endl;
	// host_Gran_type
	for (const auto& i : this->All_Gran)
	{
		auto A = i->Sosed;
		if (A == nullptr) A = i->Master;
		if (i->Master->type == Cell_type::C_1 || i->Master->type == Cell_type::C_2 || i->Master->type == Cell_type::C_3 ||
			A->type == Cell_type::C_1 || A->type == Cell_type::C_2 || A->type == Cell_type::C_3)
		{
			int value = 1;
			if (i->type == Usualy) value = 2;
			if (i->type == Inner_sphere) value = 3;
			if (i->type == Extern) value = 4;
			if (i->type == Axis) value = 5;
			file.write(reinterpret_cast<const char*>(&value), sizeof(int));
		}
	}

	vv = 122.0;
	file.write(reinterpret_cast<const char*>(&vv), sizeof(double));

	cout << "R" << endl;
	// host_Gran_POTOK
	for (const auto& i : this->All_Gran)
	{
		auto A = i->Sosed;
		if (A == nullptr) A = i->Master;
		if (i->Master->type == Cell_type::C_1 || i->Master->type == Cell_type::C_2 || i->Master->type == Cell_type::C_3 ||
			A->type == Cell_type::C_1 || A->type == Cell_type::C_2 || A->type == Cell_type::C_3)
		{
			for (int j = 0; j < 9; j++)
			{
				double value = 0.0;
				file.write(reinterpret_cast<const char*>(&value), sizeof(double));
			}
		}
	}

	cout << "S" << endl;
	// Проверяющий параметр
	vv = 123.0;
	file.write(reinterpret_cast<const char*>(&vv), sizeof(double));

}

void Setka::Read_file_for_FCMHD(void)
{
	std::ifstream file("FCMHD_1.4_out.bin", std::ios::binary);
	if (!file.is_open()) {
		throw std::runtime_error("Error opening file: FCMHD_1.8_out.bin");
	}

	double Tall = 0.0;
	file.read(reinterpret_cast<char*>(&Tall), sizeof(double));


	// host_Cell_par
	for (const auto& i : this->All_Cells)
	{
		if (i->type == Cell_type::C_1 || i->type == Cell_type::C_2 || i->type == Cell_type::C_3)
		{
			double value;
			file.read(reinterpret_cast<char*>(&value), sizeof(double));
			i->par[0].ro = value;

			file.read(reinterpret_cast<char*>(&value), sizeof(double));
			i->par[0].u = value;

			file.read(reinterpret_cast<char*>(&value), sizeof(double));
			i->par[0].v = value;

			file.read(reinterpret_cast<char*>(&value), sizeof(double));

			file.read(reinterpret_cast<char*>(&value), sizeof(double));
			i->par[0].p = value;

			file.read(reinterpret_cast<char*>(&value), sizeof(double));
			file.read(reinterpret_cast<char*>(&value), sizeof(double));
			file.read(reinterpret_cast<char*>(&value), sizeof(double));

			// атомы
			file.read(reinterpret_cast<char*>(&value), sizeof(double));
			file.read(reinterpret_cast<char*>(&value), sizeof(double));
			file.read(reinterpret_cast<char*>(&value), sizeof(double));
			file.read(reinterpret_cast<char*>(&value), sizeof(double));

			file.read(reinterpret_cast<char*>(&value), sizeof(double));
			file.read(reinterpret_cast<char*>(&value), sizeof(double));
			file.read(reinterpret_cast<char*>(&value), sizeof(double));
			file.read(reinterpret_cast<char*>(&value), sizeof(double));
		}
	}


	// Проверяющий параметр
	double vv;
	file.read(reinterpret_cast<char*>(&vv), sizeof(double));

	cout << "Proverka (321) = " << vv << endl;

	file.close();
}

void Pereinterpol(Setka* S1, Setka* S2)
{
	for (const auto& i : S2->All_Cells)
	{
		double x, y;
		int b;
		i->Get_Center(x, y);
		Cell* A = S1->Find_cell(b, x, y);
		if (b == 0)
		{
			cout << "error iergiyhgoiergger  " << x << "   " << y << endl;
			i->par[0].ro = 0.00001;
			i->par[0].p = 0.00001;
			i->par[0].u = 0.0;
			i->par[0].v = 0.0;
			i->par[0].H_n[2] = 0.0;
			i->par[0].H_n[3] = 0.0;
			//exit(-1);
		}
		else
		{
			i->par[0] = A->par[0];
		}
	}
}

void Setka::Save_surface(string name)
{
	double x, y, alp, r;
	ofstream out(name, ios::binary|ios::out);
	//ofstream out2("testik.txt");
	int n = this->Line_Inner.size();
	cout << "Save_surface START   n = " << n << endl;

	out.write((char*)&n, sizeof n); //Записываем в файл число 

	for (int i = 0; i < n; i++)
	{
		this->Line_Inner[i]->Get_Center(x, y);
		//out2 << x << " " << y << endl;
		r = sqrt(x * x + y * y);
		alp = polar_angle(x, y);
		out.write((char*)&alp, sizeof alp);
		out.write((char*)&r, sizeof r);
	}

	n = this->Line_Contact.size();
	out.write((char*)&n, sizeof n); //Записываем в файл число 
	for (int i = 0; i < n; i++)
	{
		this->Line_Contact[i]->Get_Center(x, y);
		r = sqrt(x * x + y * y);
		alp = polar_angle(x, y);

		out.write((char*)&x, sizeof x);
		out.write((char*)&y, sizeof y);
		out.write((char*)&alp, sizeof alp);
		out.write((char*)&r, sizeof r);
	}

	n = this->Line_Outer.size();
	out.write((char*)&n, sizeof n); //Записываем в файл число 
	for (int i = 0; i < n; i++)
	{
		this->Line_Outer[i]->Get_Center(x, y);
		r = sqrt(x * x + y * y);
		alp = polar_angle(x, y);

		out.write((char*)&x, sizeof x);
		out.write((char*)&y, sizeof y);
		out.write((char*)&alp, sizeof alp);
		out.write((char*)&r, sizeof r);
	}

	out.close(); //Закрываем файл
	//out2.close();

	cout << "Save_surface END" << endl;
}

void Setka::Download_surface(string name)
{
	ifstream fin(name, ios::binary | ios::in);
	if (!fin)
	{
		cout << "Net fail 678ijhgftyuwgfghjiru" << endl;
		exit(-1);
	}
	int n1, n2, n3;
	double x, y, al, r;

	fin.read((char*)&n1, sizeof n1); // Считали количество точек на TS
	TS_n = n1;
	TS_al = new double[n1];
	TS_r = new double[n1];
	for (int i = 0; i < n1; i++)
	{
		fin.read((char*)&TS_al[i], sizeof x);
		fin.read((char*)&TS_r[i], sizeof x);
	}

	fin.read((char*)&n2, sizeof n2); // Считали количество точек на HP
	HP_n = n2;
	HP_al = new double[n2];
	HP_r = new double[n2];
	HP_x = new double[n2];
	HP_y = new double[n2];
	for (int i = 0; i < n2; i++)
	{
		fin.read((char*)&HP_x[i], sizeof x);
		fin.read((char*)&HP_y[i], sizeof x);
		fin.read((char*)&HP_al[i], sizeof x);
		fin.read((char*)&HP_r[i], sizeof x);
	}


	fin.read((char*)&n3, sizeof n3); // Считали количество точек на BS
	BS_n = n3;
	BS_al = new double[n3];
	BS_r = new double[n3];
	BS_x = new double[n3];
	BS_y = new double[n3];
	for (int i = 0; i < n3; i++)
	{
		fin.read((char*)&BS_x[i], sizeof x);
		fin.read((char*)&BS_y[i], sizeof x);
		fin.read((char*)&BS_al[i], sizeof x);
		fin.read((char*)&BS_r[i], sizeof x);
	}

	for (auto& i : this->All_Points)
	{
		i->x2 = i->x;
		i->y2 = i->y;
	}

	for (int jj = 0; jj < this->A_Rails.size(); jj++)
	{
		auto i = this->A_Rails[jj];

		x = i->Key_point[0]->x;
		y = i->Key_point[0]->y;

		al = polar_angle(x, y);
		r = find_TS(al);
		i->Key_point[0]->x2 = r * cos(al);
		i->Key_point[0]->y2 = r * sin(al);

		r = find_HP_angle(al);
		i->Key_point[1]->x2 = r * cos(al);
		i->Key_point[1]->y2 = r * sin(al);

		r = find_BS_angle(al);
		i->Key_point[2]->x2 = r * cos(al);
		i->Key_point[2]->y2 = r * sin(al);
	}

	for (int jj = 0; jj < this->B_Rails.size(); jj++)
	{
		auto i = this->B_Rails[jj];
		// Подвинем ключевые точки   0 - внутренняя ударная волна

		x = i->Key_point[0]->x;
		y = i->Key_point[0]->y;

		al = polar_angle(x, y);
		r = find_TS(al);
		i->Key_point[0]->x2 = r * cos(al);
		i->Key_point[0]->y2 = r * sin(al);

		x = i->Key_point[0]->x2;
		y = find_HP_x(x);
		i->Key_point[1]->x2 = x;
		i->Key_point[1]->y2 = y;

		x = i->Key_point[0]->x2;
		y = find_BS_x(x);
		i->Key_point[2]->x2 = x;
		i->Key_point[2]->y2 = y;
	}

	for (int jj = 0; jj < this->C_Rails.size(); jj++)
	{
		auto i = this->C_Rails[jj];

		x = i->Key_point[0]->x;
		y = i->Key_point[0]->y;

		al = polar_angle(x, y);
		r = find_TS(al);

		if (jj != 0) // Этот узел уже был подвинут в B - лучах
		{
			i->Key_point[0]->x2 = r * cos(al);
			i->Key_point[0]->y2 = r * sin(al);
		}
	}

	Move_Setka_Calculate_stat();

	for (auto& i : this->All_Points)
	{
		i->x = i->x2;
		i->y = i->y2;
	}

	for (int jj = 0; jj < this->D_Rails.size(); jj++)
	{
		auto i = this->D_Rails[jj];

		i->Key_point[1]->x2 = i->Key_point[0]->x2;

		x = i->Key_point[1]->x2;
		y = find_HP_x(x);
		i->Key_point[1]->y2 = y;

		i->Key_point[2]->x2 = i->Key_point[0]->x2;
		x = i->Key_point[2]->x2;
		y = find_BS_x(x);
		i->Key_point[2]->y2 = y;


		if (jj >= n_outer_shock)                                                        // ТУТ ИЗМЕНИЛ
		{
			i->Key_point[2]->y2 = this->D_Rails[jj - 1]->Key_point[2]->y2;
		}

	}


	Move_Setka_Calculate_stat();

	for (auto& i : this->All_Points)
	{
		i->x = i->x2;
		i->y = i->y2;
	}

	delete[] BS_al;
	delete[] BS_r;
	delete[] BS_x;
	delete[] BS_y;
	delete[] HP_al;
	delete[] HP_r;
	delete[] HP_x;
	delete[] HP_y;
	delete[] TS_al;
	delete[] TS_r;
}

double Setka::find_TS(const double& alp)
{
	// Получить r_TS по углу 
	double al1, al2, r1, r2;
	int n = this->TS_n;

	if (TS_al[0] >= alp)
	{
		return TS_r[0];
	}
	else if (TS_al[n - 1] <= alp)
	{
		return TS_r[n - 1];
	}
	else
	{
		for (int i = 0; i < n; i++)
		{
			if (alp < TS_al[i])
			{
				al1 = TS_al[i - 1];
				al2 = TS_al[i];
				r1 = TS_r[i - 1];
				r2 = TS_r[i];
				break;
			}
		}
		return linear(al1, r1, al2, r2, alp);
	}

	cout << "Error 1333 yuikngtrdcvhjiolmnhgfvb" << endl;
	return 0.0;
}

double Setka::find_HP_angle(const double& alp)
{
	// Получить r_TS по углу 
	double al1, al2, r1, r2;
	int n = this->HP_n;

	if (HP_al[0] >= alp)
	{
		return HP_r[0];
	}
	else if (HP_al[n - 1] <= alp)
	{
		return HP_r[n - 1];
	}
	else
	{
		for (int i = 0; i < n; i++)
		{
			if (alp < HP_al[i])
			{
				al1 = HP_al[i - 1];
				al2 = HP_al[i];
				r1 = HP_r[i - 1];
				r2 = HP_r[i];
				break;
			}
		}
		return linear(al1, r1, al2, r2, alp);
	}

	cout << "Error 1369 iuytgfghjkoplknbvcxzaqw3edsx" << endl;
	return 0.0;
}

double Setka::find_HP_x(const double& x)
{
	// Получить r_TS по углу 
	double x1, x2, y1, y2;
	int n = this->HP_n;

	if (HP_x[0] <= x)
	{
		return HP_y[0];
	}
	else if (HP_x[n - 1] >= x)
	{
		return HP_y[n - 1];
	}
	else
	{
		for (int i = 0; i < n; i++)
		{
			if (x > HP_x[i])
			{
				x1 = HP_x[i - 1];
				x2 = HP_x[i];
				y1 = HP_y[i - 1];
				y2 = HP_y[i];
				break;
			}
		}
		return linear(x1, y1, x2, y2, x);
	}

	cout << "Error 1420 234esfgytredghuyt56" << endl;
	return 0.0;
}

double Setka::find_BS_x(const double& x)
{
	// Получить r_TS по углу 
	double x1, x2, y1, y2;
	int n = this->BS_n;

	if (BS_x[0] <= x)
	{
		return BS_y[0];
	}
	else if (BS_x[n - 1] >= x)
	{
		return BS_y[n - 1];
	}
	else
	{
		for (int i = 0; i < n; i++)
		{
			if (x > BS_x[i])
			{
				x1 = BS_x[i - 1];
				x2 = BS_x[i];
				y1 = BS_y[i - 1];
				y2 = BS_y[i];
				break;
			}
		}
		return linear(x1, y1, x2, y2, x);
	}

	cout << "Error 1420 234esfgytredghuyt56" << endl;
	return 0.0;
}

double Setka::find_BS_angle(const double& alp)
{
	// Получить r_TS по углу 
	double al1, al2, r1, r2;
	int n = this->BS_n;

	if (BS_al[0] >= alp)
	{
		return BS_r[0];
	}
	else if (BS_al[n - 1] <= alp)
	{
		return BS_r[n - 1];
	}
	else
	{
		for (int i = 0; i < n; i++)
		{
			if (alp < BS_al[i])
			{
				al1 = BS_al[i - 1];
				al2 = BS_al[i];
				r1 = BS_r[i - 1];
				r2 = BS_r[i];
				break;
			}
		}
		return linear(al1, r1, al2, r2, alp);
	}

	cout << "Error 1369 iuytgfghjkoplknbvcxzaqw3edsx" << endl;
	return 0.0;
}

void Setka::normir(int ii)
{
	if (ii == 0)// УМЕНЬШАЕМ СКОРОСТЬ
	{
		for (auto i : this->All_Cells)
		{
			if (i->par[0].Q/ i->par[0].ro < 90.0)//(i->type == C_centr || i->type == C_1 || i->type == C_2 || i->type == C_3)
			{
				// УМЕНЬШАЕМ СКОРОСТЬ
				i->par[0].u = i->par[0].u / (chi_real / chi_);       // Перенормировка
				i->par[0].v = i->par[0].v / (chi_real / chi_);
				i->par[0].ro = i->par[0].ro * kv(chi_real / chi_);
				i->par[0].Q = i->par[0].Q * kv(chi_real / chi_);

				i->par[1].u = i->par[1].u / (chi_real / chi_);       // Перенормировка
				i->par[1].v = i->par[1].v / (chi_real / chi_);
				i->par[1].ro = i->par[1].ro * kv(chi_real / chi_);
				i->par[1].Q = i->par[1].Q * kv(chi_real / chi_);
			}
		}
	}
	else
	{
		for (auto i : this->All_Cells)
		{
			if (i->par[0].Q / i->par[0].ro < 90.0)//(i->type == C_centr || i->type == C_1 || i->type == C_2 || i->type == C_3)
			{
				i->par[0].u = i->par[0].u * (chi_real / chi_);       // Перенормировка
				i->par[0].v = i->par[0].v * (chi_real / chi_);
				i->par[0].ro = i->par[0].ro / kv(chi_real / chi_);
				i->par[0].Q = i->par[0].Q / kv(chi_real / chi_);

				i->par[1].u = i->par[1].u * (chi_real / chi_);       // Перенормировка
				i->par[1].v = i->par[1].v * (chi_real / chi_);
				i->par[1].ro = i->par[1].ro / kv(chi_real / chi_);
				i->par[1].Q = i->par[1].Q / kv(chi_real / chi_);
			}
		}
	}
}

//void Setka::TVD_prepare(void)
//{
//	double n1, n2, n;
//	double x, y;
//	double x1, y1;
//	double x3, y3;
//	double a, b, N;
//	Cell* A = nullptr;
//	cout << "Setka.cpp    " << "TVD  prepare" << endl;
//	for (auto& i : this->All_Gran)
//	{
//		N = 1.0;
//		A = nullptr;
//		if (i->type == Axis)
//		{
//			i->Get_Center(x, y);
//			i->Get_normal(n1, n2);
//			n1 = -n1;
//			n2 = -n2;
//			for (auto& j : i->Master->Grans)
//			{
//				if (j->type == Usualy)
//				{
//					j->Sosed->Get_Center(x1, y1);
//
//					a = x1 - x;
//					b = y1 - y;
//					n = sqrt(kv(a) + kv(b));
//					a = a / n;
//					b = b / n;
//					//cout << "Setka.cpp    " << fabs(a * n1 + b * n2) << endl;
//					if (fabs(1.0 - fabs(a * n1 + b * n2)) < N)
//					{
//						N = fabs(1.0 - fabs(a * n1 + b * n2));
//						A = j->Sosed;
//					}
//				}
//
//			}
//			i->Sosed_down = A;
//			continue;
//		}
//		if (i->type != Usualy)
//		{
//			continue;
//		}
//		i->Sosed->Get_Center(x3, y3);
//		i->Get_Center(x, y);
//		n1 = x3 - x;                            // Направление от центра грани к центру соседа для грани
//		n2 = y3 - y;
//		n = sqrt(kv(n1) + kv(n2));
//		n1 = n1 / n;
//		n2 = n2 / n;
//
//		for (auto& j : i->Sosed->Grans)
//		{
//			if (j->type != Usualy || j->Sosed->number != i->Master->number)
//			{
//				if (j->Sosed != nullptr)
//				{
//					j->Sosed->Get_Center(x1, y1);
//				}
//				else
//				{
//					j->Get_Center(x1, y1);
//				}
//
//				a = x1 - x;
//				b = y1 - y;
//				n = sqrt(kv(a) + kv(b));
//				a = a / n;
//				b = b / n;
//				//cout << "Setka.cpp    " << fabs(a * n1 + b * n2) << endl;
//				if (fabs(1.0 - fabs(a * n1 + b * n2)) < N)
//				{
//					N = fabs(1.0 - fabs(a * n1 + b * n2));
//					A = j->Sosed;
//				}
//			}
//			
//		}
//		i->Sosed_up = A;
//	}
//
//	for (auto& i : this->All_Gran_copy)
//	{
//		A = nullptr;
//		if (i->type != Usualy)
//		{
//			continue;
//		}
//		N = 1.0;
//		i->Sosed->Get_Center(x3, y3);
//		i->Get_Center(x, y);
//		n1 = x3 - x;
//		n2 = y3 - y;
//		n = sqrt(kv(n1) + kv(n2));
//		n1 = n1 / n;
//		n2 = n2 / n;
//		for (auto& j : i->Sosed->Grans)
//		{
//			if (j->type != Usualy || j->Sosed->number != i->Master->number)
//			{
//				if (j->Sosed != nullptr)
//				{
//					j->Sosed->Get_Center(x1, y1);
//				}
//				else
//				{
//					j->Get_Center(x1, y1);
//				}
//
//				a = x1 - x;
//				b = y1 - y;
//				n = sqrt(kv(a) + kv(b));
//				a = a / n;
//				b = b / n;
//				//cout << "Setka.cpp    " << fabs(a * n1 + b * n2) << endl;
//				if (fabs(1.0 - fabs(a * n1 + b * n2)) < N)
//				{
//					N = fabs(1.0 - fabs(a * n1 + b * n2));
//					A = j->Sosed;
//				}
//			}
//
//		}
//		i->Sosed_up = A;
//	}
//
//	for (auto& i : this->All_Gran)
//	{
//		if (i->Gran_copy != nullptr)
//		{
//			i->Sosed_down = i->Gran_copy->Sosed_up;
//		}
//	}
//
//	for (auto& i : this->All_Gran_copy)
//	{
//		if (i->Gran_copy != nullptr)
//		{
//			i->Sosed_down = i->Gran_copy->Sosed_up;
//		}
//	}
//
//	cout << "Setka.cpp    " << "TVD  prepare  end" << endl;
//}

void Setka::TVD_prepare(void)
{
	double n1, n2, n;
	double x, y;
	double x1, y1;
	double x3, y3;
	double a, b, N;
	Cell* A = nullptr;
	Gran* B = nullptr;
	//cout << "Setka.cpp    " << "TVD  prepare" << endl;

	// Находим Sosed_up для этой грани
	for (auto& i : this->All_Gran)
	{
		N = -1.0;
		B = nullptr;
		if (i->Sosed == nullptr)
		{
			i->Sosed_up = nullptr;
			continue;
		}

		i->Sosed->Get_Center(x3, y3);
		i->Get_Center(x, y);
		n1 = x3 - x;                            // Направление от центра грани к центру соседа для грани
		n2 = y3 - y;
		n = sqrt(kv(n1) + kv(n2));
		n1 = n1 / n;
		n2 = n2 / n;

		for (auto& j : i->Sosed->Grans)
		{
			if (j->type == Usualy)
			{
				if (j->Gran_copy->number == i->number)
				{
					continue;
				}
			}

			j->Get_Center(x1, y1);
			a = x1 - x3;
			b = y1 - y3;
			n = sqrt(kv(a) + kv(b));
			a = a / n;
			b = b / n;
			if (a * n1 + b * n2 > N)
			{
				N = a * n1 + b * n2;
				B = j;
			}
		}

		if (B->Sosed != nullptr)
		{
			i->Sosed_up = B->Sosed;
		}
		else
		{
			i->Sosed_up = nullptr;
		}
	}

	for (auto& i : this->All_Gran_copy)
	{
		N = -1.0;
		B = nullptr;
		if (i->Sosed == nullptr)
		{
			if (i->type == Axis)
			{
				cout << "Setka.cpp    " << "Axis " << endl;
			}
			i->Sosed_up = nullptr;
			continue;
		}

		i->Sosed->Get_Center(x3, y3);
		i->Get_Center(x, y);
		n1 = x3 - x;                            // Направление от центра грани к центру соседа для грани
		n2 = y3 - y;
		n = sqrt(kv(n1) + kv(n2));
		n1 = n1 / n;
		n2 = n2 / n;

		for (auto& j : i->Sosed->Grans)
		{
			if (j->type == Usualy)
			{
				if (j->Gran_copy->number == i->number)
				{
					continue;
				}
			}

			j->Get_Center(x1, y1);
			a = x1 - x3;
			b = y1 - y3;
			n = sqrt(kv(a) + kv(b));
			a = a / n;
			b = b / n;
			if (a * n1 + b * n2 > N)
			{
				N = a * n1 + b * n2;
				B = j;
			}
		}

		if (B->Sosed != nullptr)
		{
			i->Sosed_up = B->Sosed;
		}
		else
		{
			i->Sosed_up = nullptr;
		}
	}

	for (auto& i : this->All_Gran)
	{
		if (i->Gran_copy != nullptr)
		{
			i->Sosed_down = i->Gran_copy->Sosed_up;
		}
		else
		{
			i->Sosed_down = nullptr;
		}
	}

	for (auto& i : this->All_Gran_copy)
	{
		if (i->Gran_copy != nullptr)
		{
			i->Sosed_down = i->Gran_copy->Sosed_up;
		}
		else
		{
			i->Sosed_down = nullptr;
		}
	}

	//cout << "Setka.cpp    " << "TVD  prepare  end" << endl;

	// Подготовим всё для переинтерполяции источников
	for (auto& i : this->All_Cells)
	{
		bool b1 = false;
		bool b2 = false;
		double xx, yy;
		double a1, a2;
		i->Get_Center(xx, yy);
		i->r_istoch = sqrt(kv(xx) + kv(yy));
		if (xx >= 0.0)
		{
			a1 = xx / i->r_istoch;
			a2 = yy / i->r_istoch;

			for (auto& j : i->Grans)
			{
				if (j->Sosed == nullptr)
				{
					continue;
				}
				double xx2, yy2;
				j->Sosed->Get_Center(xx2, yy2);
				double n1 = (xx2 - xx);
				double n2 = (yy2 - yy);
				double nn = sqrt(kv(n1) + kv(n2));
				n1 = n1 / nn;
				n2 = n2 / nn;
				double sk = a1 * n1 + a2 * n2;
				if (fabs(1.0 - sk) < 0.001)
				{
					i->Next = j->Sosed;
					b1 = true;
					if (sqrt(kv(xx2) + kv(yy2)) < i->r_istoch)
					{
						cout << "EROR 1459 vdvdftghdvdscsa " << endl;
						exit(-1);
					}
				}
				else if (fabs(-1.0 - sk) < 0.001)
				{
					i->Back = j->Sosed;
					b2 = true;
					if (sqrt(kv(xx2) + kv(yy2)) > i->r_istoch)
					{
						cout << "EROR 1467 vdvdftghdvdscsa " << endl;
						exit(-1);
					}
				}
			}

		}
	}

	//cout << "Setka.cpp    " << "TVD  prepare  end 2" << endl;
}

void Setka::Print_TVD(void)
{
	int ll = this->All_Gran.size();
	ofstream fout;
	fout.open("TVD.txt");
	fout << "TITLE = \"HP\" ";
	fout << " VARIABLES = \"X\", \"Y\", \"C\"  ZONE T= \"HP\", N=  " << 2 * ll;
	fout << " , E= " << ll;
	fout << " , F=FEPOINT, ET=LINESEG  " << endl;
	double x1, y1;
	double x2, y2;
	int nn = 0;
	for (auto& i : this->All_Gran)
	{
		nn++;
		i->Get_Center(x1, y1);
		if (i->Sosed_down != nullptr)
		{
			i->Sosed_down->Get_Center(x2, y2);
		}
		else
		{
			x2 = x1;
			y2 = y1;
		}
		fout << x1 << " " << y1 << " " << nn << endl;
		fout << x2 << " " << y2 << " " << nn << endl;
	}

	for (int i = 0; i < ll; i++)
	{
		fout << 2 * i + 1 << " " << 2 * i + 2 << endl;
	}

	fout.close();
}

void Setka::Print_for_Igor(void)
{
	ofstream fout;
	string name_fff = "print_for_Igor_2.txt";
	fout.open(name_fff);

	// Считаем дивергенцию
	for (int num_cell = 0; num_cell < this->All_Cells.size(); num_cell++)
	{
		//cout << num_cell << endl;
		auto K = this->All_Cells[num_cell];

		Parametr par1 = K->par[0];

		double dist;
		double x, y;
		K->Get_Center(x, y);


		double radius = sqrt(kv(x) + kv(y));


		double Volume = K->Get_Volume();
		bool np = true;
		double Vdiv = 0.0;

		for (auto& i : K->Grans)
		{
			double x2, y2, x4, y4;
			double S = i->Get_square();
			double n1, n2;
			i->Get_normal(n1, n2);
			i->Get_Center(x2, y2);
			i->Get_Center_posle(x4, y4);
			dist = sqrt(kv(x - x2) + kv(y - y2));
			Parametr par2, par11;
			par11 = par1;
			i->Get_par(par2, 0);
			double xx2, yy2;
			double dist2 = dist;

			if (i->Sosed != nullptr)
			{
				i->Sosed->Get_Center(xx2, yy2);
				dist2 = sqrt(kv(x2 - xx2) + kv(y2 - yy2));
				par2 = i->Sosed->par[0];
			}
			else
			{
				par2 = par1;
				if (i->type == Axis)
				{
					par2.v = -par2.v;
				}
			}

			if (i->Sosed != nullptr)
			{
				if (K->type == C_1 && i->Sosed->type == C_2)
				{
					par2 = par1;
				}
				if (K->type == C_2 && i->Sosed->type == C_1)
				{
					par2 = par1;
					cout << "2  ->  1  " << x << " " << y << endl;
					double kx, ky;
					i->Sosed_down->Get_Center(kx, ky);
					double dist3 = sqrt(kv(x2 - kx) + kv(y2 - ky));
					cout << kx << "  " << ky << endl;
					par2.u = linear(-dist3, i->Sosed_down->par[0].u, -dist, par1.u, dist2);
					par2.v = linear(-dist3, i->Sosed_down->par[0].v, -dist, par1.v, dist2);
				}
				if (K->type == C_1 && i->Sosed->type == C_3)
				{
					par2 = par1;
				}
				if (K->type == C_3 && i->Sosed->type == C_1)
				{
					par2 = par1;

					double kx, ky;
					i->Sosed_down->Get_Center(kx, ky);
					double dist3 = sqrt(kv(x2 - kx) + kv(y2 - ky));
					par2.u = linear(-dist3, i->Sosed_down->par[0].u, -dist, par1.u, dist2);
					par2.v = linear(-dist3, i->Sosed_down->par[0].v, -dist, par1.v, dist2);
				}
				if (K->type == C_3 && i->Sosed->type == C_4)
				{
					par2 = par1;
					double kx, ky;
					i->Sosed_down->Get_Center(kx, ky);
					double dist3 = sqrt(kv(x2 - kx) + kv(y2 - ky));
					par2.u = linear(-dist3, i->Sosed_down->par[0].u, -dist, par1.u, dist2);
					par2.v = linear(-dist3, i->Sosed_down->par[0].v, -dist, par1.v, dist2);
				}
				if (K->type == C_4 && i->Sosed->type == C_3)
				{
					par2 = par1;
				}
			}

			double sks = n1 * (par11.u * dist2 + par2.u * dist) / (dist + dist2) + n2 * (par11.v * dist2 + par2.v * dist) / (dist + dist2);

			Vdiv = Vdiv + S * sks;

		}


		double v = par1.v;

		K->par[0].divV = v / y + Vdiv / Volume;
	}

	cout << "Divergence " << endl;

	for (auto& i : this->All_Cells)
	{
		double x, y;
		i->Get_Center(x, y);
		fout << x << " " << y << //
			" " << i->par[0].ro << " " << i->par[0].p << " " //
			<< i->par[0].u << " " << i->par[0].v << " " << i->par[0].divV << " " <<//
			i->par[0].H_n[0] << " " << i->par[0].H_u[0] << " " << i->par[0].H_v[0] << " " << i->par[0].H_T[0] << " " << i->par[0].H_uu[0] << " " << i->par[0].H_uv[0] << " " << i->par[0].H_vv[0] << " " << i->par[0].H_uuu[0] << " " //
			<< i->par[0].H_n[1] << " " << i->par[0].H_u[1] << " " << i->par[0].H_v[1] << " " << i->par[0].H_T[1] << " " << i->par[0].H_uu[1] << " " << i->par[0].H_uv[1] << " " << i->par[0].H_vv[1] << " " << i->par[0].H_uuu[1] << " " //
			<< i->par[0].H_n[2] << " " << i->par[0].H_u[2] << " " << i->par[0].H_v[2] << " " << i->par[0].H_T[2] << " " << i->par[0].H_uu[2] << " " << i->par[0].H_uv[2] << " " << i->par[0].H_vv[2] << " " << i->par[0].H_uuu[2] << " " //
			<< i->par[0].H_n[3] << " " << i->par[0].H_u[3] << " " << i->par[0].H_v[3] << " " << i->par[0].H_T[3] << " " << i->par[0].H_uu[3] << " " << i->par[0].H_uv[3] << " " << i->par[0].H_vv[3] << " " << i->par[0].H_uuu[3] << " " //
			<< i->par[0].H_n[4] << " " << i->par[0].H_u[4] << " " << i->par[0].H_v[4] << " " << i->par[0].H_T[4] << " " << i->par[0].H_uu[4] << " " << i->par[0].H_uv[4] << " " << i->par[0].H_vv[4] << " " << i->par[0].H_uuu[4] << " " //
			<< i->par[0].H_n[5] << " " << i->par[0].H_u[5] << " " << i->par[0].H_v[5] << " " << i->par[0].H_T[5] << " " << i->par[0].H_uu[5] << " " << i->par[0].H_uv[5] << " " << i->par[0].H_vv[5] << " " << i->par[0].H_uuu[5] << endl;
		for (int ii = 0; ii < n_S - 1; ii++)
		{
			fout << i->S_m[ii] << " ";
		}
		fout << i->S_m[n_S - 1] << endl;

		for (int ii = 0; ii < n_S - 1; ii++)
		{
			fout << i->S_p[ii] << " ";
		}
		fout << i->S_p[n_S - 1] << endl;
	}
	fout.close();


	fout.open("TS.txt");
	for (auto& i : this->Line_Inner)
	{
		double x, y;
		double phi, r;
		x = (i->A->x + i->B->x) / 2.0;
		y = (i->A->y + i->B->y) / 2.0;
		phi = polar_angle(x, y);
		r = sqrt(x * x + y * y);
		fout << phi << " " << r << " " << i->Master->par[0].ro << " " << i->Sosed->par[0].ro << endl;
	}
	fout.close();

	fout.open("HP.txt");
	for (auto& i : this->Line_Contact)
	{
		double x, y;
		double phi, r;
		x = (i->A->x + i->B->x) / 2.0;
		y = (i->A->y + i->B->y) / 2.0;
		phi = polar_angle(x, y);
		r = sqrt(x * x + y * y);
		fout << phi << " " << r << " " << i->Master->par[0].ro << " " << i->Sosed->par[0].ro << endl;
	}
	fout.close();
}

void Setka::Print_point()
{
	ofstream fout;
	fout.open("point.txt");

	for (auto& i : this->All_Points)
	{
		fout << i->x << " " << i->y << " " << i->type << endl;
	}
}

void Setka::Print_Gran()
{
	int ll = this->Line_Contact.size() + this->Line_Inner.size() + this->Line_Outer.size();
	ofstream fout;
	fout.open("surfaces.txt");
	fout << "TITLE = \"HP\" ";
	fout << " VARIABLES = \"X\", \"Y\"  ZONE T= \"HP\", N=  " << 2 * ll;
	fout << " , E= " << ll;
	fout << " , F=FEPOINT, ET=LINESEG  " << endl;
	double x1, y1;
	double LL = 1.0;
	for (auto& i : this->Line_Contact)
	{
		fout << i->A->x/LL << " " << i->A->y / LL << endl;
		fout << i->B->x / LL << " " << i->B->y / LL << endl;
	}
	for (auto& i : this->Line_Inner)
	{
		fout << i->A->x / LL << " " << i->A->y / LL << endl;
		fout << i->B->x / LL << " " << i->B->y / LL << endl;
	}
	for (auto& i : this->Line_Outer)
	{
		fout << i->A->x / LL << " " << i->A->y / LL << endl;
		fout << i->B->x / LL << " " << i->B->y / LL << endl;
	}
	for (int i = 0; i < ll; i++)
	{
		fout << 2 * i + 1 << " " << 2 * i + 2 << endl;
	}

	fout.close();
}

void Setka::Print_Gran(string name)
{
	int ll = this->Line_Contact.size() + this->Line_Inner.size() + this->Line_Outer.size();
	ofstream fout;
	fout.open(name);
	ofstream fout2;
	fout2.open("position_" + name);
	fout << "TITLE = \"HP\" ";
	fout << " VARIABLES = \"X\", \"Y\"  ZONE T= \"HP\", N=  " << 2 * ll;
	fout << " , E= " << ll;
	fout << " , F=FEPOINT, ET=LINESEG  " << endl;
	double x1, y1;
	double LL = 1.0;
	double aa, bb;
	int kk = 0;
	for (auto& i : this->Line_Contact)
	{
		fout << i->A->x / LL << " " << i->A->y / LL << endl;
		fout << i->B->x / LL << " " << i->B->y / LL << endl;
		if (kk == 0)
		{
			kk = 1;
			aa = i->A->x / LL;
			bb = i->A->y / LL;
		}
	}

	fout2 << aa << " " << bb << endl;
	kk = 0;
	for (auto& i : this->Line_Inner)
	{
		fout << i->A->x / LL << " " << i->A->y / LL << endl;
		fout << i->B->x / LL << " " << i->B->y / LL << endl;
		if (kk == 0)
		{
			kk = 1;
			aa = i->A->x / LL;
			bb = i->A->y / LL;
		}
	}

	fout2 << aa << " " << bb << endl;
	kk = 0;
	for (auto& i : this->Line_Outer)
	{
		fout << i->A->x / LL << " " << i->A->y / LL << endl;
		fout << i->B->x / LL << " " << i->B->y / LL << endl;
		if (kk == 0)
		{
			kk = 1;
			aa = i->A->x / LL;
			bb = i->A->y / LL;
		}
	}

	fout2 << aa << " " << bb << endl;

	for (int i = 0; i < ll; i++)
	{
		fout << 2 * i + 1 << " " << 2 * i + 2 << endl;
	}

	fout.close();
}

void Setka::Print_cell_type(void)
{
	ofstream fout;
	fout.open("Cell_type.txt");
	double x, y;
	for (auto& i : this->All_Cells)
	{
		i->Get_Center(x, y);
		fout << x << " " << y << " " << i->type << endl;
	}
}

void Setka::Print_Gran_type(void)
{
	ofstream fout;
	fout.open("Gran_type.txt");
	double x, y;
	for (auto& i : this->All_Gran)
	{
		i->Get_Center(x, y);
		fout << x << " " << y << " " << i->type << endl;
	}
}

void Setka::Print_cell(void)
{
	int ll = this->All_Cells.size();
	ofstream fout;
	fout.open("Setka_all_cell.txt");
	fout << "TITLE = \"HP\" ";
	fout << " VARIABLES = \"X\", \"Y\"  ZONE T= \"HP\", N=  " << 4 * ll;
	fout << " , E= " << 4 * ll;
	fout << " , F=FEPOINT, ET=LINESEG  " << endl;
	for (auto& i : this->All_Cells)
	{
		for (auto& j : i->contour)
		{
			fout << j->x << " " << j->y << endl;
		}
	}
	for (int i = 0; i < ll; i++)
	{
		fout << 4 * i + 1 << " " << 4 * i + 2 << endl;
		fout << 4 * i + 2 << " " << 4 * i + 3 << endl;
		fout << 4 * i + 3 << " " << 4 * i + 4 << endl;
		fout << 4 * i + 4 << " " << 4 * i + 1 << endl;
	}

	fout.close();
}

void Setka::Print_cell2(string name)
{
	int ll = 0;
	for (auto& i : this->All_Cells)
	{
		ll += i->contour.size();
	}

	ofstream fout;
	fout.open(name);
	fout << "TITLE = \"HP\" ";
	fout << " VARIABLES = \"X\", \"Y\"  ZONE T= \"HP\", N=  " << this->All_Points.size();
	fout << " , E= " << ll;
	fout << " , F=FEPOINT, ET=LINESEG  " << endl;
	for (auto& i : this->All_Points)
	{
		fout << i->x << " " << i->y << endl;
	}

	for (auto& i : this->All_Cells)
	{
		for (int j = 0; j < i->contour.size() - 1; j++)
		{
			fout << i->contour[j]->number + 1 << " " << i->contour[j + 1]->number + 1 << endl;
		}
		fout << i->contour[i->contour.size() - 1]->number + 1 << " " << i->contour[0]->number + 1 << endl;
	}

	fout.close();
}

void Setka::Print_connect(void)
{
	int ll = 0;
	for (auto& i : this->All_Cells)
	{
		ll += i->Grans.size();
	}
	double x, y;
	double x2, y2;
	ofstream fout;
	fout.open("Setka_connect.txt");
	fout << "TITLE = \"HP\" ";
	fout << " VARIABLES = \"X\", \"Y\"  ZONE T= \"HP\", N=  " << ll * 2;
	fout << " , E= " << ll;
	fout << " , F=FEPOINT, ET=LINESEG  " << endl;
	for (auto& i : this->All_Cells)
	{
		for (auto& j : i->Grans)
		{
			if (false) // (j->type == Usualy)
			{
				i->Get_Center(x, y);
				j->Sosed->Get_Center(x2, y2);
				fout << x << " " << y << endl;
				fout << (x + x2) / 2.0 << " " << (y + y2) / 2.0 << endl;
			}
			else
			{
				i->Get_Center(x, y);
				j->Get_Center(x2, y2);
				fout << x << " " << y << endl;
				fout << x2 << " " << y2 << endl;
			}
		}
	}

	for (int i = 0; i < ll; i++)
	{
		fout << 2 * i + 1 << " " << 2 * i + 2 << endl;
	}

	fout.close();
}

void Setka::Print_point_connect(void)
{
	int ll = 0;
	for (auto& i : this->All_Points)
	{
		ll += i->my_cell.size();
	}
	double x, y;
	double x2, y2;
	ofstream fout;
	fout.open("Setka_connect_point.txt");
	fout << "TITLE = \"HP\" ";
	fout << " VARIABLES = \"X\", \"Y\"  ZONE T= \"HP\", N=  " << ll * 2;
	fout << " , E= " << ll;
	fout << " , F=FEPOINT, ET=LINESEG  " << endl;
	for (auto& i : this->All_Points)
	{
		x = i->x;
		y = i->y;
		for (auto& j : i->my_cell)
		{
			j->Get_Center(x2, y2);
			fout << x << " " << y << endl;
			fout << x2 << " " << y2 << endl;
		}
	}

	for (int i = 0; i < ll; i++)
	{
		fout << 2 * i + 1 << " " << 2 * i + 2 << endl;
	}

	fout.close();
}

void Setka::Print_Tecplot(void)
{
	double r_o = 1.0;    // Размер расстояния
	double ro_o = 0.06;   // Размер плотности
	double ro_o_H = 0.18;   // Размер плотности
	double p_o = 1.0;   // Размер давления
	double u_o = 10.38;   // Размер давления

	ofstream fout;
	string name_f = "2D_tecplot.txt";
	fout.open(name_f);
	fout << "TITLE = \"HP\"  VARIABLES = \"X\", \"Y\", \"r\", \"Ro\", \"P\", \"Vx\", \"Vy\", \"Max\",\"Q\",\"Ro_H1\", \"P_H1\", \"Vx_H1\", \"Vy_H1\"," << //
		"\"Ro_H2\", \"P_H2\", \"Vx_H2\", \"Vy_H2\",\"Ro_H3\", \"P_H3\", \"Vx_H3\", \"Vy_H3\",\"Ro_H4\", \"P_H4\", \"Vx_H4\", \"Vy_H4\", \"RO_H\", ZONE T = \"HP\"" << endl;
	for (auto& i : this->All_Cells)
	{
		double kk;
		if (i->type == C_centr || i->type == C_1 || i->type == C_2 || i->type == C_3)
		{
			kk = (chi_real / chi_);
		}
		else
		{
			kk = 1.0;
		}
		//double kk = 1.0;
		double Max = 0.0;
		double QQ = 0.0;
		double x, y;
		i->Get_Center(x, y);
		if (i->par[0].ro > 0.000000000001)
		{
			QQ = i->par[0].Q / i->par[0].ro;
			Max = sqrt(kvv(i->par[0].u, i->par[0].v, 0.0) / (ggg * i->par[0].p / i->par[0].ro));
		}

		fout << x * r_o << " " << y * r_o << " " << sqrt(x * r_o * x * r_o + y * r_o * y * r_o) << //
			" " << i->par[0].ro * ro_o/kv(kk) << " " << i->par[0].p * p_o << " " //
			<< i->par[0].u * u_o * kk << " " << i->par[0].v * u_o * kk << " " << Max << " " << QQ << //
			" " << i->par[0].ro_H1 * ro_o_H << " " << i->par[0].p_H1 * p_o << " " //
			<< i->par[0].u_H1 * u_o << " " << i->par[0].v_H1 * u_o << //
			" " << i->par[0].ro_H2 * ro_o_H << " " << i->par[0].p_H2 * p_o << " " //
			<< i->par[0].u_H2 * u_o << " " << i->par[0].v_H2 * u_o << //
			" " << i->par[0].ro_H3 * ro_o_H << " " << i->par[0].p_H3 * p_o << " " //
			<< i->par[0].u_H3 * u_o << " " << i->par[0].v_H3 * u_o << //
			" " << i->par[0].ro_H4 * ro_o_H << " " << i->par[0].p_H4 * p_o << " " //
			<< i->par[0].u_H4 * u_o << " " << i->par[0].v_H4 * u_o << " " <<//
			(i->par[0].ro_H1 + i->par[0].ro_H2 + i->par[0].ro_H3 + i->par[0].ro_H4) * ro_o_H << endl;

	}
	fout.close();

	name_f = "1D_tecplot.txt";
	fout.open(name_f);
	fout << "TITLE = \"HP\"  VARIABLES = \"X\", \"Y\", \"r\", \"Ro\", \"Ro_r_r\", \"P\", \"Vx\", \"Vy\", \"Max\",\"Q\",\"Ro_H1\", \"P_H1\", \"Vx_H1\", \"Vy_H1\"," << //
		"\"M_H1\", \"T_H1\", \"Ro_H2\", \"P_H2\", \"Vx_H2\", \"Vy_H2\",\"Ro_H3\", \"P_H3\", \"Vx_H3\", \"Vy_H3\",\"Ro_H4\", \"P_H4\", \"Vx_H4\", \"Vy_H4\", \"RO_H\", ZONE T = \"HP\"" << endl;
	int num = 0;
	//double ro = (389.988 * 389.988) / (chi_ * chi_);


	fout << 1 * r_o << " " << 0.0 * r_o << " " << 1.0 << //
		" " << 0.06 * (389.988 * 389.988) / (chi_real * chi_real) << " " << 0.0 << " " << 0.0 << " " //
		<< 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0  << //
		" " << 0.001666666 * ro_o_H << " " << (0.001666666 * chi_real * chi_real / (ggg * 50.0 * 50.0)) * pow(1.0 / 1.0, 2.0 * ggg) << " " //
		<< chi_real * u_o << " " << 0.0 << //
		" " << 50.0 << " " << 0.0 << " " << 0.0 << " " //
		<< 0.0 << " " << 0.0 << //
		" " << 0.0 << " " << 0.0 << " " //
		<< 0.0 << " " << 0.0 << //
		" " << 0.0 << " " << 0.0 << " " << 0.0 << " " //
		<< 0.0 << " " << 0.0 << " " <<//
		0.0 << endl;


	for (auto& i : this->All_Cells)
	{
		double kk;
		if (i->type == C_centr || i->type == C_1 || i->type == C_2 || i->type == C_3)
		{
			kk = (chi_real / chi_);
		}
		else
		{
			kk = 1.0;
		}
		//kk = 1.0;
		num++;
		if (num > this->M1 + this->M2 + this->M3 + this->M4)
		{
			break;
		}
		double Max = 0.0;
		double Max_H1 = 0.0;
		double T_H1 = 0.0;
		double QQ = 0.0;
		double x, y;
		i->Get_Center(x, y);
		if (i->par[0].ro > 0.000000000001)
		{
			QQ = i->par[0].Q / i->par[0].ro;
			Max = sqrt(kvv(i->par[0].u, i->par[0].v, 0.0) / (ggg * i->par[0].p / i->par[0].ro));
		}
		if (i->par[0].ro_H1 > 0.000000000001)
		{
			Max_H1 = sqrt(kvv(i->par[0].u_H1, i->par[0].v_H1, 0.0) / (ggg * i->par[0].p_H1 / i->par[0].ro_H1));
			T_H1 = 2.0 * (i->par[0].p_H1 / i->par[0].ro_H1) * 6527.0;
		}
		double dist = sqrt(kv(x) + kv(y));
		fout << x * r_o << " " << y * r_o << " " << sqrt(x * r_o * x * r_o + y * r_o * y * r_o) << //
			" " << i->par[0].ro * ro_o / kv(kk) << " " << 0.06 * (389.988 * 389.988) / (chi_real * chi_real)/ (dist * dist) << " " << i->par[0].p * p_o << " " //
			<< i->par[0].u * u_o * kk << " " << i->par[0].v * u_o * kk << " " << Max << " " << QQ << //
			" " << i->par[0].ro_H1 * ro_o_H << " " << i->par[0].p_H1 * p_o << " "//
			<< i->par[0].u_H1 * u_o << " " << i->par[0].v_H1 * u_o << " " << Max_H1 << " " << T_H1 << " "//
			<< i->par[0].ro_H2 * ro_o_H << " " << i->par[0].p_H2 * p_o << " " //
			<< i->par[0].u_H2 * u_o << " " << i->par[0].v_H2 * u_o << //
			" " << i->par[0].ro_H3 * ro_o_H << " " << i->par[0].p_H3 * p_o << " " //
			<< i->par[0].u_H3 * u_o << " " << i->par[0].v_H3 * u_o << //
			" " << i->par[0].ro_H4 * ro_o_H << " " << i->par[0].p_H4 * p_o << " " //
			<< i->par[0].u_H4 * u_o << " " << i->par[0].v_H4 * u_o << " " <<//
			(i->par[0].ro_H1 + i->par[0].ro_H2 + i->par[0].ro_H3 + i->par[0].ro_H4) * ro_o_H << endl;

	}
	fout.close();
}

void Setka::Print_Sourse(void)
{
	string name_f = "1D_Sourse.txt";
	ofstream fout;
	fout.open(name_f);
	fout << "TITLE = \"HP\"  VARIABLES = \"X\", \"I1_mc_1\",  \"I1_mc_2\", \"I1_mc_3\", \"I1_mc_4\", \"I2_mc_1\",  \"I2_mc_2\", \"I2_mc_3\", \"I2_mc_4\", \"I3_mc_1\",  \"I3_mc_2\", \"I3_mc_3\", \"I3_mc_4\", \"I1_mf_1\",  \"I1_mf_2\", \"I1_mf_3\", \"I1_mf_4\", \"I2_mf_1\",  \"I2_mf_2\", \"I2_mf_3\", \"I2_mf_4\", \"I3_mf_1\",  \"I3_mf_2\", \"I3_mf_3\", \"I3_mf_4\",  " << //
		"ZONE T = \"HP\"" << endl;
	int num = 0;
	for (auto& i : this->All_Cells)
	{
		num++;
		if (num > this->M1 + this->M2 + this->M3 + this->M4)
		{
			break;
		}
		
		double x, y;
		i->Get_Center(x, y);
		i->Calc_Sourse();
		double dist = sqrt(kv(x) + kv(y));
		fout << sqrt(kvv(0.0, x, y)) << //
			" " << i->par[0].I1_mc[0] << " " << i->par[0].I1_mc[1] << " " << i->par[0].I1_mc[2] << " " << i->par[0].I1_mc[3] << //
			" " << i->par[0].I2_mc[0] << " " << i->par[0].I2_mc[1] << " " << i->par[0].I2_mc[2] << " " << i->par[0].I2_mc[3] << //
			" " << i->par[0].I3_mc[0] << " " << i->par[0].I3_mc[1] << " " << i->par[0].I3_mc[2] << " " << i->par[0].I3_mc[3] << //
			" " << i->par[0].I1_mf[0] << " " << i->par[0].I1_mf[1] << " " << i->par[0].I1_mf[2] << " " << i->par[0].I1_mf[3] << //
			" " << i->par[0].I2_mf[0] << " " << i->par[0].I2_mf[1] << " " << i->par[0].I2_mf[2] << " " << i->par[0].I2_mf[3] << //
			" " << i->par[0].I3_mf[0] << " " << i->par[0].I3_mf[1] << " " << i->par[0].I3_mf[2] << " " << i->par[0].I3_mf[3] << endl;

	}
	fout.close();
}

void Setka::Print_Tecplot_MK_2d(string name0)
{
	double r_o = 1.0; // RR_;    // Размер расстояния
	double ro_o = 1.0; // 0.06;   // Размер плотности
	double ro_o_H = 1.0;    //0.18;   // Размер плотности
	double p_o = 1.0;    //1.0;   // Размер давления
	double u_o = 1.0; // 10.38;   // Размер скорости
	double T_o = 1.0;  // 6527.0
	int n;

	n = this->All_Cells.size();

	Cell* A, *B, *C, *D, *E, *M;
	E = nullptr;

	ofstream fout;
	string name_f = "Interpol_2D_tecplot_" + name0 + ".txt";
	fout.open(name_f);
	fout << "TITLE = \"HP\"  VARIABLES = \"X\", \"Y\", \"r\", \"Rho\", \"P\", \"Vx\", \"Vy\", \"M\",\"Q\",\"Rho_H1\", \"P_H1\", \"Vx_H1\", \"Vy_H1\"," << //
		"\"Rho_H2\", \"P_H2\", \"Vx_H2\", \"Vy_H2\",\"Rho_H3\", \"P_H3\", \"Vx_H3\", \"Vy_H3\",\"Rho_H4\", \"P_H4\"," << //
		" \"Vx_H4\", \"Vy_H4\", \"RO_H\", \"F_n\", \"F_u\", \"F_v\", \"F_T\", \"I_u\", \"I_v\", \"I_T\", \"II_u\", \"II_v\", \"II_T\", \"M_u\", \"M_v\", \"M_T\",\"H1_n\", \"H1_u\", \"H1_v\", \"H1_T\"," << //
		"\"H2_n\", \"H2_u\", \"H2_v\", \"H2_T\", \"H3_n\", \"H3_u\", \"H3_v\", \"H3_T\", \"H4_n\", \"H4_u\", \"H4_v\", \"H4_T\", \"H5_n\", \"H5_u\", \"H5_v\", \"H5_T\", \"H6_n\", \"H6_u\", \"H6_v\", \"H6_T\", \"k_u\", \"k_v\", \"k_T\", \"m_1\",  \"m_2\",  \"m_3\",  \"m_4\",  \"m_5\",  \"m_6\",  \"m_7\",  \"num1\", \"num2\", \"num3\", \"num4\", ZONE T = \"HP\"" << //
		", N= " << 5 * (n) + (Line_Inner.size() - 1) * 4 << ", E =  "<<  4 * (n) + 2 * (Line_Inner.size() - 1)<< ", F=FEPOINT, ET=TRIANGLE "<< endl;
	
	for (auto i : this->All_Cells)
	{
		double kk = 1.0;
		double Max = 0.0;
		double QQ = 0.0;
		double x, y;
		double x1, y1;
		double x2, y2;
		double x3, y3;
		double x4, y4;
		double x5, y5;
		double u1, v1;
		double u2, v2;
		double u3, v3;
		double u4, v4;
		double n1, n2, n3, n4;
		bool b1, b2, b3;

		b1 = b2 = b3 = false;

		A = i;
		B = i->Grans[0]->Sosed;
		C = i->Grans[1]->Sosed;
		D = i->Grans[2]->Sosed;

		if (i->Grans.size() == 3)
		{
			E = A;
		}
		else if (i->Grans.size() <= 3)
		{
			cout << "ERROR  56278476526r7duhf34r" << endl;
		}
		else
		{
			E = i->Grans[3]->Sosed;
		}

		if (B == nullptr)  B = A;
		if (C == nullptr)  C = A;
		if (D == nullptr)  D = A;
		if (E == nullptr)  E = A;

		a1:
		A->Get_Center(x1, y1);
		B->Get_Center(x2, y2);
		C->Get_Center(x3, y3);
		D->Get_Center(x4, y4);
		E->Get_Center(x5, y5);

		u1 = x2 - x1;
		v1 = y2 - y1;

		u2 = x3 - x1;
		v2 = y3 - y1;

		u3 = x4 - x1;
		v3 = y4 - y1;

		u4 = x5 - x1;
		v4 = y5 - y1;

		n1 = sqrt(u1 * u1 + v1 * v1);
		n2 = sqrt(u2 * u2 + v2 * v2);
		n3 = sqrt(u3 * u3 + v3 * v3);
		n4 = sqrt(u4 * u4 + v4 * v4);

		if (n1 < 0.000001 || n2 < 0.000001 || n3 < 0.000001 || n4 < 0.000001) goto a2;

		u1 /= n1;
		v1 /= n1;

		u2 /= n2;
		v2 /= n2;

		u3 /= n3;
		v3 /= n3;

		u4 /= n4;
		v4 /= n4;

		if (u1 * u2 + v1 * v2 < -0.95 && b1 == false)
		{
			b1 = true;
			M = B;
			B = C;
			C = D;
			D = E;
			E = M;
			goto a1;
		}
		else if(u3 * u2 + v3 * v2 < -0.95 && b2 == false)
		{
			b2 = true;
			M = C;
			C = D;
			D = E;
			E = M;
			goto a1;
		}
		else if (u3 * u4 + v3 * v4 < -0.95 && b3 == false)
		{
			b3 = true;
			M = E;
			D = E;
			E = M;
			goto a1;
		}
		a2:

		for (int ij = 1; ij <= 5; ij++)
		{
			if (ij == 1) i = A;
			if (ij == 2) i = B;
			if (ij == 3) i = C;
			if (ij == 4) i = D;
			if (ij == 5) i = E;

			i->Get_Center(x, y);
			if (i->par[0].ro > 0.000000000001)
			{
				QQ = i->par[0].Q / i->par[0].ro;
				Max = sqrt(kvv(i->par[0].u, i->par[0].v, 0.0) / (ggg * i->par[0].p / i->par[0].ro));
			}

			fout << x << " " << y << " " << sqrt(x * r_o * x * r_o + y * r_o * y * r_o) << //
				" " << i->par[0].ro * ro_o / kv(kk) << " " << i->par[0].p * p_o << " " //
				<< i->par[0].u * u_o * kk << " " << i->par[0].v * u_o * kk << " " << Max << " " << QQ << //
				" " << i->par[0].ro_H1 * ro_o_H << " " << i->par[0].p_H1 * p_o << " " //
				<< i->par[0].u_H1 * u_o << " " << i->par[0].v_H1 * u_o << //
				" " << i->par[0].ro_H2 * ro_o_H << " " << i->par[0].p_H2 * p_o << " " //
				<< i->par[0].u_H2 * u_o << " " << i->par[0].v_H2 * u_o << //
				" " << i->par[0].ro_H3 * ro_o_H << " " << i->par[0].p_H3 * p_o << " " //
				<< i->par[0].u_H3 * u_o << " " << i->par[0].v_H3 * u_o << //
				" " << i->par[0].ro_H4 * ro_o_H << " " << i->par[0].p_H4 * p_o << " " //
				<< i->par[0].u_H4 * u_o << " " << i->par[0].v_H4 * u_o << " " <<//
				(i->par[0].ro_H1 + i->par[0].ro_H2 + i->par[0].ro_H3 + i->par[0].ro_H4) * ro_o_H << //
				" " << i->par[0].F_n << " " << i->par[0].F_u << " " << i->par[0].F_v << " " << i->par[0].F_T << " " //
				<< i->par[0].I_u * 0.00001107 << " " << i->par[0].I_v * 0.00001107 << " " << i->par[0].I_T * 11.4767 << " " << //
				i->par[0].II_u * 0.00001107 << " " << i->par[0].II_v * 0.00001107 << " " << i->par[0].II_T * 11.4767 << " " << //
				i->par[0].M_u * 0.00001107 << " " << i->par[0].M_v * 0.00001107 << " " << i->par[0].M_T * 11.4767 << " " << //
				i->par[0].H_n[0] * ro_o_H << " " << i->par[0].H_u[0] * u_o << " " << i->par[0].H_v[0] * u_o << " " << i->par[0].H_T[0] * T_o << " " //
				<< i->par[0].H_n[1] * ro_o_H << " " << i->par[0].H_u[1] * u_o << " " << i->par[0].H_v[1] * u_o << " " << i->par[0].H_T[1] * T_o << " " //
				<< i->par[0].H_n[2] * ro_o_H << " " << i->par[0].H_u[2] * u_o << " " << i->par[0].H_v[2] * u_o << " " << i->par[0].H_T[2] * T_o << " " //
				<< i->par[0].H_n[3] * ro_o_H << " " << i->par[0].H_u[3] * u_o << " " << i->par[0].H_v[3] * u_o << " " << i->par[0].H_T[3] * T_o << " " //
				<< i->par[0].H_n[4] * ro_o_H << " " << i->par[0].H_u[4] * u_o << " " << i->par[0].H_v[4] * u_o << " " << i->par[0].H_T[4] * T_o << " " //
				<< i->par[0].H_n[5] * ro_o_H << " " << i->par[0].H_u[5] * u_o << " " << i->par[0].H_v[5] * u_o << " " << i->par[0].H_T[5] * T_o << //
				" " << i->par[0].k_u << " " << i->par[0].k_v << " " << i->par[0].k_T << //
				" " << i->par[0].w_m[0] << //
				" " << i->par[0].w_m[1] << //
				" " << i->par[0].w_m[2] << //
				" " << i->par[0].w_m[3] << //
				" " << i->par[0].w_m[4] << //
				" " << i->par[0].w_m[5] << //
				" " << i->par[0].w_m[6] << //
				" " << i->par[0].num_atoms[0] << //
				" " << i->par[0].num_atoms[1] << //
				" " << i->par[0].num_atoms[2] << //
				" " << i->par[0].num_atoms[3] << endl;
		}

	}

	for (int ii = 0; ii < Line_Inner.size() - 1; ii++)
	{
		double kk = 1.0;
		double Max = 0.0;
		double QQ = 0.0;
		double x, y;
		auto m1 = Line_Inner[ii];
		auto m2 = Line_Inner[ii + 1];
		Cell* i;
		i = nullptr;

		A = m1->Master;
		B = m2->Master;

		for (int ij = 1; ij <= 4; ij++)
		{
			if (ij == 1) i = A;
			if (ij == 2) i = A;
			if (ij == 3) i = B;
			if (ij == 4) i = B;

			if (ij == 1) i->Get_Center(x, y);
			if (ij == 2) m1->Get_Center(x, y);
			if (ij == 3) i->Get_Center(x, y);
			if (ij == 4) m2->Get_Center(x, y);


			if (i->par[0].ro > 0.000000000001)
			{
				QQ = i->par[0].Q / i->par[0].ro;
				Max = sqrt(kvv(i->par[0].u, i->par[0].v, 0.0) / (ggg * i->par[0].p / i->par[0].ro));
			}

			fout << x << " " << y << " " << sqrt(x * r_o * x * r_o + y * r_o * y * r_o) << //
				" " << i->par[0].ro * ro_o / kv(kk) << " " << i->par[0].p * p_o << " " //
				<< i->par[0].u * u_o * kk << " " << i->par[0].v * u_o * kk << " " << Max << " " << QQ << //
				" " << i->par[0].ro_H1 * ro_o_H << " " << i->par[0].p_H1 * p_o << " " //
				<< i->par[0].u_H1 * u_o << " " << i->par[0].v_H1 * u_o << //
				" " << i->par[0].ro_H2 * ro_o_H << " " << i->par[0].p_H2 * p_o << " " //
				<< i->par[0].u_H2 * u_o << " " << i->par[0].v_H2 * u_o << //
				" " << i->par[0].ro_H3 * ro_o_H << " " << i->par[0].p_H3 * p_o << " " //
				<< i->par[0].u_H3 * u_o << " " << i->par[0].v_H3 * u_o << //
				" " << i->par[0].ro_H4 * ro_o_H << " " << i->par[0].p_H4 * p_o << " " //
				<< i->par[0].u_H4 * u_o << " " << i->par[0].v_H4 * u_o << " " <<//
				(i->par[0].ro_H1 + i->par[0].ro_H2 + i->par[0].ro_H3 + i->par[0].ro_H4) * ro_o_H << //
				" " << i->par[0].F_n << " " << i->par[0].F_u << " " << i->par[0].F_v << " " << i->par[0].F_T << " " //
				<< i->par[0].I_u * 0.00001107 << " " << i->par[0].I_v * 0.00001107 << " " << i->par[0].I_T * 11.4767 << " " << //
				i->par[0].II_u * 0.00001107 << " " << i->par[0].II_v * 0.00001107 << " " << i->par[0].II_T * 11.4767 << " " << //
				i->par[0].M_u * 0.00001107 << " " << i->par[0].M_v * 0.00001107 << " " << i->par[0].M_T * 11.4767 << " " << //
				i->par[0].H_n[0] * ro_o_H << " " << i->par[0].H_u[0] * u_o << " " << i->par[0].H_v[0] * u_o << " " << i->par[0].H_T[0] * T_o << " " //
				<< i->par[0].H_n[1] * ro_o_H << " " << i->par[0].H_u[1] * u_o << " " << i->par[0].H_v[1] * u_o << " " << i->par[0].H_T[1] * T_o << " " //
				<< i->par[0].H_n[2] * ro_o_H << " " << i->par[0].H_u[2] * u_o << " " << i->par[0].H_v[2] * u_o << " " << i->par[0].H_T[2] * T_o << " " //
				<< i->par[0].H_n[3] * ro_o_H << " " << i->par[0].H_u[3] * u_o << " " << i->par[0].H_v[3] * u_o << " " << i->par[0].H_T[3] * T_o << " " //
				<< i->par[0].H_n[4] * ro_o_H << " " << i->par[0].H_u[4] * u_o << " " << i->par[0].H_v[4] * u_o << " " << i->par[0].H_T[4] * T_o << " " //
				<< i->par[0].H_n[5] * ro_o_H << " " << i->par[0].H_u[5] * u_o << " " << i->par[0].H_v[5] * u_o << " " << i->par[0].H_T[5] * T_o << //
				" " << i->par[0].k_u << " " << i->par[0].k_v << " " << i->par[0].k_T << //
				" " << i->par[0].w_m[0] << //
				" " << i->par[0].w_m[1] << //
				" " << i->par[0].w_m[2] << //
				" " << i->par[0].w_m[3] << //
				" " << i->par[0].w_m[4] << //
				" " << i->par[0].w_m[5] << //
				" " << i->par[0].w_m[6] << //
				" " << i->par[0].num_atoms[0] << //
				" " << i->par[0].num_atoms[1] << //
				" " << i->par[0].num_atoms[2] << //
				" " << i->par[0].num_atoms[3] << endl;
		}
	}

	for (int ij = 0; ij < n; ij++)
	{
		fout << 1 + 5 * ij << " " << 2 + 5 * ij << " " << 3 + 5 * ij << endl;
		fout << 1 + 5 * ij << " " << 3 + 5 * ij << " " << 4 + 5 * ij << endl;
		fout << 1 + 5 * ij << " " << 4 + 5 * ij << " " << 5 + 5 * ij << endl;
		fout << 1 + 5 * ij << " " << 5 + 5 * ij << " " << 2 + 5 * ij << endl;
	}


	for (int ij = 0; ij < Line_Inner.size() - 1; ij++)
	{
		fout << 5*n + 1 + 4 * ij << " " << 5 * n + 2 + 4 * ij << " " << 5 * n + 3 + 4 * ij << endl;
		fout << 5 * n + 2 + 4 * ij << " " << 5 * n + 3 + 4 * ij << " " << 5 * n + 4 + 4 * ij << endl;
	}

	fout.close();

}


void Setka::Print_Tecplot_MK(string name0)
{
	double r_o = 1.0; // RR_;    // Размер расстояния
	double ro_o = 1.0; // 0.06;   // Размер плотности
	double ro_o_H = 1.0;    //0.18;   // Размер плотности
	double p_o = 1.0;    //1.0;   // Размер давления
	double u_o = 1.0; // 10.38;   // Размер скорости
	double T_o = 1.0;  // 6527.0

	ofstream fout;
	string name_f = "2D_tecplot_" + name0 + ".txt";
	fout.open(name_f);
	fout << "TITLE = \"HP\"  VARIABLES = \"X\", \"Y\", \"r\", \"Rho\", \"P\", \"Vx\", \"Vy\", \"M\",\"|V|\",\"Rho_H1\", \"P_H1\", \"Vx_H1\", \"Vy_H1\"," << //
		"\"Rho_H2\", \"P_H2\", \"Vx_H2\", \"Vy_H2\",\"Rho_H3\", \"P_H3\", \"Vx_H3\", \"Vy_H3\",\"Rho_H4\", \"P_H4\"," << //
		" \"Vx_H4\", \"Vy_H4\", \"RO_H\", \"F_n\", \"F_u\", \"F_v\", \"F_T\", \"I_u\", \"I_v\", \"I_T\", \"II_u\", \"II_v\", \"II_T\", \"M_u\", \"M_v\", \"M_T\",\"H1_n\", \"H1_u\", \"H1_v\", \"H1_T\"," << //
		"\"H2_n\", \"H2_u\", \"H2_v\", \"H2_T\", \"H3_n\", \"H3_u\", \"H3_v\", \"H3_T\", \"H4_n\", \"H4_u\", \"H4_v\", \"H4_T\", \"H5_n\", \"H5_u\", \"H5_v\", \"H5_T\", \"H6_n\", \"H6_u\", \"H6_v\", \"H6_T\", \"k_u\", \"k_v\", \"k_T\", \"m_1\",  \"m_2\",  \"m_3\",  \"m_4\",  \"m_5\",  \"m_6\",  \"m_7\",  \"num1\", \"num2\", \"num3\", \"num4\", ZONE T = \"HP\"" << endl;
	for (auto& i : this->All_Cells)
	{
		double kk = 1.0;
		/*if (i->type == C_centr || i->type == C_1 || i->type == C_2 || i->type == C_3)
		{
			kk = (chi_real / chi_);
		}
		else
		{
			kk = 1.0;
		}*/
		double Max = 0.0;
		double QQ = 0.0;
		double x, y;
		i->Get_Center(x, y);
		if (i->par[0].ro > 0.000000000001)
		{
			Max = sqrt(kvv(i->par[0].u, i->par[0].v, 0.0) / (ggg * i->par[0].p / i->par[0].ro));
		}

		QQ = sqrt(kvv(i->par[0].u, i->par[0].v, 0.0));

		fout << x << " " << y << " " << sqrt(x * r_o * x * r_o + y * r_o * y * r_o) << //
			" " << i->par[0].ro * ro_o / kv(kk) << " " << i->par[0].p * p_o << " " //
			<< i->par[0].u * u_o * kk << " " << i->par[0].v * u_o * kk << " " << Max << " " << QQ << //
			" " << i->par[0].ro_H1 * ro_o_H << " " << i->par[0].p_H1 * p_o << " " //
			<< i->par[0].u_H1 * u_o << " " << i->par[0].v_H1 * u_o << //
			" " << i->par[0].ro_H2 * ro_o_H << " " << i->par[0].p_H2 * p_o << " " //
			<< i->par[0].u_H2 * u_o << " " << i->par[0].v_H2 * u_o << //
			" " << i->par[0].ro_H3 * ro_o_H << " " << i->par[0].p_H3 * p_o << " " //
			<< i->par[0].u_H3 * u_o << " " << i->par[0].v_H3 * u_o << //
			" " << i->par[0].ro_H4 * ro_o_H << " " << i->par[0].p_H4 * p_o << " " //
			<< i->par[0].u_H4 * u_o << " " << i->par[0].v_H4 * u_o << " " <<//
			(i->par[0].ro_H1 + i->par[0].ro_H2 + i->par[0].ro_H3 + i->par[0].ro_H4) * ro_o_H << //
			" " << i->par[0].F_n << " " << i->par[0].F_u << " " << i->par[0].F_v << " " << i->par[0].F_T << " " //
			<< i->par[0].I_u * 0.00001107 << " " << i->par[0].I_v * 0.00001107 << " " << i->par[0].I_T * 11.4767 << " " << //
			i->par[0].II_u * 0.00001107 << " " << i->par[0].II_v * 0.00001107 << " " << i->par[0].II_T * 11.4767 << " " << //
			i->par[0].M_u * 0.00001107 << " " << i->par[0].M_v * 0.00001107 << " " << i->par[0].M_T * 11.4767 << " " << //
			i->par[0].H_n[0] * ro_o_H << " " << i->par[0].H_u[0] * u_o << " " << i->par[0].H_v[0] * u_o << " " << i->par[0].H_T[0] * T_o << " " //
			<< i->par[0].H_n[1] * ro_o_H << " " << i->par[0].H_u[1] * u_o << " " << i->par[0].H_v[1] * u_o << " " << i->par[0].H_T[1] * T_o << " " //
			<< i->par[0].H_n[2] * ro_o_H << " " << i->par[0].H_u[2] * u_o << " " << i->par[0].H_v[2] * u_o << " " << i->par[0].H_T[2] * T_o << " " //
			<< i->par[0].H_n[3] * ro_o_H << " " << i->par[0].H_u[3] * u_o << " " << i->par[0].H_v[3] * u_o << " " << i->par[0].H_T[3] * T_o << " " //
			<< i->par[0].H_n[4] * ro_o_H << " " << i->par[0].H_u[4] * u_o << " " << i->par[0].H_v[4] * u_o << " " << i->par[0].H_T[4] * T_o << " " //
			<< i->par[0].H_n[5] * ro_o_H << " " << i->par[0].H_u[5] * u_o << " " << i->par[0].H_v[5] * u_o << " " << i->par[0].H_T[5] * T_o << //
			" " << i->par[0].k_u << " " << i->par[0].k_v << " " << i->par[0].k_T << //
			" " << i->par[0].w_m[0] << //
			" " << i->par[0].w_m[1] << //
			" " << i->par[0].w_m[2] << //
			" " << i->par[0].w_m[3] << //
			" " << i->par[0].w_m[4] << //
			" " << i->par[0].w_m[5] << //
			" " << i->par[0].w_m[6] << //
			" " << i->par[0].num_atoms[0] << //
			" " << i->par[0].num_atoms[1] << //
			" " << i->par[0].num_atoms[2] << //
			" " << i->par[0].num_atoms[3] << endl;

	}
	fout.close();

	name_f = "1D_tecplot_new_" + name0 + ".txt";
	fout.open(name_f);
	fout << "TITLE = \"HP\"  VARIABLES = \"X\", \"Y\", \"r\", \"Ro\", \"Ro_r_r\", \"P\", \"Vx\", \"Vy\", \"T\", \"Max\",\"Q\",\"Ro_H1\", \"T_H1\", \"Vx_H1\", \"Vy_H1\"," << //
		" \"Ro_H2\", \"T_H2\", \"Vx_H2\", \"Vy_H2\",\"Ro_H3\", \"T_H3\", \"Vx_H3\", \"Vy_H3\",\"Ro_H4\", \"T_H4\", \"Vx_H4\", \"Vy_H4\", \"RO_H\", \"T_H\", \"I_u\", \"I_v\", \"I_T\"," << //
		"\"zone\", ZONE T = \"HP\"" << endl;
	int num = 0;
	//double ro = (389.988 * 389.988) / (chi_ * chi_);


	for (auto& i : this->All_Cells)
	{
		num++;
		if (num > this->M1 + this->M2 + this->M3 + this->M4)
		{
			break;
		}
		double Max = 0.0;
		double Max_H1 = 0.0;
		double T_H1 = 0.0;
		double QQ = 0.0;
		double x, y;
		i->Get_Center(x, y);
		if (i->par[0].ro > 0.000000000001)
		{
			QQ = i->par[0].Q / i->par[0].ro;
			Max = sqrt(kvv(i->par[0].u, i->par[0].v, 0.0) / (ggg * i->par[0].p / i->par[0].ro));
		}
		if (i->par[0].ro_H1 > 0.000000000001)
		{
			// Max_H1 = sqrt(kvv(i->par[0].u_H1, i->par[0].v_H1, 0.0) / (ggg * i->par[0].p_H1 / i->par[0].ro_H1));
			// T_H1 = 2.0 * (i->par[0].p_H1 / i->par[0].ro_H1) * 6527.0;
		}
		double dist = sqrt(kv(x) + kv(y));
		fout << x * r_o << " " << y * r_o << " " << sqrt(x * r_o * x * r_o + y * r_o * y * r_o) << //
			" " << i->par[0].ro * ro_o  << " " << 0.06 * (389.988 * 389.988) / (chi_real * chi_real) / (dist * dist) << " " << i->par[0].p * p_o << " " //
			<< i->par[0].u * u_o  << " " << i->par[0].v * u_o << " " << i->par[0].p/ i->par[0].ro << " " << Max << " " << QQ << //
			" " << i->par[0].H_n[0] * ro_o_H << " " << i->par[0].H_T[0] * T_o << " "//
			<< i->par[0].H_u[0] * u_o << " " << i->par[0].H_v[0] * u_o << " " //
			<< i->par[0].H_n[1] * ro_o_H << " " << i->par[0].H_T[1] * T_o << " " //
			<< i->par[0].H_u[1] * u_o << " " << i->par[0].H_v[1] * u_o << //
			" " << i->par[0].H_n[2] * ro_o_H << " " << i->par[0].H_T[2] * T_o << " " //
			<< i->par[0].H_u[2] * u_o << " " << i->par[0].H_v[2] * u_o << //
			" " << i->par[0].H_n[3] * ro_o_H << " " << i->par[0].H_T[3] * T_o << " " //
			<< i->par[0].H_u[3] * u_o << " " << i->par[0].H_v[3] * u_o << " " <<//
			(i->par[0].H_n[0] + i->par[0].H_n[1] + i->par[0].H_n[2] + i->par[0].H_n[3]) * ro_o_H << //
			" " << (i->par[0].H_n[0] * i->par[0].H_T[0] + i->par[0].H_n[1] * i->par[0].H_T[1] + i->par[0].H_n[2] * i->par[0].H_T[2] + i->par[0].H_n[3] * i->par[0].H_T[3])/ (i->par[0].H_n[0] + i->par[0].H_n[1] + i->par[0].H_n[2] + i->par[0].H_n[3]) << //
			" " << i->par[0].I_u<< " " << i->par[0].I_v << " " << i->par[0].I_T << " " << i->type << endl;

	}
	fout.close();


	name_f = "1D_tecplot_discontinity_" + name0 + ".txt";
	fout.open(name_f);
	fout << "TITLE = \"HP\"  VARIABLES = \"X\", \"Ro\", " << //
		"\"zone\", ZONE T = \"HP\"" << endl;

	for(int nn = 0; nn < this->M1 + this->M2 + this->M3 + this->M4; nn++)
	{
		auto i = this->All_Cells[nn];
		

		double x, y;
		i->Get_Center(x, y);

		fout << x << " " <<  i->par[0].ro << " " << i->type << endl;

		if (nn == this->M1 + this->M2 + this->M3 - 1)
		{
			double x2, y2;
			this->All_Cells[nn + 1]->Get_Center(x2, y2);

			double x3, y3;
			this->All_Cells[nn - 1]->Get_Center(x3, y3);

			double x4, y4;
			this->All_Cells[nn + 2]->Get_Center(x4, y4);


			fout << 0.5 * (x + x2) << " " << linear(x3, this->All_Cells[nn - 1]->par[0].ro, x, i->par[0].ro, 0.5 * (x + x2)) << " " << 3 << endl;
			double ddd = linear(x2, this->All_Cells[nn + 1]->par[0].ro, x4, this->All_Cells[nn + 2]->par[0].ro, 0.5 * (x + x2));
			fout << 0.5 * (x + x2) << " " << ddd << " " << 4 << endl;
		}

		

	}
	fout.close();


	name_f = "1D_tecplot_source_" + name0 + ".txt";
	fout.open(name_f);
	fout << "TITLE = \"HP\"  VARIABLES = \"X\", \"M_21\", \"M_22\", \"M_3\", \"I_21\", \"I_22\", \"I_3\", " << //
		"ZONE T = \"HP\"" << endl;

	for (int nn = 0; nn < this->M1 + this->M2 + this->M3 + this->M4; nn++)
	{
		auto i = this->All_Cells[nn];
		double x, y;
		i->Get_Center(x, y);
		double a1, a2, a3, u, v, ro, p;
		ro = i->par[0].ro;
		p = i->par[0].p;
		u = i->par[0].u;
		v = i->par[0].v;
		//i->Get_Sourse_MK1(a1, a2, a3, u, v, ro, p, true);

		fout << x << " " << a1 << " " << a2 << " " << a3 << " " << i->par[0].I_u << " " << i->par[0].I_v << " " << i->par[0].I_T << endl;

	}
	fout.close();


	// Выводим каждый источник отдельно  -------------------------------------=-=-========================================================================

	name_f = "1D_SOURSE_Multifluid_" + name0 + ".txt";
	fout.open(name_f);
	fout << "TITLE = \"HP\"  VARIABLES = \"X\", \"I1_u\", \"I1_v\", \"I1_T\", \"I2_u\", \"I2_v\", \"I2_T\", \"I3_u\", \"I3_v\", \"I3_T\", \"I4_u\", \"I4_v\", \"I4_T\", \"I5_u\", \"I5_v\", \"I5_T\", \"I6_u\", \"I6_v\", \"I6_T\"" << //
		"ZONE T = \"HP\"" << endl;

	double u_H[pop_atom];
	double v_H[pop_atom];
	double ro_H[pop_atom];
	double p_H[pop_atom];
	double T_H[pop_atom];
	double U_M_H[pop_atom];
	double U_H[pop_atom];
	double sigma_H[pop_atom];
	double nu_H[pop_atom];

	num = 0;
	for (auto& KK : this->All_Cells)
	{
		num++;
		if (num > this->M1 + this->M2 + this->M3 + this->M4)
		{
			break;
		}

		double x, y;
		KK->Get_Center(x, y);

		double u, v, p, ro;
		u = KK->par[0].u;
		v = KK->par[0].v;
		p = KK->par[0].p;
		ro = KK->par[0].ro;

		for (int i = 0; i < pop_atom; i++)
		{
			u_H[i] = KK->par[0].H_u[i];
			v_H[i] = KK->par[0].H_v[i];
			ro_H[i] = KK->par[0].H_n[i];
			T_H[i] = KK->par[0].H_T[i];
			p_H[i] = 0.5 * T_H[i] * ro_H[i];
		}

		for (int i = 0; i < pop_atom; i++)
		{
			if (ro_H[i] <= 0.000000000001)
			{
				ro_H[i] = 0.0000001;
				p_H[i] = 0.0;
			}
		}


		for (int i = 0; i < pop_atom; i++)
		{
			U_M_H[i] = sqrt(kv(u - u_H[i]) + kv(v - v_H[i]) + (64.0 / (9.0 * pi_)) //
				* (p / ro + 2.0 * p_H[i] / ro_H[i]));

			U_H[i] = sqrt(kv(u - u_H[i]) + kv(v - v_H[i]) + (4.0 / (pi_)) //
				* (p / ro + 2.0 * p_H[i] / ro_H[i]));

			sigma_H[i] = kv(1.0 - a_2 * log(U_M_H[i]));

			nu_H[i] = ro * ro_H[i] * U_M_H[i] * sigma_H[i];
		}


		fout << x << " ";

		for (int i = 0; i < pop_atom; i++)
		{
			fout << (n_p_LISM_ / Kn_) * nu_H[i] * (u_H[i] - u) << " " << (n_p_LISM_ / Kn_) * nu_H[i] * (v_H[i] - v) << " " << (n_p_LISM_ / Kn_) * nu_H[i] * ((kv(u_H[i]) + kv(v_H[i]) - kv(u) - kv(v)) / 2.0 + //
				(U_H[i] / U_M_H[i]) * (2.0 * p_H[i] / ro_H[i] - p / ro)) << " ";
		}

		fout << endl;
	}

	fout.close();

	// --------------------------------------------------------------------------------------------------------------


	ofstream f_weyght;
	f_weyght.open("weyght.txt");
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < I_ + 1; j++)
		{
				f_weyght << i << " " << j <<  " " << kv(Rmax_/ Ri[j]) * Mu_stat[i][j]/(this->AllNumber) << " " << I_stat[i][j] << endl;
		}
	}
	f_weyght.close();
}

void Setka::Proverka(void)
{
	//cout << "Setka.cpp    " << "Proverka  Start" << endl;
	// Выписываем ячейки у которых не 4 соседа - особые
	for (auto& i : this->All_Cells)  
	{
		if (i->Grans.size() != 4)
		{
			double x, y;
			i->Get_Center(x, y);
			cout << "Setka.cpp    " << "Gran  82753   in  point:  " << x << " " <<  y << "   Have  sosed: " << i->Grans.size() << endl;
		}
	}


	// Проверяем правильное определение нормалей к граням
	for (auto& i : this->All_Gran)
	{
		if (i->type == Usualy)
		{
			double xc, yc;
			double n1, n2;
			double x2, y2;
			double x, y;
			i->Master->Get_Center(x, y);
			i->Get_Center(xc, yc);
			i->Sosed->Get_Center(x2, y2);
			i->Get_normal(n1, n2);
			if (n1 * (x2 - xc) + n2 * (y2 - yc) < 0)
			{
				cout << "Setka.cpp    " << "ERROR  gran  2345656435434543234:  " << xc << " " << yc << endl;
			}
			if (n1 * (x - xc) + n2 * (y - yc) > 0)
			{
				cout << "Setka.cpp    " << "ERROR  gran  wecwcefwcfwcfwcwcfwf:  " << xc << " " << yc << endl;
			}
		}

	}

	//cout << "Setka.cpp    " << "Proverka TVD" << endl;

	if (false)
	{
		for (auto& i : this->All_Gran)
		{
			if (i->type == Usualy)
			{
				Parametr par;
				i->Get_par_TVD(par, 0);
				if ((par.ro < i->Master->par[0].ro && par.ro < i->Sosed->par[0].ro) || (par.ro > i->Master->par[0].ro && par.ro > i->Sosed->par[0].ro))
				{
					cout << "Setka.cpp    " << "PROBLEM BIG" << endl;
				}
				if ((par.p < i->Master->par[0].p && par.p < i->Sosed->par[0].p) || (par.p > i->Master->par[0].p && par.p > i->Sosed->par[0].p))
				{
					cout << "Setka.cpp    " << "PROBLEM BIG" << endl;
				}
				if ((par.u < i->Master->par[0].u && par.u < i->Sosed->par[0].u) || (par.u > i->Master->par[0].u && par.u > i->Sosed->par[0].u))
				{
					cout << "Setka.cpp    " << "PROBLEM BIG" << endl;
				}
				if ((par.v < i->Master->par[0].v && par.v < i->Sosed->par[0].v) || (par.v > i->Master->par[0].v && par.v > i->Sosed->par[0].v))
				{
					cout << "Setka.cpp    " << "PROBLEM BIG" << endl;
				}
				if ((par.Q < i->Master->par[0].Q && par.Q < i->Sosed->par[0].Q) || (par.Q > i->Master->par[0].Q && par.Q > i->Sosed->par[0].Q))
				{
					cout << "Setka.cpp    " << "PROBLEM BIG" << endl;
				}
			}

		}
	}

	//cout << "Setka.cpp    " << "Proverka TVD 2" << endl;

	if (false)
	{
		for (auto& i : this->All_Gran_copy)
		{
			if (i->type == Usualy)
			{
				Parametr par;
				i->Get_par_TVD(par, 0);
				if ((par.ro < i->Master->par[0].ro && par.ro < i->Sosed->par[0].ro) || (par.ro > i->Master->par[0].ro && par.ro > i->Sosed->par[0].ro))
				{
					cout << "Setka.cpp    " << "PROBLEM BIG" << endl;
				}
			}

		}
	}
	
	//cout << "Setka.cpp    " << "Proverka TVD  3" << endl;

	// Проверяем расположение точек в ячейке по кругу!  (как грани расположены тоже надо будет проверить)
	for (auto& i : this->All_Cells)
	{
		bool t = false;
		for (int j = 0; j < i->contour.size() - 1; j++)
		{
			auto A = i->contour[j];
			auto B = i->contour[j + 1];
			for (auto& k : i->Grans)
			{
				if ((k->A == A && k->B == B) || (k->A == B && k->B == A))
				{
					t = true;
				}
			}
		}

		if (t == false)
		{
			cout << "Setka.cpp    " << "EROR hygwyfwuhlduwhguydcw" << endl;
			for (auto& j: i->contour)
			{
				cout << "Setka.cpp    " << j->x << " " << j->y << endl;
			}
			cout << "Setka.cpp    " << "EROR hygwyfwuhlduwhguydcw" << endl;
		}
	}


	//cout << "Setka.cpp    " << "Proverka  End" << endl;
}

void Setka::Copy(Setka* S)
{
	for (auto& i : this->A_Rails)
	{
		for (auto& j : S->A_Rails)
		{
			if (i->s <= j->s)
			{
				//cout << "Setka.cpp    " << i->s << " " << j->s << endl;
				i->Key_point[0]->Vx = j->Key_point[0]->x - i->Key_point[0]->x;
				i->Key_point[0]->Vy = j->Key_point[0]->y - i->Key_point[0]->y;
				i->Key_point[1]->Vx = j->Key_point[1]->x - i->Key_point[1]->x;
				i->Key_point[1]->Vy = j->Key_point[1]->y - i->Key_point[1]->y;
				i->Key_point[2]->Vx = j->Key_point[2]->x - i->Key_point[2]->x;
				i->Key_point[2]->Vy = j->Key_point[2]->y - i->Key_point[2]->y;
				break;
			}
		}
	}

	bool jk = false;
	Rail* DD = S->B_Rails[0];
	for (auto& i : this->B_Rails)
	{
		jk = false;
		//cout << "Setka.cpp    " << i->s << endl;
		for (auto& j : S->B_Rails)
		{
			if (i->s - 0.1 <= j->s)
			{
				jk = true;
				DD = j;
				//cout << "Setka.cpp    " << i->s << " " << j->s << endl;
				i->Key_point[0]->Vx = j->Key_point[0]->x - i->Key_point[0]->x;
				i->Key_point[0]->Vy = j->Key_point[0]->y - i->Key_point[0]->y;
				i->Key_point[1]->Vx = j->Key_point[1]->x - i->Key_point[1]->x;
				i->Key_point[1]->Vy = j->Key_point[1]->y - i->Key_point[1]->y;
				i->Key_point[2]->Vx = j->Key_point[2]->x - i->Key_point[2]->x;
				i->Key_point[2]->Vy = j->Key_point[2]->y - i->Key_point[2]->y;
				break;
			}
		}

		if (jk == false)
		{
			i->Key_point[0]->Vx = DD->Key_point[0]->x - i->Key_point[0]->x;
			i->Key_point[0]->Vy = DD->Key_point[0]->y - i->Key_point[0]->y;
			i->Key_point[1]->Vx = DD->Key_point[1]->x - i->Key_point[1]->x;
			i->Key_point[1]->Vy = DD->Key_point[1]->y - i->Key_point[1]->y;
			i->Key_point[2]->Vx = DD->Key_point[2]->x - i->Key_point[2]->x;
			i->Key_point[2]->Vy = DD->Key_point[2]->y - i->Key_point[2]->y;
		}
	}

	for (auto& i : this->C_Rails)
	{
		for (auto& j : S->C_Rails)
		{
			if (i->s - 0.001 <= j->s)
			{
				//cout << "Setka.cpp    " << i->s << " " << j->s << endl;
				i->Key_point[0]->Vx = j->Key_point[0]->x - i->Key_point[0]->x;
				i->Key_point[0]->Vy = j->Key_point[0]->y - i->Key_point[0]->y;
				break;
			}
		}
	}

	for (auto& i : this->D_Rails)
	{
		//cout << "Setka.cpp    " << i->s << endl;
		for (auto& j : S->D_Rails)
		{
			if (i->s + 1.0 >= j->s)
			{
				//cout << "Setka.cpp    " << i->s << " " << j->s << endl;
				i->Key_point[1]->Vx = j->Key_point[1]->x - i->Key_point[1]->x;
				i->Key_point[1]->Vy = j->Key_point[1]->y - i->Key_point[1]->y;
				i->Key_point[2]->Vx = j->Key_point[2]->x - i->Key_point[2]->x;
				i->Key_point[2]->Vy = j->Key_point[2]->y - i->Key_point[2]->y;
				break;
			}
		}
	}

	this->Move_Setka_Calculate_2(1.0);

	for (auto& i : this->All_Points)
	{
		i->x = i->x2;
		i->y = i->y2;
		i->Vx = 0.0;
		i->Vy = 0.0;
		i->count = 0;
	}

	this->Line_Outer[0]->A->Vx = this->Line_Outer[0]->B->Vx + (this->Line_Outer[0]->B->x - this->Line_Outer[0]->A->x);
	this->Line_Outer[this->Line_Outer.size() - 1]->B->Vy = this->Line_Outer[this->Line_Outer.size() - 1]->A->y - this->Line_Outer[this->Line_Outer.size() - 1]->B->y;
	this->Line_Contact[0]->A->Vx = this->Line_Contact[0]->B->Vx + (this->Line_Contact[0]->B->x - this->Line_Contact[0]->A->x);

	this->Line_Contact[this->Line_Contact.size() - 1]->B->Vy = this->Line_Contact[this->Line_Contact.size() - 1]->A->Vy + //
		(this->Line_Contact[this->Line_Contact.size() - 1]->A->y - this->Line_Contact[this->Line_Contact.size() - 1]->B->y);
	this->Move_Setka_Calculate_2(1.0);
	
	for (auto& i : this->All_Points)
	{
		i->x = i->x2;
		i->y = i->y2;
		i->Vx = 0.0;
		i->Vy = 0.0;
		i->count = 0;
	}

	cout << "Setka.cpp    " << "Filling data" << endl;
	double x, y;
	for (auto& i : this->All_Cells)
	{
		i->Get_Center(x, y);
		for (auto& j : S->All_Cells)
		{
			if (j->belong(x, y) == true)
			{
				i->par[0] = j->par[0];
				i->par[1] = j->par[1];
				break;
			}
		}
	}
}

void Setka::Save_G_D(void)
{
	ofstream fout;
	fout.open("Setka_all.txt");

	for (auto& i : this->All_Cells)
	{
		fout << i->par[0].ro << " " << i->par[0].p << " " << i->par[0].u << " " << i->par[0].v << " " << i->par[0].Q << endl;
	}
}

void Setka::Save_G_D_5_komponent(void)
{
	ofstream fout;
	fout.open("Setka_all_5_component.txt");

	for (auto& i : this->All_Cells)
	{
		fout << i->par[0].ro << " " << i->par[0].p << " " << i->par[0].u << " " << i->par[0].v << " " << i->par[0].Q << " " << //
			i->par[0].ro_H1 << " " << i->par[0].p_H1 << " " << i->par[0].u_H1 << " " << i->par[0].v_H1 << " " <<//
			i->par[0].ro_H2 << " " << i->par[0].p_H2 << " " << i->par[0].u_H2 << " " << i->par[0].v_H2 << " " <<//
			i->par[0].ro_H3 << " " << i->par[0].p_H3 << " " << i->par[0].u_H3 << " " << i->par[0].v_H3 << " " <<//
			i->par[0].ro_H4 << " " << i->par[0].p_H4 << " " << i->par[0].u_H4 << " " << i->par[0].v_H4 << " " <<//
			endl;
	}
}

void Setka::Save_Source_MK(string name)
{
	ofstream fout;
	fout.open(name);

	for (auto& i : this->All_Cells)
	{
		fout << i->par[0].I_u << " " << i->par[0].I_v << " " << i->par[0].I_T << endl;
		fout << i->par[0].II_u << " " << i->par[0].II_v << " " << i->par[0].II_T << endl;
		for (int j = 0; j < pop_atom; j++)
		{
			fout << i->par[0].H_n[j] << " " << i->par[0].H_u[j] << " " << i->par[0].H_v[j] << " " << i->par[0].H_T[j] << " " //
				<< i->par[0].H_uu[j] << " " << i->par[0].H_uv[j] << " " << i->par[0].H_vv[j] << " " << i->par[0].H_uuu[j] <<  endl;
		}
		fout << i->par[0].k_u << " " << i->par[0].k_v << " " << i->par[0].k_T << endl;
	}

	fout.close();
}

void Setka::Download_Source_MK(string name)
{
	ifstream fout;
	fout.open(name);

	if (fout.is_open() == false)
	{
		cout << "ERROR open  " << name << endl;
		exit(-100);
	}

	double a, b, c, d, e, f, g, h;

	for (auto& i : this->All_Cells)
	{
		fout >> a >> b >> c;
		i->par[0].I_u = a;
		i->par[0].I_v = b;
		i->par[0].I_T = c;
		//cout << "Istok = " << a << " " << b << " " << c << endl;
		fout >> a >> b >> c;
		i->par[0].II_u = a;
		i->par[0].II_v = b;
		i->par[0].II_T = c;

		for (int j = 0; j < pop_atom; j++)
		{
			i->par[0].H_n[j] = 0.0;
			i->par[0].H_u[j] = 0.0;
			i->par[0].H_v[j] = 0.0;
			i->par[0].H_T[j] = 0.0;
			i->par[0].H_uu[j] = 0.0;
			i->par[0].H_uv[j] = 0.0;
			i->par[0].H_vv[j] = 0.0;
			i->par[0].H_uuu[j] = 0.0;
		}

		for (int j = 0; j < pop_atom; j++)
		{
			fout >> a >> b >> c >> d >> e >> f >> g >> h;
			i->par[0].H_n[j] += a;
			i->par[0].H_u[j] += b;
			i->par[0].H_v[j] += c;
			i->par[0].H_T[j] += d;
			i->par[0].H_uu[j] += e;
			i->par[0].H_uv[j] += f;
			i->par[0].H_vv[j] += g;
			i->par[0].H_uuu[j] += h;
		}

		fout >> a >> b >> c;
		//fout >> i->par[0].k_u >> i->par[0].k_v >> i->par[0].k_T;
		i->par[0].k_u = a;
		i->par[0].k_v = b;
		i->par[0].k_T = c;
		//cout << "koeff = " << a << " " << b << " " << c << endl;
	}

	fout.close();
}

void Setka::Save_Setka_ALL_ALPHA(string name)
{
	int num = 0;
	for (auto& i : this->All_Gran)  
	{
		i->number = num;
		num++;
	}
	num = 0;
	for (auto& i : this->All_Gran_copy)  
	{
		i->number = num;
		num++;
	}


	ofstream fout;
	fout.open(name);

	fout << this->N1 << " " << this->N2 << " " << this->N3 << " " << this->N4 << " " << this->M1 << " " << this->M2 << " " << this->M3 << " " << this->M4 << " " << endl;
	
	fout << this->All_Points.size() << endl;
	for (auto& i : this->All_Points)
	{
		fout << i->x << " " << i->y << " " << i->type << endl;
	}

	fout << this->All_Gran.size() << endl;
	for (auto& i : this->All_Gran)
	{
		fout << i->A->number << " " << i->B->number << " " << i->type << endl;
	}

	fout << this->All_Gran_copy.size() << endl;
	for (auto& i : this->All_Gran_copy)
	{
		fout << i->A->number << " " << i->B->number << " " << i->type << endl;
	}

	fout << this->A_Rails.size() << endl;
	for (auto& i : this->A_Rails)
	{
		fout << i->M1 << " " << i->M2 << " " << i->M3 << " " << i->M4 << " " << i->s << " " << i->type << endl;
		fout << i->All_point.size() << endl;
		for (auto& j : i->All_point)
		{
			fout << j->number << endl;
		}

		fout << i->Key_point.size() << endl;
		for (auto& j : i->Key_point)
		{
			fout << j->number << endl;
		}
	}

	fout << this->B_Rails.size() << endl;
	for (auto& i : this->B_Rails)
	{
		fout << i->M1 << " " << i->M2 << " " << i->M3 << " " << i->M4 << " " << i->s << " " << i->type << endl;
		fout << i->All_point.size() << endl;
		for (auto& j : i->All_point)
		{
			fout << j->number << endl;
		}

		fout << i->Key_point.size() << endl;
		for (auto& j : i->Key_point)
		{
			fout << j->number << endl;
		}
	}

	fout << this->C_Rails.size() << endl;
	for (auto& i : this->C_Rails)
	{
		fout << i->M1 << " " << i->M2 << " " << i->M3 << " " << i->M4 << " " << i->s << " " << i->type << endl;
		fout << i->All_point.size() << endl;
		for (auto& j : i->All_point)
		{
			fout << j->number << endl;
		}

		fout << i->Key_point.size() << endl;
		for (auto& j : i->Key_point)
		{
			fout << j->number << endl;
		}
	}

	fout << this->D_Rails.size() << endl;
	for (auto& i : this->D_Rails)
	{
		fout << i->M1 << " " << i->M2 << " " << i->M3 << " " << i->M4 << " " << i->s << " " << i->type << endl;
		fout << i->All_point.size() << endl;
		for (auto& j : i->All_point)
		{
			fout << j->number << endl;
		}

		fout << i->Key_point.size() << endl;
		for (auto& j : i->Key_point)
		{
			fout << j->number << endl;
		}
	}

	fout << this->All_Cells.size() << endl;
	for (auto& i : this->All_Cells)
	{
		fout << i->type << endl;
		fout << i->contour.size() << endl;
		for (auto& j : i->contour)
		{
			fout << j->number << endl;
		}

		fout << i->Grans.size() << endl;
		for (auto& j : i->Grans)
		{
			fout << j->number << " " << j->main_gran << endl;
		}
	}

	fout << this->All_Gran.size() << endl;
	for (auto& i : this->All_Gran)
	{
		int i1, i2, i3;
		if (i->Sosed != nullptr)
		{
			i1 = i->Sosed->number;
		}
		else
		{
			i1 = -1;
		}

		if (i->Gran_copy != nullptr)
		{
			i2 = i->Gran_copy->number;
		}
		else
		{
			i2 = -1;
		}

		fout << i->Master->number << " " << i1 << " " << i2 << endl;
	}

	fout << this->All_Gran_copy.size() << endl;
	for (auto& i : this->All_Gran_copy)
	{
		int i_1, i_2, i_3;
		if (i->Master != nullptr)
		{
			i_1 = i->Master->number;
		}
		else
		{
			i_1 = -1;
		}

		if (i->Sosed != nullptr)
		{
			i_2 = i->Sosed->number;
		}
		else
		{
			i_2 = -1;
		}

		fout << i_1 << " " << i_2 << " " << i->Gran_copy->number << endl;
	}

	fout << this->All_Points.size() << endl;
	for (auto& i : this->All_Points)
	{
		fout << i->my_cell.size() << endl;
		for (auto& j : i->my_cell)
		{
			fout << j->number << endl;
		}
	}

	fout << this->Line_Contact.size() << endl;
	for (auto& i : this->Line_Contact)
	{
		fout << i->number << endl;
	}

	fout << this->Line_Inner.size() << endl;
	for (auto& i : this->Line_Inner)
	{
		fout << i->number << endl;
	}

	fout << this->Line_Outer.size() << endl;
	for (auto& i : this->Line_Outer)
	{
		fout << i->number << endl;
	}

	fout << this->All_Cells.size() << endl;
	for (auto& i : this->All_Cells)
	{
		fout << i->par[0].ro << " " << i->par[0].p << " " << i->par[0].u << " " << i->par[0].v << " " << i->par[0].Q << endl;
		fout << i->par[0].ro_H1 << " " << i->par[0].p_H1 << " " << i->par[0].u_H1 << " " << i->par[0].v_H1 << endl;
		fout << i->par[0].ro_H2 << " " << i->par[0].p_H2 << " " << i->par[0].u_H2 << " " << i->par[0].v_H2 << endl;
		fout << i->par[0].ro_H3 << " " << i->par[0].p_H3 << " " << i->par[0].u_H3 << " " << i->par[0].v_H3 << endl;
		fout << i->par[0].ro_H4 << " " << i->par[0].p_H4 << " " << i->par[0].u_H4 << " " << i->par[0].v_H4 << endl;
		fout << i->par[0].I_u << " " << i->par[0].I_v << " " << i->par[0].I_T << endl;
		for (int j = 0; j < 4; j++)
		{
			fout << i->par[0].H_n[j] << " " << i->par[0].H_u[j] << " " << i->par[0].H_v[j] << " " << i->par[0].H_T[j] << endl;
		}
		fout << i->par[0].k_u << " " << i->par[0].k_v << " " << i->par[0].k_T << endl;
	}

	fout.close();
}

void Setka::Download_Setka_ALL_ALPHA(string name)
{
	int n, m, k, type;
	double x, y;
	int num;
	int m1, m2, m3, m4;
	double s;
	bool b;
	ifstream fout;
	fout.open(name);

	fout >> this->N1 >> this->N2 >> this->N3 >> this->N4 >> this->M1 >> this->M2 >> this->M3 >> this->M4;

	fout >> n;
	this->All_Points.reserve(n);
	for (int i = 0; i < n; i++)
	{
		fout >> x >> y; 
		auto P = new Point(x, y);
		fout >> type;
		P->type = static_cast<Point_type>(type);
		this->All_Points.push_back(P);
	}

	fout >> n;
	this->All_Gran.reserve(n);
	for (int i = 0; i < n; i++)
	{
		fout >> m >> k;
		auto G = new Gran(this->All_Points[m], this->All_Points[k]);
		G->main_gran = true;
		fout >> type;
		G->type = static_cast<Gran_type>(type);
		this->All_Gran.push_back(G);
	}

	fout >> n;
	this->All_Gran_copy.reserve(n);
	for (int i = 0; i < n; i++)
	{
		fout >> m >> k;
		auto G = new Gran(this->All_Points[m], this->All_Points[k]);
		G->main_gran = false;
		fout >> type;
		G->type = static_cast<Gran_type>(type);
		this->All_Gran_copy.push_back(G);
	}

	fout >> n;
	this->A_Rails.reserve(n);
	for (int i = 0; i < n; i++)
	{
		fout >> m1 >> m2 >> m3 >> m4 >> s >> type;
		auto R = new Rail(s);
		R->M1 = m1;
		R->M2 = m2;
		R->M3 = m3;
		R->M4 = m4;
		R->type = static_cast<Rail_type>(type);

		fout >> m;
		R->All_point.reserve(m);
		for (int j = 0; j < m; j++)
		{
			fout >> k;
			R->All_point.push_back(this->All_Points[k]);
		}

		fout >> m;
		R->Key_point.reserve(m);
		for (int j = 0; j < m; j++)
		{
			fout >> k;
			R->Key_point.push_back(this->All_Points[k]);
		}
		this->A_Rails.push_back(R);
	}

	fout >> n;
	this->B_Rails.reserve(n);
	for (int i = 0; i < n; i++)
	{
		fout >> m1 >> m2 >> m3 >> m4 >> s >> type;
		auto R = new Rail(s);
		R->M1 = m1;
		R->M2 = m2;
		R->M3 = m3;
		R->M4 = m4;
		R->type = static_cast<Rail_type>(type);

		fout >> m;
		R->All_point.reserve(m);
		for (int j = 0; j < m; j++)
		{
			fout >> k;
			R->All_point.push_back(this->All_Points[k]);
		}

		fout >> m;
		R->Key_point.reserve(m);
		for (int j = 0; j < m; j++)
		{
			fout >> k;
			R->Key_point.push_back(this->All_Points[k]);
		}
		this->B_Rails.push_back(R);
	}

	fout >> n;
	this->C_Rails.reserve(n);
	for (int i = 0; i < n; i++)
	{
		fout >> m1 >> m2 >> m3 >> m4 >> s >> type;
		auto R = new Rail(s);
		R->M1 = m1;
		R->M2 = m2;
		R->M3 = m3;
		R->M4 = m4;
		R->type = static_cast<Rail_type>(type);

		fout >> m;
		R->All_point.reserve(m);
		for (int j = 0; j < m; j++)
		{
			fout >> k;
			R->All_point.push_back(this->All_Points[k]);
		}

		fout >> m;
		R->Key_point.reserve(m);
		for (int j = 0; j < m; j++)
		{
			fout >> k;
			R->Key_point.push_back(this->All_Points[k]);
		}
		this->C_Rails.push_back(R);
	}

	fout >> n;
	this->D_Rails.reserve(n);
	for (int i = 0; i < n; i++)
	{
		fout >> m1 >> m2 >> m3 >> m4 >> s >> type;
		auto R = new Rail(s);
		R->M1 = m1;
		R->M2 = m2;
		R->M3 = m3;
		R->M4 = m4;
		R->type = static_cast<Rail_type>(type);

		fout >> m;
		R->All_point.reserve(m);
		for (int j = 0; j < m; j++)
		{
			fout >> k;
			R->All_point.push_back(this->All_Points[k]);
		}

		fout >> m;
		R->Key_point.reserve(m);
		for (int j = 0; j < m; j++)
		{
			fout >> k;
			R->Key_point.push_back(this->All_Points[k]);
		}
		this->D_Rails.push_back(R);
	}

	fout >> n;
	this->All_Cells.reserve(n);
	for (int i = 0; i < n; i++)
	{
		auto C = new Cell();
		fout >> type;
		C->type = static_cast<Cell_type>(type);

		fout >> m;
		C->contour.reserve(m);
		for (int j = 0; j < m; j++)
		{
			fout >> k;
			C->contour.push_back(this->All_Points[k]);
		}

		fout >> m;
		C->Grans.reserve(m);
		for (int j = 0; j < m; j++)
		{
			fout >> k >> b;
			if (b == true)
			{
				C->Grans.push_back(this->All_Gran[k]);
			}
			else
			{
				C->Grans.push_back(this->All_Gran_copy[k]);
			}
		}

		this->All_Cells.push_back(C);
	}

	fout >> n;
	for (auto& i : this->All_Gran)
	{
		fout >> m1 >> m2 >> m3;
		i->Master = this->All_Cells[m1];
		if (m2 != -1)
		{
			i->Sosed = this->All_Cells[m2];
			i->Gran_copy = this->All_Gran_copy[m3];
		}
		else
		{
			i->Sosed = nullptr;
			i->Gran_copy = nullptr;
		}
	}

	fout >> n;
	for (auto& i : this->All_Gran_copy)
	{
		fout >> m1 >> m2 >> m3;
		i->Master = this->All_Cells[m1];
		if (m2 != -1)
		{
			i->Sosed = this->All_Cells[m2];
			i->Gran_copy = this->All_Gran[m3];
		}
		else
		{
			cout << "Setka.cpp    " << "Syda ne dolgny popadat  hrgrfwgydwfy2e2443" << endl;
			i->Sosed = nullptr;
			i->Gran_copy = nullptr;
		}
	}


	fout >> n;
	for (auto& i : this->All_Points)
	{
		fout >> m;
		i->my_cell.reserve(m);
		for (int j = 0; j < m; j++)
		{
			fout >> k;
			i->my_cell.push_back(this->All_Cells[k]);
		}
	}

	fout >> n;
	this->Line_Contact.reserve(n);
	for (int i = 0; i < n; i++)
	{
		fout >> k;
		this->Line_Contact.push_back(this->All_Gran[k]);
	}

	fout >> n;
	this->Line_Inner.reserve(n);
	for (int i = 0; i < n; i++)
	{
		fout >> k;
		this->Line_Inner.push_back(this->All_Gran[k]);
	}

	fout >> n;
	this->Line_Outer.reserve(n);
	for (int i = 0; i < n; i++)
	{
		fout >> k;
		this->Line_Outer.push_back(this->All_Gran[k]);
	}

	fout >> n;
	for (auto& i : this->All_Cells)
	{
		fout >> i->par[0].ro >> i->par[0].p >> i->par[0].u >> i->par[0].v >> i->par[0].Q;
		fout >> i->par[0].ro_H1 >> i->par[0].p_H1 >> i->par[0].u_H1 >> i->par[0].v_H1;
		fout >> i->par[0].ro_H2 >> i->par[0].p_H2 >> i->par[0].u_H2 >> i->par[0].v_H2;
		fout >> i->par[0].ro_H3 >> i->par[0].p_H3 >> i->par[0].u_H3 >> i->par[0].v_H3;
		fout >> i->par[0].ro_H4 >> i->par[0].p_H4 >> i->par[0].u_H4 >> i->par[0].v_H4;

		i->par[1] = i->par[0];
	}

	// Нумерация всего
	num = 0;
	for (auto& i : this->All_Gran)  // Нужно для выделения контактной поверхности
	{
		i->number = num;
		num++;
	}

	num = 0;
	for (auto& i : this->All_Gran_copy)  // Нужно для выделения контактной поверхности
	{
		i->number = num;
		num++;
	}

	num = 0;
	for (auto& i : this->All_Cells)
	{
		i->number = num;
		num++;
	}

	num = 0;
	// Нумеруем точки
	for (auto& i : this->All_Points)
	{
		i->number = num;
		num++;
	}


	cout << "Setka.cpp    " << "Vsego tochek = " << this->All_Points.size() << endl;
	cout << "Setka.cpp    " << "Vsego Yacheek = " << this->All_Cells.size() << endl;

	fout.close();
}

void Setka::Download_Setka_ALL_ALPHA_2_0(string name)
{
	int n, m, k, type;
	double x, y;
	int num;
	int m1, m2, m3, m4;
	double s;
	bool b;
	ifstream fout;
	fout.open(name);
	if (fout.is_open() == false)
	{
		cout << "ERROR open  " << name << endl;
		exit(-100);
	}
	cout << "Setka.cpp    " << "Downloading start" << endl;
	fout >> this->N1 >> this->N2 >> this->N3 >> this->N4 >> this->M1 >> this->M2 >> this->M3 >> this->M4;
	double LL = 1.0; // / RR_; // 0.00256505;
	fout >> n;
	this->All_Points.reserve(n);
	for (int i = 0; i < n; i++)
	{
		fout >> x >> y;

		if (std::fpclassify(x) != FP_NORMAL && std::fpclassify(x) != FP_ZERO)
		{
			cout << x << "    ERROR  3109 hsugfugufwe" << endl;
			exit(-1);
		}

		if (std::fpclassify(y) != FP_NORMAL && std::fpclassify(y) != FP_ZERO)
		{
			cout << y << "    ERROR  3109 hsugfugufwesrfrsf" << endl;
			exit(-1);
		}

		auto P = new Point(x * LL, y * LL);
		fout >> type;
		P->type = static_cast<Point_type>(type);
		this->All_Points.push_back(P);
	}

	fout >> n;
	this->All_Gran.reserve(n);
	for (int i = 0; i < n; i++)
	{
		fout >> m >> k;
		auto G = new Gran(this->All_Points[m], this->All_Points[k]);
		G->main_gran = true;
		fout >> type;
		G->type = static_cast<Gran_type>(type);
		this->All_Gran.push_back(G);
	}

	fout >> n;
	this->All_Gran_copy.reserve(n);
	for (int i = 0; i < n; i++)
	{
		fout >> m >> k;
		auto G = new Gran(this->All_Points[m], this->All_Points[k]);
		G->main_gran = false;
		fout >> type;
		G->type = static_cast<Gran_type>(type);
		this->All_Gran_copy.push_back(G);
	}

	fout >> n;
	this->A_Rails.reserve(n);
	for (int i = 0; i < n; i++)
	{
		fout >> m1 >> m2 >> m3 >> m4 >> s >> type;
		auto R = new Rail(s);
		R->M1 = m1;
		R->M2 = m2;
		R->M3 = m3;
		R->M4 = m4;
		R->type = static_cast<Rail_type>(type);

		fout >> m;
		R->All_point.reserve(m);
		for (int j = 0; j < m; j++)
		{
			fout >> k;
			R->All_point.push_back(this->All_Points[k]);
		}

		fout >> m;
		R->Key_point.reserve(m);
		for (int j = 0; j < m; j++)
		{
			fout >> k;
			R->Key_point.push_back(this->All_Points[k]);
		}
		this->A_Rails.push_back(R);
	}

	fout >> n;
	this->B_Rails.reserve(n);
	for (int i = 0; i < n; i++)
	{
		fout >> m1 >> m2 >> m3 >> m4 >> s >> type;
		auto R = new Rail(s);
		R->M1 = m1;
		R->M2 = m2;
		R->M3 = m3;
		R->M4 = m4;
		R->type = static_cast<Rail_type>(type);

		fout >> m;
		R->All_point.reserve(m);
		for (int j = 0; j < m; j++)
		{
			fout >> k;
			R->All_point.push_back(this->All_Points[k]);
		}

		fout >> m;
		R->Key_point.reserve(m);
		for (int j = 0; j < m; j++)
		{
			fout >> k;
			R->Key_point.push_back(this->All_Points[k]);
		}
		this->B_Rails.push_back(R);
	}

	fout >> n;
	this->C_Rails.reserve(n);
	for (int i = 0; i < n; i++)
	{
		fout >> m1 >> m2 >> m3 >> m4 >> s >> type;
		auto R = new Rail(s);
		R->M1 = m1;
		R->M2 = m2;
		R->M3 = m3;
		R->M4 = m4;
		R->type = static_cast<Rail_type>(type);

		fout >> m;
		R->All_point.reserve(m);
		for (int j = 0; j < m; j++)
		{
			fout >> k;
			R->All_point.push_back(this->All_Points[k]);
		}

		fout >> m;
		R->Key_point.reserve(m);
		for (int j = 0; j < m; j++)
		{
			fout >> k;
			R->Key_point.push_back(this->All_Points[k]);
		}
		this->C_Rails.push_back(R);
	}

	fout >> n;
	this->D_Rails.reserve(n);
	for (int i = 0; i < n; i++)
	{
		fout >> m1 >> m2 >> m3 >> m4 >> s >> type;
		auto R = new Rail(s);
		R->M1 = m1;
		R->M2 = m2;
		R->M3 = m3;
		R->M4 = m4;
		R->type = static_cast<Rail_type>(type);

		fout >> m;
		R->All_point.reserve(m);
		for (int j = 0; j < m; j++)
		{
			fout >> k;
			R->All_point.push_back(this->All_Points[k]);
		}

		fout >> m;
		R->Key_point.reserve(m);
		for (int j = 0; j < m; j++)
		{
			fout >> k;
			R->Key_point.push_back(this->All_Points[k]);
		}
		this->D_Rails.push_back(R);
	}
	cout << "Setka.cpp    " << "AAA1" << endl;
	fout >> n;
	this->All_Cells.reserve(n);
	for (int i = 0; i < n; i++)
	{
		auto C = new Cell();
		fout >> type;
		C->type = static_cast<Cell_type>(type);

		fout >> m;
		C->contour.reserve(m);
		for (int j = 0; j < m; j++)
		{
			fout >> k;
			C->contour.push_back(this->All_Points[k]);
		}

		fout >> m;
		C->Grans.reserve(m);
		for (int j = 0; j < m; j++)
		{
			fout >> k >> b;
			if (b == true)
			{
				C->Grans.push_back(this->All_Gran[k]);
			}
			else
			{
				if (this->All_Gran_copy.size() - 1 < k)
				{
					cout << "Setka.cpp    " << "ERROR  1741 " << endl;
					cout << "Setka.cpp    " << k << endl;
					exit(-2);
				}
				C->Grans.push_back(this->All_Gran_copy[k]);
			}
		}

		this->All_Cells.push_back(C);
	}
	cout << "Setka.cpp    " << "AAA2" << endl;
	fout >> n;
	for (auto& i : this->All_Gran)
	{
		fout >> m1 >> m2 >> m3;
		i->Master = this->All_Cells[m1];
		if (m2 != -1)
		{
			i->Sosed = this->All_Cells[m2];
		}
		else
		{
			i->Sosed = nullptr;
		}

		if (m3 != -1)
		{
			i->Gran_copy = this->All_Gran_copy[m3];
		}
		else
		{
			i->Gran_copy = nullptr;
		}
	}

	fout >> n;
	cout << "Setka.cpp    " << "AAA3" << endl;
	for (auto& i : this->All_Gran_copy)
	{
		fout >> m1 >> m2 >> m3;
		if (m1 != -1)
		{
			i->Master = this->All_Cells[m1];
		}
		else
		{
			i->Master = nullptr;
		}
		if (m2 != -1)
		{
			i->Sosed = this->All_Cells[m2];
			i->Gran_copy = this->All_Gran[m3];
		}
		else
		{
			cout << "Setka.cpp    " << "Syda ne dolgny popadat  hrgrfwgydwfy2e2443" << endl;
			i->Sosed = nullptr;
			i->Gran_copy = nullptr;
		}
	}


	fout >> n;
	for (auto& i : this->All_Points)
	{
		fout >> m;
		i->my_cell.reserve(m);
		for (int j = 0; j < m; j++)
		{
			fout >> k;
			i->my_cell.push_back(this->All_Cells[k]);
		}
	}

	fout >> n;
	this->Line_Contact.reserve(n);
	for (int i = 0; i < n; i++)
	{
		fout >> k;
		this->Line_Contact.push_back(this->All_Gran[k]);
	}

	fout >> n;
	this->Line_Inner.reserve(n);
	for (int i = 0; i < n; i++)
	{
		fout >> k;
		this->Line_Inner.push_back(this->All_Gran[k]);
	}

	fout >> n;
	this->Line_Outer.reserve(n);
	for (int i = 0; i < n; i++)
	{
		fout >> k;
		this->Line_Outer.push_back(this->All_Gran[k]);
	}

	fout >> n;
	for (auto& i : this->All_Cells)
	{
		fout >> i->par[0].ro >> i->par[0].p >> i->par[0].u >> i->par[0].v >> i->par[0].Q;
		fout >> i->par[0].ro_H1 >> i->par[0].p_H1 >> i->par[0].u_H1 >> i->par[0].v_H1;
		fout >> i->par[0].ro_H2 >> i->par[0].p_H2 >> i->par[0].u_H2 >> i->par[0].v_H2;
		fout >> i->par[0].ro_H3 >> i->par[0].p_H3 >> i->par[0].u_H3 >> i->par[0].v_H3;
		fout >> i->par[0].ro_H4 >> i->par[0].p_H4 >> i->par[0].u_H4 >> i->par[0].v_H4;

		fout >> i->par[0].I_u >> i->par[0].I_v >> i->par[0].I_T;

		for (int j = 0; j < 4; j++)
		{
			fout >> i->par[0].H_n[j] >> i->par[0].H_u[j] >> i->par[0].H_v[j] >> i->par[0].H_T[j];
		}

		fout >> i->par[0].k_u >> i->par[0].k_v >> i->par[0].k_T;

		i->par[1] = i->par[0];

		if (std::fpclassify(i->par[0].ro) != FP_NORMAL)
		{
			cout << i->par[0].ro << "    ERROR  3109 hsugfugufwesrfrsf i->par[0].ro" << endl;
			exit(-1);
		}
	}




	// Нумерация всего
	num = 0;
	for (auto& i : this->All_Gran)  // Нужно для выделения контактной поверхности
	{
		i->number = num;
		num++;
	}

	for (auto& i : this->All_Gran_copy)  // Нужно для выделения контактной поверхности
	{
		i->number = num;
		num++;
	}

	num = 0;
	for (auto& i : this->All_Cells)
	{
		i->number = num;
		num++;
	}

	num = 0;
	// Нумеруем точки
	for (auto& i : this->All_Points)
	{
		i->number = num;
		num++;
	}


	cout << "Setka.cpp    " << "Vsego tochek = " << this->All_Points.size() << endl;
	cout << "Setka.cpp    " << "Vsego Yacheek = " << this->All_Cells.size() << endl;
	cout << "Setka.cpp    " << "Downloading end" << endl;
	fout.close();

	for (auto& i : this->All_Cells)
	{
		i->renew();
	}
}

void Setka::Download_G_D_5_komponent(void)
{
	ifstream fout;
	fout.open("Setka_all_5_component.txt");

	for (auto& i : this->All_Cells)
	{
		fout >> i->par[0].ro >>i->par[0].p >>i->par[0].u >>i->par[0].v >>i->par[0].Q >>//
			i->par[0].ro_H1 >>i->par[0].p_H1 >>i->par[0].u_H1 >>i->par[0].v_H1 >>//
			i->par[0].ro_H2 >>i->par[0].p_H2 >>i->par[0].u_H2 >>i->par[0].v_H2 >>//
			i->par[0].ro_H3 >>i->par[0].p_H3 >>i->par[0].u_H3 >>i->par[0].v_H3 >>//
			i->par[0].ro_H4 >>i->par[0].p_H4 >>i->par[0].u_H4 >>i->par[0].v_H4;
	}
}

void Setka::Download_G_D(void)
{
	ifstream fout;
	fout.open("Setka_all.txt");

	for (auto& i : this->All_Cells)
	{
		fout >> i->par[0].ro >> i->par[0].p >> i->par[0].u >> i->par[0].v >> i->par[0].Q;
		i->par[1] = i->par[0];
	}
}

void Setka::Move_surface_hand(void)
{
	/*for (auto& i : this->Line_Contact)
	{
		i->A->Vx = 0.0;
		i->A->Vy = 0.0;
		i->B->Vx = 0.0;
		i->B->Vy = 0.0;
	}*/
	
	double x, y, d, dist;

	for (auto& i : this->Line_Contact)
	{
		d = 1000;
		x = i->A->x;
		y = i->A->y;
		for (auto& j : this->Contact)
		{
			dist = sqrt(kv(j->x - x) + kv(j->y - y));
			if (dist < d)
			{
				d = dist;
				i->A->Vx = (j->x - x) * 0.1;
				i->A->Vy = (j->y - y) * 0.1;
			}
		}
	}

	this->Line_Contact[0]->A->Vx = this->Line_Contact[0]->B->Vx + (this->Line_Contact[0]->B->x - this->Line_Contact[0]->A->x);

	this->Line_Contact[this->Line_Contact.size() - 1]->B->Vy = this->Line_Contact[this->Line_Contact.size() - 1]->A->Vy + //
		(this->Line_Contact[this->Line_Contact.size() - 1]->A->y - this->Line_Contact[this->Line_Contact.size() - 1]->B->y);

	double VY;

	for (auto& i : this->Line_Outer)
	{
		d = 1000;
		x = i->A->x;
		y = i->A->y;
		if (x > 0)
		{
			for (auto& j : this->Outer)
			{
				dist = sqrt(kv(j->x - x) + kv(j->y - y));
				if (dist < d)
				{
					d = dist;
					i->A->Vx = (j->x - x) * 0.1;
					i->A->Vy = (j->y - y) * 0.1;
					VY = i->A->Vy;
				}
			}
		}
		else
		{
			i->A->Vx = 0.0;
			i->A->Vy = VY;
		}
	}

	this->Line_Outer[0]->A->Vx = this->Line_Outer[0]->B->Vx + (this->Line_Outer[0]->B->x - this->Line_Outer[0]->A->x);

	this->Line_Outer[this->Line_Outer.size() - 1]->B->Vy = this->Line_Outer[this->Line_Outer.size() - 1]->A->Vy + //
		(this->Line_Outer[this->Line_Outer.size() - 1]->A->y - this->Line_Outer[this->Line_Outer.size() - 1]->B->y);

	for (auto& i : this->Line_Inner)
	{
		d = 1000;
		x = i->A->x;
		y = i->A->y;
		for (auto& j : this->Inner)
		{
			dist = sqrt(kv(j->x - x) + kv(j->y - y));
			if (dist < d)
			{
				d = dist;
				i->A->Vx = (j->x - x) * 0.1;
				i->A->Vy = (j->y - y) * 0.1;
			}
		}
	}

	this->Line_Inner[0]->A->Vx = this->Line_Inner[0]->B->Vx + (this->Line_Inner[0]->B->x - this->Line_Inner[0]->A->x);

	this->Line_Inner[this->Line_Inner.size() - 1]->B->Vx = this->Line_Inner[this->Line_Inner.size() - 1]->A->Vx + //
		(this->Line_Inner[this->Line_Inner.size() - 1]->A->x - this->Line_Inner[this->Line_Inner.size() - 1]->B->x);


}

void Setka::Move_surface(int ii, const double& dt = 1.0)
{
	double koef1 = 0.01 * 1.0; // 0.08; // 0.3;
	double koef2 = 0.1 * 0.05; // 0.005;
	double koef3 = 0.1 * 1.0; // 0.1;
	// Разбираемся с контактом

	//for (int j = 0; j < this->Line_Contact.size(); j++)  // Вычисляем скорость контакта
	//{
	//	auto i = this->Line_Contact[j];
	//	i->A->Vx = 0.3;
	//	i->A->Vy = 0.0;
	//	i->B->Vx = 0.3;
	//	i->B->Vy = 0.0;

	//}
	
	// Контакт
	if (true)
	{
		//int bb = -1;
		for (int j = 0; j < this->Line_Contact.size(); j++)  // Вычисляем скорость контакта
		{
			if (j == 0)
			{
				continue;
			}
			auto S = Solvers();
			auto i = this->Line_Contact[j];
			Parametr par1;// = i->Master->par[ii];
			Parametr par2;// = i->Sosed->par[ii];
			i->Get_par_TVD(par1, ii);
			i->Gran_copy->Get_par_TVD(par2, ii);
			double n1, n2;
			vector<double> qqq1(5);
			vector<double> qqq2(5);
			vector<double> qqq(5);
			vector<double> n(3);
			i->Get_normal(n1, n2);
			double VV, Vl, Vp;
			double P[4];
			double PQ;

			n[0] = n1;
			n[1] = n2;
			n[2] = 0.0;

			qqq1[0] = par1.ro;
			qqq1[1] = par1.p;
			qqq1[2] = par1.u;
			qqq1[3] = par1.v;
			qqq1[4] = 0.0;

			qqq2[0] = par2.ro;
			qqq2[1] = par2.p;
			qqq2[2] = par2.u;
			qqq2[3] = par2.v;
			qqq2[4] = 0.0;

			/*double x, y;

			i->Get_Center(x, y);*/

			//this->HLLC_2d_Korolkov_b_s(par1.ro, par1.Q, par1.p, par1.u, par1.v, par2.ro, par2.Q, //
			//	par2.p, par2.u, par2.v, 0.0, P, PQ, n1, n2, 1.0, 1, Vl, VV, Vp, false);

			
			S.Godunov_Solver_Alexashov(qqq1, qqq2, n, qqq, Vl, Vp, VV);
			

			double Max = max(sqrt((kv(par1.u) + kv(par1.v)) / (ggg * par1.p / par1.ro)), sqrt((kv(par2.u) + kv(par2.v)) / (ggg * par2.p / par2.ro)));
			VV = VV * koef1;// * 0.3;

			//VV = 0.2;

			/*if (i->A->x < -400)
			{
				VV = 0.0;
			}*/

			/*if (i->A->x < -200)
			{
				VV = par1.u * n1 + par1.v * n2;
			}*/



			/*if ( (i->A->x < -200 && i->A->y < 200 && VV <= 0)|| (i->B->x < -200 && i->B->y < 200 && VV <= 0) )
			{
				VV = 0.0;
			}*/

			double t1 = -n2;
			double t2 = n1;

			if (j == 1)//(Max < 1)
			{
				i->B->Vx += VV * n1;
				i->B->Vy += VV * n2;
				i->A->Vx += VV * n1;
				i->A->Vy += VV * n2;
				i->B->count++;
				i->A->count++;
			}
			else
			{
				i->B->Vx += VV * n1;
				i->B->Vy += VV * n2;
				i->B->count++;
			}

			/*if (par1.u * t1 + par1.v * t2 > 0)
			{
				i->B->Vx += VV * n1;
				i->B->Vy += VV * n2;
				i->B->count++;
			}
			else if (par1.u * t1 + par1.v * t2 < 0)
			{
				i->A->Vx += VV * n1;
				i->A->Vy += VV * n2;
				i->A->count++;
			}
			else
			{
				i->B->Vx += VV * n1;
				i->B->Vy += VV * n2;
				i->A->Vx += VV * n1;
				i->A->Vy += VV * n2;
				i->B->count++;
				i->A->count++;
			}*/

			/*if (i->A->x < -200 && i->A->y < 400 && i->B->x < -200 && i->B->y < 400)
			{
				VV = 0.05;
			}*/
			//cout << "Setka.cpp    " << i->A->x << " " << i->A->y << " " << VV << " " << n1 << " " << n2 <<  endl;
			/*if (Max > 1)
			{
				if (bb == -1)
				{
					bb = j;
				}
				i->B->Vx += VV * n1;
				i->B->Vy += VV * n2;
			}
			else
			{
				i->A->Vx += VV * n1 * 0.5;
				i->A->Vy += VV * n2 * 0.5;
				i->B->Vx += VV * n1 * 0.5;
				i->B->Vy += VV * n2 * 0.5;
			}*/
		}

		for (int j = 0; j < this->Line_Contact.size(); j++)  // Вычисляем скорость контакта
		{
			auto i = this->Line_Contact[j];
			if (i->A->count > 0)
			{
				i->A->Vx /= koef1 * (1.0 * i->A->count);
				i->A->Vy /= koef1 * (1.0 * i->A->count);
			}

		}

		this->Line_Contact[0]->A->Vx = this->Line_Contact[0]->B->Vx + (this->Line_Contact[0]->B->x - this->Line_Contact[0]->A->x)/dt;

		this->Line_Contact[this->Line_Contact.size() - 1]->B->Vy = this->Line_Contact[this->Line_Contact.size() - 1]->A->Vy + //
			(this->Line_Contact[this->Line_Contact.size() - 1]->A->y - this->Line_Contact[this->Line_Contact.size() - 1]->B->y)/dt;

		/*if (bb >= 0)
		{
			this->Line_Contact[bb]->A->Vx *= 2.0;
			this->Line_Contact[bb]->A->Vy *= 2.0;
		}*/

		// Делаем поверхостное натяжение
		for (int j = 0; j < this->Line_Contact.size(); j++)  // Вычисляем скорость контакта
		{
			if (j <= 0 || j >= this->Line_Contact.size() - 2)
			{
				continue;
			}

			auto i = this->Line_Contact[j];
			auto i1 = this->Line_Contact[j + 1];
			auto i2 = this->Line_Contact[j + 2];
			auto i3 = this->Line_Contact[j - 1];

			double x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, xx, yy, zz;

			x1 = i2->B->x;
			y1 = i2->B->y;
			x2 = i1->B->x;
			y2 = i1->B->y;
			x3 = i->B->x;
			y3 = i->B->y;
			x4 = i->A->x;
			y4 = i->A->y;
			x5 = i3->A->x;
			y5 = i3->A->y;

			//this->Smooth_kvadr3(x1, y1, 0.0, x2, y2, 0.0, x3, y3, 0.0, x4, y4, 0.0, x5, y5, 0.0, xx, yy, zz);
			xx = (x2  + x4) / 2.0;
			yy = (y2  + y4) / 2.0;
			i->B->Vx += koef1 * 0.1 * (xx - x3)/ dt;
			i->B->Vy += koef1 * 0.1 * (yy - y3) / dt;

		}
	}

	// Движение внутренней ударной волны
	if (true)
	{
		// Движение внутренней ударной волны
		for (int j = 0; j < this->Line_Inner.size(); j++)
		{
			/*if (j > this->N1 + n_inner_shock)
			{
				continue;
			}*/
			auto i = this->Line_Inner[j];
			Parametr par1 = i->Master->par[ii];
			double x, y, x2, y2;
			i->Master->Get_Center(x, y);
			double radius = sqrt(kv(x) + kv(y));
			i->Get_Center(x2, y2);
			double dis = sqrt(kv(x2) + kv(y2));
			Parametr par2;// = i->Sosed->par[ii];
			//par2 = i->Gran_copy->Master->par[ii];
			i->Gran_copy->Get_par_TVD(par2, ii);
			double n1, n2;
			i->Get_normal(n1, n2);
			double VV, Vl, Vp;
			double P[4];
			double PQ;
			par1.ro = par1.ro * kv(radius) / kv(dis);
			par1.p = par1.p * pow(radius / dis, 2.0 * ggg);
			polar_perenos(x, y, x2, y2, par1.u, par1.v);

			auto S = Solvers();
			vector<double> qqq1(5);
			vector<double> qqq2(5);
			vector<double> qqq(5);
			vector<double> n(3);

			n[0] = n1;
			n[1] = n2;
			n[2] = 0.0;

			qqq1[0] = par1.ro;
			qqq1[1] = par1.p;
			qqq1[2] = par1.u;
			qqq1[3] = par1.v;
			qqq1[4] = 0.0;

			qqq2[0] = par2.ro;
			qqq2[1] = par2.p;
			qqq2[2] = par2.u;
			qqq2[3] = par2.v;
			qqq2[4] = 0.0;

			S.Godunov_Solver_Alexashov(qqq1, qqq2, n, qqq, Vl, Vp, VV);

			//this->HLLC_2d_Korolkov_b_s(par1.ro, par1.Q, par1.p, par1.u, par1.v, par2.ro, par2.Q, //
			//	par2.p, par2.u, par2.v, 0.0, P, PQ, n1, n2, 1.0, 1, Vl, VV, Vp);
			Vl = Vl * koef2;  // 0.1
			double t1 = -n2;
			double t2 = n1;

			if (true)//(fabs(par2.u * t1 + par2.v * t2) < 0.1) // 0.1
			{
				//i->A->Vx += Vl * n1;
				//i->A->Vy += Vl * n2;
				//i->A->count++;
				i->B->Vx += Vl * n1;
				i->B->Vy += Vl * n2;
				i->B->count++;
			}
			else
			{
				i->A->Vx += Vl * n1;
				i->A->Vy += Vl * n2;
				i->A->count++;
			}

			/*if ( fabs(par2.u * t1 + par2.v * t2) < 0.01)
			{
				i->A->Vx += Vl * n1;
				i->A->Vy += Vl * n2;
				i->A->count++;
				i->B->Vx += Vl * n1;
				i->B->Vy += Vl * n2;
				i->B->count++;
			}
			else if (par2.u * t1 + par2.v * t2 > 0)
			{
				i->B->Vx += Vl * n1;
				i->B->Vy += Vl * n2;
				i->B->count++;
			}
			else
			{
				i->A->Vx += Vl * n1;
				i->A->Vy += Vl * n2;
				i->A->count++;
			}*/

			/*if (i->A->x < 0)
			{
				i->B->Vx += Vl * n1;
				i->B->Vy += Vl * n2;
				i->B->count++;
			}
			else
			{
				i->B->Vx += Vl * n1;
				i->B->Vy += Vl * n2;
				i->A->Vx += Vl * n1;
				i->A->Vy += Vl * n2;
				i->B->count++;
				i->A->count++;
			}*/

			//if (true)// (j > 12)
			//{
			//	if (par1.u * t1 + par1.v * t2 > 0)
			//	{
			//		i->B->Vx += Vl * n1;
			//		i->B->Vy += Vl * n2;
			//		i->B->count++;
			//	}
			//	else if (par1.u * t1 + par1.v * t2 < 0)
			//	{
			//		i->A->Vx += Vl * n1;
			//		i->A->Vy += Vl * n2;
			//		i->A->count++;
			//	}
			//	else
			//	{
			//		i->B->Vx += Vl * n1;
			//		i->B->Vy += Vl * n2;
			//		i->A->Vx += Vl * n1;
			//		i->A->Vy += Vl * n2;
			//		i->B->count++;
			//		i->A->count++;
			//	}
			//}
			/*else if(j < 12)
			{
				if (par1.u * t1 + par1.v * t2 < 0)
				{
					i->B->Vx += Vl * n1;
					i->B->Vy += Vl * n2;
					i->B->count++;
				}
				else if (par1.u * t1 + par1.v * t2 > 0)
				{
					i->A->Vx += Vl * n1;
					i->A->Vy += Vl * n2;
					i->A->count++;
				}
				else
				{
					i->B->Vx += Vl * n1;
					i->B->Vy += Vl * n2;
					i->A->Vx += Vl * n1;
					i->A->Vy += Vl * n2;
					i->B->count++;
					i->A->count++;
				}
			}
			else
			{
				i->B->Vx += Vl * n1;
				i->B->Vy += Vl * n2;
				i->A->Vx += Vl * n1;
				i->A->Vy += Vl * n2;
				i->B->count++;
				i->A->count++;
			}*/


		}

		for (int j = 0; j < this->Line_Inner.size(); j++)  // Вычисляем скорость контакта
		{
			auto i = this->Line_Inner[j];
			// Можно что-то ещё придумать с узлом где count = 0
			if (i->A->count > 0)
			{
				i->A->Vx /= (1.0 * i->A->count);
				i->A->Vy /= (1.0 * i->A->count);
			}
			/*else if(j > 0 && j < this->Line_Inner.size() - 1)
			{
				i->A->Vx = (i->B->count + this->Line_Contact[j - 1]->A->Vx)/2.0;
			}*/
		}

		this->Line_Inner[0]->A->Vx = this->Line_Inner[0]->B->Vx + (this->Line_Inner[0]->B->x - this->Line_Inner[0]->A->x)/dt;

		this->Line_Inner[this->Line_Inner.size() - 1]->B->Vx = this->Line_Inner[this->Line_Inner.size() - 1]->A->Vx + //
			(this->Line_Inner[this->Line_Inner.size() - 1]->A->x - this->Line_Inner[this->Line_Inner.size() - 1]->B->x)/dt;


		// Поверхностное натяжение
		// Делаем поверхостное натяжение
		for (int j = 0; j < this->Line_Inner.size(); j++)  // Вычисляем скорость контакта
		{
			if (j <= 0 || j >= this->Line_Inner.size() - 2)
			{
				continue;
			}

			auto i = this->Line_Inner[j];
			auto i1 = this->Line_Inner[j + 1];
			auto i2 = this->Line_Inner[j + 2];
			auto i3 = this->Line_Inner[j - 1];

			double x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, xx, yy, zz;

			x1 = i2->B->x;
			y1 = i2->B->y;
			x2 = i1->B->x;
			y2 = i1->B->y;
			x3 = i->B->x;
			y3 = i->B->y;
			x4 = i->A->x;
			y4 = i->A->y;
			x5 = i3->A->x;
			y5 = i3->A->y;

			//this->Smooth_kvadr3(x1, y1, 0.0, x2, y2, 0.0, x3, y3, 0.0, x4, y4, 0.0, x5, y5, 0.0, xx, yy, zz);
			xx = (x2 + x4) / 2.0;
			yy = (y2 + y4) / 2.0;
			i->B->Vx += koef2 * 0.0001 * (xx - x3) / dt;
			i->B->Vy += koef2 * 0.0001 * (yy - y3) / dt;

		}

	}

	// Двигаем ли внешнюю волну
	if (true)
	{
		//double Vconst = 0.0;
		for (int j = 0; j < this->Line_Outer.size(); j++)
		{
			if (j >= n_outer_shock + this->N1 + this->N2 + 1)
			{
				continue;
			}
			auto i = this->Line_Outer[j];
			Parametr par1;// = i->Master->par[ii];
			Parametr par2;// = i->Sosed->par[ii];
			i->Get_par_TVD(par1, ii);
			i->Gran_copy->Get_par_TVD(par2, ii);
			double n1, n2;
			i->Get_normal(n1, n2);
			double VV, Vl, Vp;
			double P[4];
			double PQ;

			vector<double> qqq1(5);
			vector<double> qqq2(5);
			vector<double> qqq(5);
			vector<double> n(3);
			n[0] = n1;
			n[1] = n2;
			n[2] = 0.0;

			qqq1[0] = par1.ro;
			qqq1[1] = par1.p;
			qqq1[2] = par1.u;
			qqq1[3] = par1.v;
			qqq1[4] = 0.0;

			qqq2[0] = par2.ro;
			qqq2[1] = par2.p;
			qqq2[2] = par2.u;
			qqq2[3] = par2.v;
			qqq2[4] = 0.0;

			//auto S = Solvers();
			//S.Godunov_Solver_Alexashov(qqq1, qqq2, n, qqq, Vl, Vp, VV);

			this->HLLC_2d_Korolkov_b_s(par1.ro, par1.Q, par1.p, par1.u, par1.v, par2.ro, par2.Q, //
				par2.p, par2.u, par2.v, 0.0, P, PQ, n1, n2, 1.0, 1, Vl, VV, Vp);
			Vp = Vp * koef3;  // 0.1
			//cout << "Setka.cpp    " << Vp << endl;

			double t1 = -n2;
			double t2 = n1;

			if (kvv(par1.u, par1.v, 0.0) / (ggg * par1.p / par1.ro) > 1.0)
			{
				i->B->Vx += Vp * n1;
				i->B->Vy += Vp * n2;
				i->B->count++;
				/*i->A->Vx += Vp * n1;
				i->A->Vy += Vp * n2;
				i->A->count++;*/
			}
			else
			{
				if (fabs(par1.u * t1 + par1.v * t2) < 0.01)
				{
					i->A->Vx += Vp * n1;
					i->A->Vy += Vp * n2;
					i->A->count++;
					i->B->Vx += Vp * n1;
					i->B->Vy += Vp * n2;
					i->B->count++;
				}
				else if (par1.u * t1 + par1.v * t2 > 0)
				{
					i->B->Vx += Vp * n1;
					i->B->Vy += Vp * n2;
					i->B->count++;
				}
				else
				{
					i->A->Vx += Vp * n1;
					i->A->Vy += Vp * n2;
					i->A->count++;
				}
			}
			
			
		}

		for (int j = 0; j < this->Line_Outer.size(); j++)  // Вычисляем скорость контакта
		{
			auto i = this->Line_Contact[j];
			if (i->A->count > 0)
			{
				i->A->Vx /= (1.0 * i->A->count);
				i->A->Vy /= (1.0 * i->A->count);
			}
		}

		this->Line_Outer[0]->A->Vx = this->Line_Outer[0]->B->Vx + (this->Line_Outer[0]->B->x - this->Line_Outer[0]->A->x)/dt;

		this->Line_Outer[this->Line_Outer.size() - 1]->B->Vy = this->Line_Outer[this->Line_Outer.size() - 1]->A->Vy + //
			(this->Line_Outer[this->Line_Outer.size() - 1]->A->y - this->Line_Outer[this->Line_Outer.size() - 1]->B->y)/dt;
	}
}

void Setka::Move_Setka_Calculate_3(const double& dt)
{
	// Для тестовой задачи с Дмитрием Борисовичем
	double Vx, Vy, V;
	double R2, r, R3, R4;

	//double y_Outer = this->D_Rails[16]->Key_point[1]->y;

	for (int jj = 0; jj < this->A_Rails.size(); jj++)
	{
		auto i = this->A_Rails[jj];
		// Подвинем ключевые точки

		double r = 90.0/RR_;

		i->Key_point[0]->x2 = r * cos(i->s);
		i->Key_point[0]->y2 = r * sin(i->s);

		//cout << "Setka.cpp    " << jj << " " << i->Key_point[1]->Vx << " " << i->Key_point[1]->Vy << endl;
		r = 130.0 / RR_ * 0.8 * (cos(i->s) + 1.5) / (cos(i->s) + 1.0);

		i->Key_point[1]->x2 = r * cos(i->s);
		i->Key_point[1]->y2 = r * sin(i->s);

		r = 245.0 / RR_ * 2.0 * (1.0) / (cos(i->s) + 1.0);

		i->Key_point[2]->x2 = r * cos(i->s);
		i->Key_point[2]->y2 = r * sin(i->s);

		i->Key_point[3]->x2 = i->Key_point[3]->x;
		i->Key_point[3]->y2 = i->Key_point[3]->y;


		R2 = sqrt(kv(i->Key_point[0]->x2) + kv(i->Key_point[0]->y2));

		//double x;

		//x = pow(R11_ / R1_, 1.0 / n_inner);
		//for (int j = 0; j <= n_inner; j++) 
		//{
		//	r = R1_ * pow(x, j);
		//	//r = R1_ + (R11_ - R1_) * j / (n_inner);
		//	i->All_point[j]->x2 = r * cos(i->s);
		//	i->All_point[j]->y2 = r * sin(i->s);
		//}


		//for (int j = n_inner + 1; j <= m_inner; j++) 
		//{
		//	r = R11_ + (R111_ - R11_) * (j - n_inner) / (m_inner - n_inner);
		//	i->All_point[j]->x2 = r * cos(i->s);
		//	i->All_point[j]->y2 = r * sin(i->s);
		//}

		//for (int j = m_inner + 1; j < i->M1 - 1; j++)  // Передвинули точки до ударной волны
		//{
		//	//r = R11_ * pow(x, j - 7);
		//	r = R111_ + (R2 - R111_) * (j - m_inner) / (i->M1 - 1 - m_inner);
		//	i->All_point[j]->x2 = r * cos(i->s);
		//	i->All_point[j]->y2 = r * sin(i->s);
		//}

		double x;

		//x = pow(R2 / R1_, 1.0 / (i->M1 - 1));
		for (int j = 0; j < n_inner; j++)
		{
			x = j * 1.0 / (n_inner);
			r = R1_ + pow(x, 1.3) * (R111_ - R1_);
			//r = R1_ + (R11_ - R1_) * j / (n_inner);
			i->All_point[j]->x2 = r * cos(i->s);
			i->All_point[j]->y2 = r * sin(i->s);
		}

		for (int j = n_inner; j < i->M1 - 1; j++)
		{
			x = (j - n_inner) * 1.0 / (i->M1 - 1 - n_inner);
			r = R111_ + pow(x, 1.05) * (R2 - R111_);
			//r = R1_ + (R11_ - R1_) * j / (n_inner);
			i->All_point[j]->x2 = r * cos(i->s);
			i->All_point[j]->y2 = r * sin(i->s);
		}



		//R3 = sqrt(kv(i->Key_point[1]->x2) + kv(i->Key_point[1]->y2)) - zazor;
		R3 = sqrt(kv(i->Key_point[1]->x2) + kv(i->Key_point[1]->y2));

		for (int j = 0; j < i->M2 - 1; j++)  // Передвинули точки до контакта
		{
			r = R2 + (R3 - R2) * (j + 1) / (i->M2);         // Равномерное распределение точек
			i->All_point[i->M1 + j]->x2 = r * cos(i->s);
			i->All_point[i->M1 + j]->y2 = r * sin(i->s);
		}

		//i->All_point[i->M1 + i->M2 - 2]->x2 = R3 * cos(i->s);   // Устанавливаем точку до контакта 
		//i->All_point[i->M1 + i->M2 - 2]->y2 = R3 * sin(i->s);
		R3 = sqrt(kv(i->Key_point[1]->x2) + kv(i->Key_point[1]->y2));

		//i->All_point[i->M1 + i->M2]->x2 = R3 * cos(i->s);   // Устанавливаем точку после контакта 
		//i->All_point[i->M1 + i->M2]->y2 = R3 * sin(i->s);

		R4 = sqrt(kv(i->Key_point[2]->x2) + kv(i->Key_point[2]->y2));


		for (int j = 0; j < i->M3 - 1; j++)// Передвинули точки до внешней волны
		{
			//x = (j) / (1.0 * (i->M3 - 1));
			//r = R3 + (s_k * x + 3.0 * x * x - 3.0 * s_k * x * x - 2.0 * x * x * x + 2.0 * s_k * x * x * x) * (R4 - R3);
			r = R3 + (R4 - R3) * (j + 1) / (i->M3);         // Равномерное распределение точек
			i->All_point[i->M1 + i->M2 + j]->x2 = r * cos(i->s);
			i->All_point[i->M1 + i->M2 + j]->y2 = r * sin(i->s);
		}

		R4 = sqrt(kv(i->Key_point[2]->x2) + kv(i->Key_point[2]->y2));

		//x = log(R5_ / R4) / (log(1.0 * i->M4 + 1) * (i->M4));
		//double dd = 2.0 * (R5_ - R4) / (i->M4 * (i->M4 + 1.0));
		//r = R4;
		for (int j = 0; j < i->M4; j++)
		{
			//x = (j + 1.0) / (1.0 * i->M4);
			//r = R4 + (kv(x) * kv(x) + (1.0 - x) * 0.03 * x) * (R5_ - R4);
			//r = R4 *  pow( 1.0 * (j + 2), x * (j + 1));  //R4 + (R5_ - R4) * (j + 1) / (i->M4);
			//r = R4 + (R5_ - R4) * (j) / (i->M4 - 1);   // равномерно
			x = (j + 1) * 1.0 / (i->M4);
			double dd = 1.5708 / sqrt(1.0 - 1.0 / 1.76);
			r = R4 + pow(x, 1.76 * (1.0 - kv(i->s / dd))) * (R5_ - R4);
			//r = r + (j + 1) * dd;
			i->All_point[i->M1 + i->M2 + i->M3 + j]->x2 = r * cos(i->s);
			i->All_point[i->M1 + i->M2 + i->M3 + j]->y2 = r * sin(i->s);
		}
	}


	//auto i_ = this->A_Rails[this->A_Rails.size() - 1];
	//double R_distant_shok_1 = sqrt(kv(i_->Key_point[0]->x2) + kv(i_->Key_point[0]->y2));  // n_inner_shock

	//auto i_ = this->B_Rails[this->B_Rails.size() - 1];
	//double R_distant_shok_1 = sqrt(kv(i_->Key_point[0]->x) + kv(i_->Key_point[0]->y));
	//cout << "Setka.cpp    " << R_distant_shok_1 << endl;

	for (int jj = 0; jj < this->B_Rails.size(); jj++)
	{
		auto i = this->B_Rails[jj];
		// Подвинем ключевые точки   0 - внутренняя ударная волна

		double r = 90.0 / RR_;
		i->Key_point[0]->x2 = r * cos(i->s);
		i->Key_point[0]->y2 = r * sin(i->s);

		double xx = i->Key_point[0]->x2;
		double yy = 0.0;
		double y1, y2;
		y1 = 0.0;
		y2 = 3000 / RR_;
		double y0;
		double t = 0.0;
		while (fabs(y1 - y2) > 0.001 / RR_)
		{
			y0 = (y1 + y2) / 2.0;
			t = polar_angle(xx, y0);
			r = 130.0 / RR_ * 0.8 * (cos(t) + 1.5) / (cos(t) + 1.0);
			if (kv(xx) + kv(y0) - kv(r) > 0.0)
			{
				y2 = y0;
			}
			else
			{
				y1 = y0;
			}
		}
		yy = y0;


		i->Key_point[1]->x2 = xx;
		i->Key_point[1]->y2 = yy;

		y1 = 0.0;
		y2 = 3000 / RR_;
		while (fabs(y1 - y2) > 0.001 / RR_)
		{
			y0 = (y1 + y2) / 2.0;
			t = polar_angle(xx, y0);
			r = 245.0 / RR_ * 2.0 * (1.0) / (cos(t) + 1.0);
			if (kv(xx) + kv(y0) - kv(r) > 0.0)
			{
				y2 = y0;
			}
			else
			{
				y1 = y0;
			}
		}
		yy = y0;


		i->Key_point[2]->x2 = xx;
		i->Key_point[2]->y2 = yy;

		i->Key_point[3]->x2 = i->Key_point[0]->x2;
		i->Key_point[3]->y2 = i->Key_point[3]->y;

		R2 = sqrt(kv(i->Key_point[0]->x2) + kv(i->Key_point[0]->y2));

		//double x;
		//x = pow(R11_ / R1_, 1.0 / n_inner);
		//for (int j = 0; j <= n_inner; j++)  // Передвинули точки до ударной волны
		//{
		//	r = R1_ * pow(x, j);
		//	//r = R1_ + (R11_ - R1_) * j / (n_inner);
		//	i->All_point[j]->x2 = r * cos(i->s);
		//	i->All_point[j]->y2 = r * sin(i->s);
		//}


		//for (int j = n_inner + 1; j <= m_inner; j++)  // Передвинули точки до ударной волны
		//{
		//	r = R11_ + (R111_ - R11_) * (j - n_inner) / (m_inner - n_inner);
		//	i->All_point[j]->x2 = r * cos(i->s);
		//	i->All_point[j]->y2 = r * sin(i->s);
		//}

		//for (int j = m_inner + 1; j < i->M1 - 1; j++)  // Передвинули точки до ударной волны
		//{
		//	//r = R11_ * pow(x, j - 7);
		//	r = R111_ + (R2 - R111_) * (j - m_inner) / (i->M1 - 1 - m_inner);
		//	i->All_point[j]->x2 = r * cos(i->s);
		//	i->All_point[j]->y2 = r * sin(i->s);
		//}

		double x;

		for (int j = 0; j < n_inner; j++)
		{
			x = j * 1.0 / (n_inner);
			r = R1_ + pow(x, 1.3) * (R111_ - R1_);
			//r = R1_ + (R11_ - R1_) * j / (n_inner);
			i->All_point[j]->x2 = r * cos(i->s);
			i->All_point[j]->y2 = r * sin(i->s);
		}

		for (int j = n_inner; j < i->M1 - 1; j++)
		{
			x = (j - n_inner) * 1.0 / (i->M1 - 1 - n_inner);
			r = R111_ + pow(x, 1.05) * (R2 - R111_);
			//r = R1_ + (R11_ - R1_) * j / (n_inner);
			i->All_point[j]->x2 = r * cos(i->s);
			i->All_point[j]->y2 = r * sin(i->s);
		}



		R2 = i->Key_point[0]->y2;
		R3 = i->Key_point[1]->y2;

		for (int j = 0; j < i->M2 - 1; j++)
		{
			r = R2 + (R3 - R2) * (j + 1) / (i->M2);
			i->All_point[i->M1 + j]->x2 = i->Key_point[0]->x2;
			i->All_point[i->M1 + j]->y2 = r;
		}

		//i->All_point[i->M1 + i->M2 - 2]->x2 = i->Key_point[0]->x2;   // Устанавливаем точку до контакта 
		//i->All_point[i->M1 + i->M2 - 2]->y2 = R3 ;
		R3 = i->Key_point[1]->y2;

		//i->All_point[i->M1 + i->M2]->x2 = i->Key_point[0]->x2;   // Устанавливаем точку после контакта 
		//i->All_point[i->M1 + i->M2]->y2 = R3 ;

		R4 = i->Key_point[2]->y2;

		for (int j = 0; j < i->M3 - 1; j++)
		{
			//r = R3 + (R4 - R3) * (j) / (i->M3 - 1);
			//x = (j) / (1.0 * (i->M3 - 1));
			//r = R3 + (s_k * x + 3.0 * x * x - 3.0 * s_k * x * x - 2.0 * x * x * x + 2.0 * s_k * x * x * x) * (R4 - R3);
			r = R3 + (R4 - R3) * (j + 1) / (i->M3);
			i->All_point[i->M1 + i->M2 + j]->x2 = i->Key_point[0]->x2;
			i->All_point[i->M1 + i->M2 + j]->y2 = r;
		}

		//x = pow(R5_ / R4, 1.0 / (2.0 * this->M4));
		//double dd = 2.0 * (R5_ - R4) / (i->M4 * (i->M4 + 1.0));
		//r = R4;

		for (int j = 0; j < i->M4; j++)
		{
			//x = (j + 1.0) / (1.0 * i->M4);
			//r = R4 + (kv(x) * kv(x) + (1.0 - x) * (0.03 + 0.19 * jj / (this->B_Rails.size() - 1.0)) * x) * (R5_ - R4);
			//r = R4 * pow(kv(x), j + 1);  //R4 + (R5_ - R4) * (j + 1) / (i->M4);
			r = R4 + (R5_ - R4) * (j + 1) / (i->M4);
			i->All_point[i->M1 + i->M2 + i->M3 + j]->x2 = i->Key_point[0]->x2;
			i->All_point[i->M1 + i->M2 + i->M3 + j]->y2 = r;
		}

	}

	for (auto& i : this->C_Rails)
	{

		double r = 90.0 / RR_;
		i->Key_point[0]->x2 = r * cos(i->s);
		i->Key_point[0]->y2 = r * sin(i->s);

		i->Key_point[1]->x2 = Left_; // i->Key_point[1]->x;
		i->Key_point[1]->y2 = i->Key_point[1]->y;

		R2 = sqrt(kv(i->Key_point[0]->x2) + kv(i->Key_point[0]->y2));

		//double x;
		//x = pow(R11_ / R1_, 1.0 / n_inner);
		//for (int j = 0; j <= n_inner; j++)  // Передвинули точки до ударной волны
		//{
		//	r = R1_ * pow(x, j);
		//	//r = R1_ + (R11_ - R1_) * j / (n_inner);
		//	i->All_point[j]->x2 = r * cos(i->s);
		//	i->All_point[j]->y2 = r * sin(i->s);
		//}


		//for (int j = n_inner + 1; j <= m_inner; j++)  // Передвинули точки до ударной волны
		//{
		//	r = R11_ + (R111_ - R11_) * (j - n_inner) / (m_inner - n_inner);
		//	i->All_point[j]->x2 = r * cos(i->s);
		//	i->All_point[j]->y2 = r * sin(i->s);
		//}

		//for (int j = m_inner + 1; j < i->M1 - 1; j++)  // Передвинули точки до ударной волны
		//{
		//	//r = R11_ * pow(x, j - 7);
		//	r = R111_ + (R2 - R111_) * (j - m_inner) / (i->M1 - 1 - m_inner);
		//	i->All_point[j]->x2 = r * cos(i->s);
		//	i->All_point[j]->y2 = r * sin(i->s);
		//}

		double x;

		for (int j = 0; j < n_inner; j++)
		{
			x = j * 1.0 / (n_inner);
			r = R1_ + pow(x, 1.3) * (R111_ - R1_);
			//r = R1_ + (R11_ - R1_) * j / (n_inner);
			i->All_point[j]->x2 = r * cos(i->s);
			i->All_point[j]->y2 = r * sin(i->s);
		}

		for (int j = n_inner; j < i->M1 - 1; j++)
		{
			x = (j - n_inner) * 1.0 / (i->M1 - 1 - n_inner);
			r = R111_ + pow(x, 1.05) * (R2 - R111_);
			//r = R1_ + (R11_ - R1_) * j / (n_inner);
			i->All_point[j]->x2 = r * cos(i->s);
			i->All_point[j]->y2 = r * sin(i->s);
		}

		R2 = i->Key_point[0]->x2;

		R3 = i->Key_point[1]->x2;

		x = pow(R3 / (R2), 1.0 / (i->M2));

		for (int j = i->M1; j < i->All_point.size(); j++)
		{
			r = R2 * pow(x, j + 1 - i->M1);
			i->All_point[j]->x2 = r;
			i->All_point[j]->y2 = i->Key_point[0]->y2;
		}
	}

	for (int jj = 0; jj < this->D_Rails.size(); jj++)
	{
		double x;
		// Подвинем ключевые точки

		auto i = this->D_Rails[jj];

		double xx = i->Key_point[0]->x2;
		double yy = 0.0;
		double y1, y2;
		y1 = 0.0;
		y2 = 3000 / RR_;
		double y0;
		double t = 0.0;
		while (fabs(y1 - y2) > 0.001 / RR_)
		{
			y0 = (y1 + y2) / 2.0;
			t = polar_angle(xx, y0);
			r = 130.0 / RR_ * 0.8 * (cos(t) + 1.5) / (cos(t) + 1.0);
			if (kv(xx) + kv(y0) - kv(r) > 0.0)
			{
				y2 = y0;
			}
			else
			{
				y1 = y0;
			}
		}
		yy = y0;


		i->Key_point[1]->x2 = xx;
		i->Key_point[1]->y2 = yy;


		y1 = 0.0;
		y2 = 3000 / RR_;
		while (fabs(y1 - y2) > 0.001 / RR_)
		{
			y0 = (y1 + y2) / 2.0;
			t = polar_angle(xx, y0);
			r = 245.0 / RR_ * 2.0 * (1.0) / (cos(t) + 1.0);
			if (kv(xx) + kv(y0) - kv(r) > 0.0)
			{
				y2 = y0;
			}
			else
			{
				y1 = y0;
			}
		}
		yy = y0;


		i->Key_point[2]->x2 = xx;
		i->Key_point[2]->y2 = yy;


		i->Key_point[3]->x2 = i->Key_point[0]->x2;
		i->Key_point[3]->y2 = i->Key_point[3]->y;

		R2 = i->Key_point[0]->y2;

		R3 = i->Key_point[1]->y2;

		for (int j = 0; j < this->M2 - 1; j++)
		{
			r = R2 + (R3 - R2) * (j + 1) / (i->M2);
			i->All_point[j + 1]->x2 = i->Key_point[0]->x2;
			i->All_point[j + 1]->y2 = r;
		}

		//i->All_point[1 + i->M2 - 2]->x2 = i->Key_point[0]->x2;   // Устанавливаем точку до контакта 
		//i->All_point[1 + i->M2 - 2]->y2 = R3;
		R3 = i->Key_point[1]->y2;

		//i->All_point[1 + i->M2]->x2 = i->Key_point[0]->x2;   // Устанавливаем точку после контакта 
		//i->All_point[1 + i->M2]->y2 = R3;

		R4 = i->Key_point[2]->y2;

		for (int j = 0; j < this->M3 - 1; j++)
		{
			//r = R3 + (R4 - R3) * (j) / (i->M3 - 1);
			//x = (j) / (1.0 * (i->M3 - 1));
			//r = R3 + (s_k * x + 3.0 * x * x - 3.0 * s_k * x * x - 2.0 * x * x * x + 2.0 * s_k * x * x * x) * (R4 - R3);
			r = R3 + (R4 - R3) * (j + 1) / (i->M3);
			i->All_point[this->M2 + j + 1]->x2 = i->Key_point[0]->x2;
			i->All_point[this->M2 + j + 1]->y2 = r;
		}

		//double x = pow(R5_ / R4, 1.0 / (2.0 * this->M4));
		//double dd = 2.0 * (R5_ - R4) / (i->M4 * (i->M4 + 1.0));
		//r = R4;
		for (int j = 0; j < this->M4; j++)
		{
			//double x = (j + 1.0) / (1.0 * i->M4);
			//r = R4 + (kv(x) * kv(x) + (1.0 - x) * 0.22 * x) * (R5_ - R4);
			//r = R4 * pow(kv(x), j + 1);  //R4 + (R5_ - R4) * (j + 1) / (i->M4);
			r = R4 + (R5_ - R4) * (j + 1) / (i->M4);
			i->All_point[this->M2 + this->M3 + j + 1]->x2 = i->Key_point[0]->x2;
			i->All_point[this->M2 + this->M3 + j + 1]->y2 = r;
		}


	}
}

void Setka::Move_Setka_Calculate_2(const double& dt)
{
	double Vx, Vy, V;
	double R2, r, R3, R4;

	//double y_Outer = this->D_Rails[16]->Key_point[1]->y;

	for (int jj = 0; jj < this->A_Rails.size(); jj++)
	{
		auto i = this->A_Rails[jj];
		// Подвинем ключевые точки
		V = i->Key_point[0]->Vx * cos(i->s) + i->Key_point[0]->Vy * sin(i->s);
		i->Key_point[0]->x2 = i->Key_point[0]->x + dt * V * cos(i->s);
		i->Key_point[0]->y2 = i->Key_point[0]->y + dt * V * sin(i->s);
		//i->Key_point[0]->x2 = i->Key_point[0]->x/1.3;
		//i->Key_point[0]->y2 = i->Key_point[0]->y/1.3;

		//cout << "Setka.cpp    " << jj << " " << i->Key_point[1]->Vx << " " << i->Key_point[1]->Vy << endl;
		V = i->Key_point[1]->Vx * cos(i->s) + i->Key_point[1]->Vy * sin(i->s);
		i->Key_point[1]->x2 = i->Key_point[1]->x + dt * V * cos(i->s);
		i->Key_point[1]->y2 = i->Key_point[1]->y + dt * V * sin(i->s);
		//i->Key_point[1]->x2 = i->Key_point[1]->x/1.3;
		//i->Key_point[1]->y2 = i->Key_point[1]->y/1.3;

		V = i->Key_point[2]->Vx * cos(i->s) + i->Key_point[2]->Vy * sin(i->s);
		i->Key_point[2]->x2 = i->Key_point[2]->x + dt * V * cos(i->s);
		i->Key_point[2]->y2 = i->Key_point[2]->y + dt * V * sin(i->s);
		//i->Key_point[2]->x2 = i->Key_point[2]->x / 1.3;
		//i->Key_point[2]->y2 = i->Key_point[2]->y / 1.3;

		i->Key_point[3]->x2 = i->Key_point[3]->x;
		i->Key_point[3]->y2 = i->Key_point[3]->y;


		R2 = sqrt(kv(i->Key_point[0]->x2) + kv(i->Key_point[0]->y2));

		//double x;

		//x = pow(R11_ / R1_, 1.0 / n_inner);
		//for (int j = 0; j <= n_inner; j++) 
		//{
		//	r = R1_ * pow(x, j);
		//	//r = R1_ + (R11_ - R1_) * j / (n_inner);
		//	i->All_point[j]->x2 = r * cos(i->s);
		//	i->All_point[j]->y2 = r * sin(i->s);
		//}


		//for (int j = n_inner + 1; j <= m_inner; j++) 
		//{
		//	r = R11_ + (R111_ - R11_) * (j - n_inner) / (m_inner - n_inner);
		//	i->All_point[j]->x2 = r * cos(i->s);
		//	i->All_point[j]->y2 = r * sin(i->s);
		//}

		//for (int j = m_inner + 1; j < i->M1 - 1; j++)  // Передвинули точки до ударной волны
		//{
		//	//r = R11_ * pow(x, j - 7);
		//	r = R111_ + (R2 - R111_) * (j - m_inner) / (i->M1 - 1 - m_inner);
		//	i->All_point[j]->x2 = r * cos(i->s);
		//	i->All_point[j]->y2 = r * sin(i->s);
		//}

		double x;

		//x = pow(R2 / R1_, 1.0 / (i->M1 - 1));
		for (int j = 0; j < n_inner; j++)
		{
			x = j * 1.0 / (n_inner);
			r = R1_ +  pow(x, 1.3) * (R111_ - R1_);
			//r = R1_ + (R11_ - R1_) * j / (n_inner);
			i->All_point[j]->x2 = r * cos(i->s);
			i->All_point[j]->y2 = r * sin(i->s);
		}

		for (int j = n_inner; j < i->M1 - 1; j++)
		{
			x = (j - n_inner) * 1.0 / (i->M1 - 1 - n_inner);
			r = R111_ + pow(x, 1.05) * (R2 - R111_);
			//r = R1_ + (R11_ - R1_) * j / (n_inner);
			i->All_point[j]->x2 = r * cos(i->s);
			i->All_point[j]->y2 = r * sin(i->s);
		}



		//R3 = sqrt(kv(i->Key_point[1]->x2) + kv(i->Key_point[1]->y2)) - zazor;
		R3 = sqrt(kv(i->Key_point[1]->x2) + kv(i->Key_point[1]->y2));

		for (int j = 0; j < i->M2 - 1; j++)  // Передвинули точки до контакта
		{
			r = R2 + (R3 - R2) * (j + 1) / (i->M2);         // Равномерное распределение точек
			i->All_point[i->M1 + j]->x2 = r * cos(i->s);
			i->All_point[i->M1 + j]->y2 = r * sin(i->s);
		}

		//i->All_point[i->M1 + i->M2 - 2]->x2 = R3 * cos(i->s);   // Устанавливаем точку до контакта 
		//i->All_point[i->M1 + i->M2 - 2]->y2 = R3 * sin(i->s);
		R3 = sqrt(kv(i->Key_point[1]->x2) + kv(i->Key_point[1]->y2));

		//i->All_point[i->M1 + i->M2]->x2 = R3 * cos(i->s);   // Устанавливаем точку после контакта 
		//i->All_point[i->M1 + i->M2]->y2 = R3 * sin(i->s);

		R4 = sqrt(kv(i->Key_point[2]->x2) + kv(i->Key_point[2]->y2));


		for (int j = 0; j < i->M3 - 1; j++)// Передвинули точки до внешней волны
		{
			//x = (j) / (1.0 * (i->M3 - 1));
			//r = R3 + (s_k * x + 3.0 * x * x - 3.0 * s_k * x * x - 2.0 * x * x * x + 2.0 * s_k * x * x * x) * (R4 - R3);
			r = R3 + (R4 - R3) * (j + 1) / (i->M3);         // Равномерное распределение точек
			i->All_point[i->M1 + i->M2 + j]->x2 = r * cos(i->s);
			i->All_point[i->M1 + i->M2 + j]->y2 = r * sin(i->s);
		}

		R4 = sqrt(kv(i->Key_point[2]->x2) + kv(i->Key_point[2]->y2));

		//x = log(R5_ / R4) / (log(1.0 * i->M4 + 1) * (i->M4));
		//double dd = 2.0 * (R5_ - R4) / (i->M4 * (i->M4 + 1.0));
		//r = R4;
		for (int j = 0; j < i->M4; j++)
		{
			//x = (j + 1.0) / (1.0 * i->M4);
			//r = R4 + (kv(x) * kv(x) + (1.0 - x) * 0.03 * x) * (R5_ - R4);
			//r = R4 *  pow( 1.0 * (j + 2), x * (j + 1));  //R4 + (R5_ - R4) * (j + 1) / (i->M4);
			//r = R4 + (R5_ - R4) * (j) / (i->M4 - 1);   // равномерно
			x = (j + 1) * 1.0 / (i->M4);
			double dd = 1.5708/sqrt(1.0 - 1.0/1.76);
			r = R4 + pow(x, 1.76 * (1.0 - kv(i->s/dd))) * (R5_ - R4);
			//r = r + (j + 1) * dd;
			i->All_point[i->M1 + i->M2 + i->M3 + j]->x2 = r * cos(i->s);
			i->All_point[i->M1 + i->M2 + i->M3 + j]->y2 = r * sin(i->s);
		}
	}


	//auto i_ = this->A_Rails[this->A_Rails.size() - 1];
	//double R_distant_shok_1 = sqrt(kv(i_->Key_point[0]->x2) + kv(i_->Key_point[0]->y2));  // n_inner_shock

	//auto i_ = this->B_Rails[this->B_Rails.size() - 1];
	//double R_distant_shok_1 = sqrt(kv(i_->Key_point[0]->x) + kv(i_->Key_point[0]->y));
	//cout << "Setka.cpp    " << R_distant_shok_1 << endl;

	for (int jj = 0; jj < this->B_Rails.size(); jj++)
	{
		auto i = this->B_Rails[jj];
		// Подвинем ключевые точки   0 - внутренняя ударная волна
		if (true)//(jj <= n_inner_shock)
		{
			V = i->Key_point[0]->Vx * cos(i->s) + i->Key_point[0]->Vy * sin(i->s);
			i->Key_point[0]->x2 = i->Key_point[0]->x + dt * V * cos(i->s);
			i->Key_point[0]->y2 = i->Key_point[0]->y + dt * V * sin(i->s);
			//i->Key_point[0]->x2 = i->Key_point[0]->x / 1.3;
			//i->Key_point[0]->y2 = i->Key_point[0]->y / 1.3;
		}
		else
		{
			//V = i->Key_point[0]->Vx * cos(i->s) + i->Key_point[0]->Vy * sin(i->s);
			//i->Key_point[0]->x2 = R_distant_shok_1 * cos(i->s);
			//i->Key_point[0]->y2 = R_distant_shok_1 * sin(i->s);
		}

		i->Key_point[0]->Vx = i->Key_point[0]->Vy = 0.0;

		// Особое движение для стыковой точки
		/*if (jj == this->B_Rails.size() - 1)
		{
			r = sqrt(kv(this->B_Rails[jj - 1]->Key_point[0]->x2) + kv(this->B_Rails[jj - 1]->Key_point[0]->y2));
			i->Key_point[0]->x2 = r * cos(i->s);
			i->Key_point[0]->y2 = r * sin(i->s);
		}*/

		V = i->Key_point[1]->Vx * 0.0 + i->Key_point[1]->Vy * 1.0;
		i->Key_point[1]->x2 = i->Key_point[0]->x2;
		i->Key_point[1]->y2 = i->Key_point[1]->y + dt * V;
		//i->Key_point[1]->x2 = i->Key_point[0]->x / 1.3;
		//i->Key_point[1]->y2 = i->Key_point[1]->y / 1.3;


		V = i->Key_point[2]->Vx * 0.0 + i->Key_point[2]->Vy * 1.0;
		i->Key_point[2]->x2 = i->Key_point[0]->x2;
		i->Key_point[2]->y2 = i->Key_point[2]->y + dt * V;
		//i->Key_point[2]->x2 = i->Key_point[0]->x / 1.3;
		//i->Key_point[2]->y2 = i->Key_point[2]->y / 1.3;

		i->Key_point[3]->x2 = i->Key_point[0]->x2;
		i->Key_point[3]->y2 = i->Key_point[3]->y;

		R2 = sqrt(kv(i->Key_point[0]->x2) + kv(i->Key_point[0]->y2));

		//double x;
		//x = pow(R11_ / R1_, 1.0 / n_inner);
		//for (int j = 0; j <= n_inner; j++)  // Передвинули точки до ударной волны
		//{
		//	r = R1_ * pow(x, j);
		//	//r = R1_ + (R11_ - R1_) * j / (n_inner);
		//	i->All_point[j]->x2 = r * cos(i->s);
		//	i->All_point[j]->y2 = r * sin(i->s);
		//}


		//for (int j = n_inner + 1; j <= m_inner; j++)  // Передвинули точки до ударной волны
		//{
		//	r = R11_ + (R111_ - R11_) * (j - n_inner) / (m_inner - n_inner);
		//	i->All_point[j]->x2 = r * cos(i->s);
		//	i->All_point[j]->y2 = r * sin(i->s);
		//}

		//for (int j = m_inner + 1; j < i->M1 - 1; j++)  // Передвинули точки до ударной волны
		//{
		//	//r = R11_ * pow(x, j - 7);
		//	r = R111_ + (R2 - R111_) * (j - m_inner) / (i->M1 - 1 - m_inner);
		//	i->All_point[j]->x2 = r * cos(i->s);
		//	i->All_point[j]->y2 = r * sin(i->s);
		//}

		double x;

		for (int j = 0; j < n_inner; j++)
		{
			x = j * 1.0 / (n_inner);
			r = R1_ + pow(x, 1.3) * (R111_ - R1_);
			//r = R1_ + (R11_ - R1_) * j / (n_inner);
			i->All_point[j]->x2 = r * cos(i->s);
			i->All_point[j]->y2 = r * sin(i->s);
		}

		for (int j = n_inner; j < i->M1 - 1; j++)
		{
			x = (j - n_inner) * 1.0 / (i->M1 - 1 - n_inner);
			r = R111_ + pow(x, 1.05) * (R2 - R111_);
			//r = R1_ + (R11_ - R1_) * j / (n_inner);
			i->All_point[j]->x2 = r * cos(i->s);
			i->All_point[j]->y2 = r * sin(i->s);
		}

		
		
		R2 = i->Key_point[0]->y2;
		R3 = i->Key_point[1]->y2;

		for (int j = 0; j < i->M2 - 1; j++)
		{
			r = R2 + (R3 - R2) * (j + 1) / (i->M2);
			i->All_point[i->M1 + j]->x2 = i->Key_point[0]->x2;
			i->All_point[i->M1 + j]->y2 = r;
		}

		//i->All_point[i->M1 + i->M2 - 2]->x2 = i->Key_point[0]->x2;   // Устанавливаем точку до контакта 
		//i->All_point[i->M1 + i->M2 - 2]->y2 = R3 ;
		R3 = i->Key_point[1]->y2;

		//i->All_point[i->M1 + i->M2]->x2 = i->Key_point[0]->x2;   // Устанавливаем точку после контакта 
		//i->All_point[i->M1 + i->M2]->y2 = R3 ;

		R4 = i->Key_point[2]->y2;

		for (int j = 0; j < i->M3 - 1; j++)
		{
			//r = R3 + (R4 - R3) * (j) / (i->M3 - 1);
			//x = (j) / (1.0 * (i->M3 - 1));
			//r = R3 + (s_k * x + 3.0 * x * x - 3.0 * s_k * x * x - 2.0 * x * x * x + 2.0 * s_k * x * x * x) * (R4 - R3);
			r = R3 + (R4 - R3) * (j + 1) / (i->M3);
			i->All_point[i->M1 + i->M2 + j]->x2 = i->Key_point[0]->x2;
			i->All_point[i->M1 + i->M2 + j]->y2 = r;
		}

		//x = pow(R5_ / R4, 1.0 / (2.0 * this->M4));
		//double dd = 2.0 * (R5_ - R4) / (i->M4 * (i->M4 + 1.0));
		//r = R4;

		for (int j = 0; j < i->M4; j++)
		{
			//x = (j + 1.0) / (1.0 * i->M4);
			//r = R4 + (kv(x) * kv(x) + (1.0 - x) * (0.03 + 0.19 * jj / (this->B_Rails.size() - 1.0)) * x) * (R5_ - R4);
			//r = R4 * pow(kv(x), j + 1);  //R4 + (R5_ - R4) * (j + 1) / (i->M4);
			r = R4 + (R5_ - R4) * (j + 1) / (i->M4);
			i->All_point[i->M1 + i->M2 + i->M3 + j]->x2 = i->Key_point[0]->x2;
			i->All_point[i->M1 + i->M2 + i->M3 + j]->y2 = r;
		}

	}

	for (int jj = 0; jj < this->C_Rails.size(); jj++)
	{
		auto i = this->C_Rails[jj];
		/*V = i->Key_point[0]->Vx * cos(i->s) + i->Key_point[0]->Vy * sin(i->s);
		i->Key_point[0]->x2 = R_distant_shok_1 * cos(i->s);
		i->Key_point[0]->y2 = R_distant_shok_1 * sin(i->s);*/

		V = i->Key_point[0]->Vx * cos(i->s) + i->Key_point[0]->Vy * sin(i->s);

		if (jj != 0) // Этот узел уже был подвинут в B - лучах
		{
			i->Key_point[0]->x2 = i->Key_point[0]->x + dt * V * cos(i->s);
			i->Key_point[0]->y2 = i->Key_point[0]->y + dt * V * sin(i->s);

			/*if (jj == 1)
			{
				r = sqrt(kv(this->C_Rails[jj - 1]->Key_point[0]->x2) + kv(this->C_Rails[jj - 1]->Key_point[0]->y2));
				i->Key_point[0]->x2 = r * cos(i->s);
				i->Key_point[0]->y2 = r * sin(i->s);
			}*/
		}
		//i->Key_point[0]->x2 = i->Key_point[0]->x / 1.3;
		//i->Key_point[0]->y2 = i->Key_point[0]->y / 1.3;

		i->Key_point[1]->x2 = Left_; // i->Key_point[1]->x;
		i->Key_point[1]->y2 = i->Key_point[1]->y;
		//i->Key_point[1]->x2 = Left_;
		//i->Key_point[1]->y2 = i->Key_point[1]->y / 1.3;

		R2 = sqrt(kv(i->Key_point[0]->x2) + kv(i->Key_point[0]->y2));

		//double x;
		//x = pow(R11_ / R1_, 1.0 / n_inner);
		//for (int j = 0; j <= n_inner; j++)  // Передвинули точки до ударной волны
		//{
		//	r = R1_ * pow(x, j);
		//	//r = R1_ + (R11_ - R1_) * j / (n_inner);
		//	i->All_point[j]->x2 = r * cos(i->s);
		//	i->All_point[j]->y2 = r * sin(i->s);
		//}


		//for (int j = n_inner + 1; j <= m_inner; j++)  // Передвинули точки до ударной волны
		//{
		//	r = R11_ + (R111_ - R11_) * (j - n_inner) / (m_inner - n_inner);
		//	i->All_point[j]->x2 = r * cos(i->s);
		//	i->All_point[j]->y2 = r * sin(i->s);
		//}

		//for (int j = m_inner + 1; j < i->M1 - 1; j++)  // Передвинули точки до ударной волны
		//{
		//	//r = R11_ * pow(x, j - 7);
		//	r = R111_ + (R2 - R111_) * (j - m_inner) / (i->M1 - 1 - m_inner);
		//	i->All_point[j]->x2 = r * cos(i->s);
		//	i->All_point[j]->y2 = r * sin(i->s);
		//}

		double x;

		for (int j = 0; j < n_inner; j++)
		{
			x = j * 1.0 / (n_inner);
			r = R1_ + pow(x, 1.3) * (R111_ - R1_);
			//r = R1_ + (R11_ - R1_) * j / (n_inner);
			i->All_point[j]->x2 = r * cos(i->s);
			i->All_point[j]->y2 = r * sin(i->s);
		}

		for (int j = n_inner; j < i->M1 - 1; j++)
		{
			x = (j - n_inner) * 1.0 / (i->M1 - 1 - n_inner);
			r = R111_ + pow(x, 1.05) * (R2 - R111_);
			//r = R1_ + (R11_ - R1_) * j / (n_inner);
			i->All_point[j]->x2 = r * cos(i->s);
			i->All_point[j]->y2 = r * sin(i->s);
		}

		R2 = i->Key_point[0]->x2;

		R3 = i->Key_point[1]->x2;

		x = pow(R3 / (R2), 1.0 / (i->M2));

		for (int j = i->M1; j < i->All_point.size(); j++)
		{
			r = R2  * pow(x, j + 1 - i->M1);
			i->All_point[j]->x2 = r;
			i->All_point[j]->y2 = i->Key_point[0]->y2;
		}
	}

	for (int jj = 0; jj < this->D_Rails.size(); jj++)
	{
		double x;
		// Подвинем ключевые точки
		auto i = this->D_Rails[jj];
		V = i->Key_point[1]->Vy;
		i->Key_point[1]->x2 = i->Key_point[0]->x2;
		/*if (jj > 16)
		{
			i->Key_point[1]->y2 = y_Outer + (247.0 - y_Outer) * (jj - 16.0) / (this->D_Rails.size() - 1.0 - 16.0);
		}
		else
		{
			i->Key_point[1]->y2 = i->Key_point[1]->y + dt * V;
		}*/
		i->Key_point[1]->y2 = i->Key_point[1]->y + dt * V;
		//i->Key_point[1]->y2 = i->Key_point[1]->y/1.3;
		/*if (jj > 22)
		{
			i->Key_point[1]->y2 = this->D_Rails[22]->Key_point[1]->y2;
		}*/

		V = i->Key_point[2]->Vy;
		i->Key_point[2]->x2 = i->Key_point[0]->x2;
		i->Key_point[2]->y2 = i->Key_point[2]->y + dt * V;
		//i->Key_point[2]->y2 = i->Key_point[2]->y/1.3;

		if (jj >= n_outer_shock)                                                        // ТУТ ИЗМЕНИЛ
		{
			i->Key_point[2]->y2 = this->D_Rails[jj - 1]->Key_point[2]->y2;
			//cout << i->Key_point[2]->x << " " << i->Key_point[2]->y << endl;
			//cout << "Setka.cpp    " << this->D_Rails[jj - 1]->Key_point[2]->y2 << endl;
			//cout << "Setka.cpp    " << this->D_Rails[jj - 2]->Key_point[2]->y2 << endl;
			//exit(-1);
		}

		i->Key_point[3]->x2 = i->Key_point[0]->x2;
		i->Key_point[3]->y2 = i->Key_point[3]->y;

		R2 = i->Key_point[0]->y2;

		R3 = i->Key_point[1]->y2;

		for (int j = 0; j < this->M2 - 1; j++)
		{
			r = R2 + (R3 - R2) * (j + 1) / (i->M2);
			i->All_point[j + 1]->x2 = i->Key_point[0]->x2;
			i->All_point[j + 1]->y2 = r;
		}

		//i->All_point[1 + i->M2 - 2]->x2 = i->Key_point[0]->x2;   // Устанавливаем точку до контакта 
		//i->All_point[1 + i->M2 - 2]->y2 = R3;
		R3 = i->Key_point[1]->y2;

		//i->All_point[1 + i->M2]->x2 = i->Key_point[0]->x2;   // Устанавливаем точку после контакта 
		//i->All_point[1 + i->M2]->y2 = R3;

		R4 = i->Key_point[2]->y2;

		for (int j = 0; j < this->M3 - 1; j++)
		{
			//r = R3 + (R4 - R3) * (j) / (i->M3 - 1);
			//x = (j) / (1.0 * (i->M3 - 1));
			//r = R3 + (s_k * x + 3.0 * x * x - 3.0 * s_k * x * x - 2.0 * x * x * x + 2.0 * s_k * x * x * x) * (R4 - R3);
			r = R3 + (R4 - R3) * (j + 1) / (i->M3);
			i->All_point[this->M2 + j + 1]->x2 = i->Key_point[0]->x2;
			i->All_point[this->M2 + j + 1]->y2 = r;
		}

		//double x = pow(R5_ / R4, 1.0 / (2.0 * this->M4));
		//double dd = 2.0 * (R5_ - R4) / (i->M4 * (i->M4 + 1.0));
		//r = R4;
		for (int j = 0; j < this->M4; j++)
		{
			//double x = (j + 1.0) / (1.0 * i->M4);
			//r = R4 + (kv(x) * kv(x) + (1.0 - x) * 0.22 * x) * (R5_ - R4);
			//r = R4 * pow(kv(x), j + 1);  //R4 + (R5_ - R4) * (j + 1) / (i->M4);
			r = R4 + (R5_ - R4) * (j + 1) / (i->M4);
			i->All_point[this->M2 + this->M3 + j + 1]->x2 = i->Key_point[0]->x2;
			i->All_point[this->M2 + this->M3 + j + 1]->y2 = r;
		}


	}
}

void Setka::Move_Setka_Calculate_stat()
{
	double Vx, Vy, V;
	double R2, r, R3, R4;

	//double y_Outer = this->D_Rails[16]->Key_point[1]->y;

	for (int jj = 0; jj < this->A_Rails.size(); jj++)
	{
		auto i = this->A_Rails[jj];
		


		R2 = sqrt(kv(i->Key_point[0]->x2) + kv(i->Key_point[0]->y2));

		double x;

		for (int j = 0; j < i->M1 - 1; j++)
		{
			x = (j) * 1.0 / (i->M1 - 1);
			r = R1_ + pow(x, 1.05) * (R2 - R1_);
			//r = R1_ + (R11_ - R1_) * j / (n_inner);
			i->All_point[j]->x2 = r * cos(i->s);
			i->All_point[j]->y2 = r * sin(i->s);
		}



		//R3 = sqrt(kv(i->Key_point[1]->x2) + kv(i->Key_point[1]->y2)) - zazor;
		R3 = sqrt(kv(i->Key_point[1]->x2) + kv(i->Key_point[1]->y2));

		for (int j = 0; j < i->M2 - 1; j++)  // Передвинули точки до контакта
		{
			r = R2 + (R3 - R2) * (j + 1) / (i->M2);         // Равномерное распределение точек
			i->All_point[i->M1 + j]->x2 = r * cos(i->s);
			i->All_point[i->M1 + j]->y2 = r * sin(i->s);
		}

		//i->All_point[i->M1 + i->M2 - 2]->x2 = R3 * cos(i->s);   // Устанавливаем точку до контакта 
		//i->All_point[i->M1 + i->M2 - 2]->y2 = R3 * sin(i->s);
		R3 = sqrt(kv(i->Key_point[1]->x2) + kv(i->Key_point[1]->y2));

		//i->All_point[i->M1 + i->M2]->x2 = R3 * cos(i->s);   // Устанавливаем точку после контакта 
		//i->All_point[i->M1 + i->M2]->y2 = R3 * sin(i->s);

		R4 = sqrt(kv(i->Key_point[2]->x2) + kv(i->Key_point[2]->y2));


		for (int j = 0; j < i->M3 - 1; j++)// Передвинули точки до внешней волны
		{
			//x = (j) / (1.0 * (i->M3 - 1));
			//r = R3 + (s_k * x + 3.0 * x * x - 3.0 * s_k * x * x - 2.0 * x * x * x + 2.0 * s_k * x * x * x) * (R4 - R3);
			r = R3 + (R4 - R3) * (j + 1) / (i->M3);         // Равномерное распределение точек
			i->All_point[i->M1 + i->M2 + j]->x2 = r * cos(i->s);
			i->All_point[i->M1 + i->M2 + j]->y2 = r * sin(i->s);
		}

		R4 = sqrt(kv(i->Key_point[2]->x2) + kv(i->Key_point[2]->y2));

		//x = log(R5_ / R4) / (log(1.0 * i->M4 + 1) * (i->M4));
		//double dd = 2.0 * (R5_ - R4) / (i->M4 * (i->M4 + 1.0));
		//r = R4;
		for (int j = 0; j < i->M4; j++)
		{
			//x = (j + 1.0) / (1.0 * i->M4);
			//r = R4 + (kv(x) * kv(x) + (1.0 - x) * 0.03 * x) * (R5_ - R4);
			//r = R4 *  pow( 1.0 * (j + 2), x * (j + 1));  //R4 + (R5_ - R4) * (j + 1) / (i->M4);
			//r = R4 + (R5_ - R4) * (j) / (i->M4 - 1);   // равномерно
			x = (j + 1) * 1.0 / (i->M4);
			double dd = 1.5708 / sqrt(1.0 - 1.0 / 1.76);
			r = R4 + pow(x, 1.76 * (1.0 - kv(i->s / dd))) * (R5_ - R4);
			//r = r + (j + 1) * dd;
			i->All_point[i->M1 + i->M2 + i->M3 + j]->x2 = r * cos(i->s);
			i->All_point[i->M1 + i->M2 + i->M3 + j]->y2 = r * sin(i->s);
		}
	}

	for (int jj = 0; jj < this->B_Rails.size(); jj++)
	{
		auto i = this->B_Rails[jj];
		// Подвинем ключевые точки   0 - внутренняя ударная волна

		R2 = sqrt(kv(i->Key_point[0]->x2) + kv(i->Key_point[0]->y2));


		double x;

		for (int j = 0; j < i->M1 - 1; j++)
		{
			x = (j) * 1.0 / (i->M1 - 1);
			r = R1_ + pow(x, 1.05) * (R2 - R1_);
			//r = R1_ + (R11_ - R1_) * j / (n_inner);
			i->All_point[j]->x2 = r * cos(i->s);
			i->All_point[j]->y2 = r * sin(i->s);
		}



		R2 = i->Key_point[0]->y2;
		R3 = i->Key_point[1]->y2;

		for (int j = 0; j < i->M2 - 1; j++)
		{
			r = R2 + (R3 - R2) * (j + 1) / (i->M2);
			i->All_point[i->M1 + j]->x2 = i->Key_point[0]->x2;
			i->All_point[i->M1 + j]->y2 = r;
		}

		//i->All_point[i->M1 + i->M2 - 2]->x2 = i->Key_point[0]->x2;   // Устанавливаем точку до контакта 
		//i->All_point[i->M1 + i->M2 - 2]->y2 = R3 ;
		R3 = i->Key_point[1]->y2;

		//i->All_point[i->M1 + i->M2]->x2 = i->Key_point[0]->x2;   // Устанавливаем точку после контакта 
		//i->All_point[i->M1 + i->M2]->y2 = R3 ;

		R4 = i->Key_point[2]->y2;

		for (int j = 0; j < i->M3 - 1; j++)
		{
			//r = R3 + (R4 - R3) * (j) / (i->M3 - 1);
			//x = (j) / (1.0 * (i->M3 - 1));
			//r = R3 + (s_k * x + 3.0 * x * x - 3.0 * s_k * x * x - 2.0 * x * x * x + 2.0 * s_k * x * x * x) * (R4 - R3);
			r = R3 + (R4 - R3) * (j + 1) / (i->M3);
			i->All_point[i->M1 + i->M2 + j]->x2 = i->Key_point[0]->x2;
			i->All_point[i->M1 + i->M2 + j]->y2 = r;
		}

		//x = pow(R5_ / R4, 1.0 / (2.0 * this->M4));
		//double dd = 2.0 * (R5_ - R4) / (i->M4 * (i->M4 + 1.0));
		//r = R4;

		for (int j = 0; j < i->M4; j++)
		{
			//x = (j + 1.0) / (1.0 * i->M4);
			//r = R4 + (kv(x) * kv(x) + (1.0 - x) * (0.03 + 0.19 * jj / (this->B_Rails.size() - 1.0)) * x) * (R5_ - R4);
			//r = R4 * pow(kv(x), j + 1);  //R4 + (R5_ - R4) * (j + 1) / (i->M4);
			r = R4 + (R5_ - R4) * (j + 1) / (i->M4);
			i->All_point[i->M1 + i->M2 + i->M3 + j]->x2 = i->Key_point[0]->x2;
			i->All_point[i->M1 + i->M2 + i->M3 + j]->y2 = r;
		}

	}

	for (int jj = 0; jj < this->C_Rails.size(); jj++)
	{
		auto i = this->C_Rails[jj];

		R2 = sqrt(kv(i->Key_point[0]->x2) + kv(i->Key_point[0]->y2));

		double x;

		for (int j = 0; j < i->M1 - 1; j++)
		{
			x = (j) * 1.0 / (i->M1 - 1);
			r = R1_ + pow(x, 1.05) * (R2 - R1_);
			//r = R1_ + (R11_ - R1_) * j / (n_inner);
			i->All_point[j]->x2 = r * cos(i->s);
			i->All_point[j]->y2 = r * sin(i->s);
		}

		R2 = i->Key_point[0]->x2;

		R3 = i->Key_point[1]->x2;

		x = pow(R3 / (R2), 1.0 / (i->M2));

		for (int j = i->M1; j < i->All_point.size(); j++)
		{
			r = R2 * pow(x, j + 1 - i->M1);
			i->All_point[j]->x2 = r;
			i->All_point[j]->y2 = i->Key_point[0]->y2;
		}
	}

	for (int jj = 0; jj < this->D_Rails.size(); jj++)
	{
		double x;
		// Подвинем ключевые точки
		auto i = this->D_Rails[jj];

		R2 = i->Key_point[0]->y2;

		R3 = i->Key_point[1]->y2;

		for (int j = 0; j < this->M2 - 1; j++)
		{
			r = R2 + (R3 - R2) * (j + 1) / (i->M2);
			i->All_point[j + 1]->x2 = i->Key_point[0]->x2;
			i->All_point[j + 1]->y2 = r;
		}

		//i->All_point[1 + i->M2 - 2]->x2 = i->Key_point[0]->x2;   // Устанавливаем точку до контакта 
		//i->All_point[1 + i->M2 - 2]->y2 = R3;
		R3 = i->Key_point[1]->y2;

		//i->All_point[1 + i->M2]->x2 = i->Key_point[0]->x2;   // Устанавливаем точку после контакта 
		//i->All_point[1 + i->M2]->y2 = R3;

		R4 = i->Key_point[2]->y2;

		for (int j = 0; j < this->M3 - 1; j++)
		{
			//r = R3 + (R4 - R3) * (j) / (i->M3 - 1);
			//x = (j) / (1.0 * (i->M3 - 1));
			//r = R3 + (s_k * x + 3.0 * x * x - 3.0 * s_k * x * x - 2.0 * x * x * x + 2.0 * s_k * x * x * x) * (R4 - R3);
			r = R3 + (R4 - R3) * (j + 1) / (i->M3);
			i->All_point[this->M2 + j + 1]->x2 = i->Key_point[0]->x2;
			i->All_point[this->M2 + j + 1]->y2 = r;
		}

		//double x = pow(R5_ / R4, 1.0 / (2.0 * this->M4));
		//double dd = 2.0 * (R5_ - R4) / (i->M4 * (i->M4 + 1.0));
		//r = R4;
		for (int j = 0; j < this->M4; j++)
		{
			//double x = (j + 1.0) / (1.0 * i->M4);
			//r = R4 + (kv(x) * kv(x) + (1.0 - x) * 0.22 * x) * (R5_ - R4);
			//r = R4 * pow(kv(x), j + 1);  //R4 + (R5_ - R4) * (j + 1) / (i->M4);
			r = R4 + (R5_ - R4) * (j + 1) / (i->M4);
			i->All_point[this->M2 + this->M3 + j + 1]->x2 = i->Key_point[0]->x2;
			i->All_point[this->M2 + this->M3 + j + 1]->y2 = r;
		}


	}
}

void Setka::Move_Setka_Calculate(const double& dt)
{
	double Vx, Vy, V;
	double R2, r, R3, R4;

	//double y_Outer = this->D_Rails[16]->Key_point[1]->y;

	for (int jj = 0; jj < this->A_Rails.size(); jj++)
	{
		auto i = this->A_Rails[jj];
		// Подвинем ключевые точки
		V = i->Key_point[0]->Vx * cos(i->s) + i->Key_point[0]->Vy * sin(i->s);
		i->Key_point[0]->x2 = i->Key_point[0]->x + dt * V * cos(i->s);
		i->Key_point[0]->y2 = i->Key_point[0]->y + dt * V * sin(i->s);

		V = i->Key_point[1]->Vx * cos(i->s) + i->Key_point[1]->Vy * sin(i->s);
		i->Key_point[1]->x2 = i->Key_point[1]->x + dt * V * cos(i->s);
		i->Key_point[1]->y2 = i->Key_point[1]->y + dt * V * sin(i->s);

		V = i->Key_point[2]->Vx * cos(i->s) + i->Key_point[2]->Vy * sin(i->s);
		i->Key_point[2]->x2 = i->Key_point[2]->x + dt * V * cos(i->s);
		i->Key_point[2]->y2 = i->Key_point[2]->y + dt * V * sin(i->s);

		i->Key_point[3]->x2 = i->Key_point[3]->x;
		i->Key_point[3]->y2 = i->Key_point[3]->y;


		R2 = sqrt(kv(i->Key_point[0]->x2) + kv(i->Key_point[0]->y2));

		//double x;

		//x = pow(R11_ / R1_, 1.0 / n_inner);
		//for (int j = 0; j <= n_inner; j++) 
		//{
		//	r = R1_ * pow(x, j);
		//	//r = R1_ + (R11_ - R1_) * j / (n_inner);
		//	i->All_point[j]->x2 = r * cos(i->s);
		//	i->All_point[j]->y2 = r * sin(i->s);
		//}


		//for (int j = n_inner + 1; j <= m_inner; j++) 
		//{
		//	r = R11_ + (R111_ - R11_) * (j - n_inner) / (m_inner - n_inner);
		//	i->All_point[j]->x2 = r * cos(i->s);
		//	i->All_point[j]->y2 = r * sin(i->s);
		//}

		//for (int j = m_inner + 1; j < i->M1 - 1; j++)  // Передвинули точки до ударной волны
		//{
		//	//r = R11_ * pow(x, j - 7);
		//	r = R111_ + (R2 - R111_) * (j - m_inner) / (i->M1 - 1 - m_inner);
		//	i->All_point[j]->x2 = r * cos(i->s);
		//	i->All_point[j]->y2 = r * sin(i->s);
		//}

		double x;

		//x = pow(R2 / R1_, 1.0 / (i->M1 - 1));
		for (int j = 0; j < n_inner; j++)
		{
			x = j * 1.0 / (n_inner);
			r = R1_ + pow(x, 1.3) * (R111_ - R1_);
			//r = R1_ + (R11_ - R1_) * j / (n_inner);
			i->All_point[j]->x2 = r * cos(i->s);
			i->All_point[j]->y2 = r * sin(i->s);
		}

		for (int j = n_inner; j < i->M1 - 1; j++)
		{
			x = (j - n_inner) * 1.0 / (i->M1 - 1 - n_inner);
			r = R111_ + pow(x, 1.05) * (R2 - R111_);
			//r = R1_ + (R11_ - R1_) * j / (n_inner);
			i->All_point[j]->x2 = r * cos(i->s);
			i->All_point[j]->y2 = r * sin(i->s);
		}



		//R3 = sqrt(kv(i->Key_point[1]->x2) + kv(i->Key_point[1]->y2)) - zazor;
		R3 = sqrt(kv(i->Key_point[1]->x2) + kv(i->Key_point[1]->y2));

		for (int j = 0; j < i->M2 - 1; j++)  // Передвинули точки до контакта
		{
			r = R2 + (R3 - R2) * (j + 1) / (i->M2);         // Равномерное распределение точек
			i->All_point[i->M1 + j]->x2 = r * cos(i->s);
			i->All_point[i->M1 + j]->y2 = r * sin(i->s);
		}

		//i->All_point[i->M1 + i->M2 - 2]->x2 = R3 * cos(i->s);   // Устанавливаем точку до контакта 
		//i->All_point[i->M1 + i->M2 - 2]->y2 = R3 * sin(i->s);
		R3 = sqrt(kv(i->Key_point[1]->x2) + kv(i->Key_point[1]->y2));

		//i->All_point[i->M1 + i->M2]->x2 = R3 * cos(i->s);   // Устанавливаем точку после контакта 
		//i->All_point[i->M1 + i->M2]->y2 = R3 * sin(i->s);

		R4 = sqrt(kv(i->Key_point[2]->x2) + kv(i->Key_point[2]->y2));


		for (int j = 0; j < i->M3 - 1; j++)// Передвинули точки до внешней волны
		{
			//x = (j) / (1.0 * (i->M3 - 1));
			//r = R3 + (s_k * x + 3.0 * x * x - 3.0 * s_k * x * x - 2.0 * x * x * x + 2.0 * s_k * x * x * x) * (R4 - R3);
			r = R3 + (R4 - R3) * (j + 1) / (i->M3);         // Равномерное распределение точек
			i->All_point[i->M1 + i->M2 + j]->x2 = r * cos(i->s);
			i->All_point[i->M1 + i->M2 + j]->y2 = r * sin(i->s);
		}

		R4 = sqrt(kv(i->Key_point[2]->x2) + kv(i->Key_point[2]->y2));

		//x = log(R5_ / R4) / (log(1.0 * i->M4 + 1) * (i->M4));
		//double dd = 2.0 * (R5_ - R4) / (i->M4 * (i->M4 + 1.0));
		//r = R4;
		for (int j = 0; j < i->M4; j++)
		{
			//x = (j + 1.0) / (1.0 * i->M4);
			//r = R4 + (kv(x) * kv(x) + (1.0 - x) * 0.03 * x) * (R5_ - R4);
			//r = R4 *  pow( 1.0 * (j + 2), x * (j + 1));  //R4 + (R5_ - R4) * (j + 1) / (i->M4);
			//r = R4 + (R5_ - R4) * (j) / (i->M4 - 1);   // равномерно
			x = (j + 1) * 1.0 / (i->M4);
			double dd = 1.5708 / sqrt(1.0 - 1.0 / 1.76);
			r = R4 + pow(x, 1.76 * (1.0 - kv(i->s / dd))) * (R5_ - R4);
			//r = r + (j + 1) * dd;
			i->All_point[i->M1 + i->M2 + i->M3 + j]->x2 = r * cos(i->s);
			i->All_point[i->M1 + i->M2 + i->M3 + j]->y2 = r * sin(i->s);
		}
	}

	for (int jj = 0; jj < this->B_Rails.size(); jj++)
	{
		auto i = this->B_Rails[jj];
		// Подвинем ключевые точки
		V = i->Key_point[0]->Vx * cos(i->s) + i->Key_point[0]->Vy * sin(i->s);
		i->Key_point[0]->x2 = i->Key_point[0]->x + dt * V * cos(i->s);
		i->Key_point[0]->y2 = i->Key_point[0]->y + dt * V * sin(i->s);

		V = i->Key_point[1]->Vx * 0.0 + i->Key_point[1]->Vy * 1.0;
		i->Key_point[1]->x2 = i->Key_point[0]->x2;
		i->Key_point[1]->y2 = i->Key_point[1]->y + dt * V;

		V = i->Key_point[2]->Vx * 0.0 + i->Key_point[2]->Vy * 1.0;
		i->Key_point[2]->x2 = i->Key_point[0]->x2;
		i->Key_point[2]->y2 = i->Key_point[2]->y + dt * V;

		i->Key_point[3]->x2 = i->Key_point[0]->x2;
		i->Key_point[3]->y2 = i->Key_point[3]->y;

		R2 = sqrt(kv(i->Key_point[0]->x2) + kv(i->Key_point[0]->y2));

		//double x;
		//x = pow(R11_ / R1_, 1.0 / n_inner);
		//for (int j = 0; j <= n_inner; j++)  // Передвинули точки до ударной волны
		//{
		//	r = R1_ * pow(x, j);
		//	//r = R1_ + (R11_ - R1_) * j / (n_inner);
		//	i->All_point[j]->x2 = r * cos(i->s);
		//	i->All_point[j]->y2 = r * sin(i->s);
		//}


		//for (int j = n_inner + 1; j <= m_inner; j++)  // Передвинули точки до ударной волны
		//{
		//	r = R11_ + (R111_ - R11_) * (j - n_inner) / (m_inner - n_inner);
		//	i->All_point[j]->x2 = r * cos(i->s);
		//	i->All_point[j]->y2 = r * sin(i->s);
		//}

		//for (int j = m_inner + 1; j < i->M1 - 1; j++)  // Передвинули точки до ударной волны
		//{
		//	//r = R11_ * pow(x, j - 7);
		//	r = R111_ + (R2 - R111_) * (j - m_inner) / (i->M1 - 1 - m_inner);
		//	i->All_point[j]->x2 = r * cos(i->s);
		//	i->All_point[j]->y2 = r * sin(i->s);
		//}

		double x;

		for (int j = 0; j < n_inner; j++)
		{
			x = j * 1.0 / (n_inner);
			r = R1_ + pow(x, 1.3) * (R111_ - R1_);
			//r = R1_ + (R11_ - R1_) * j / (n_inner);
			i->All_point[j]->x2 = r * cos(i->s);
			i->All_point[j]->y2 = r * sin(i->s);
		}

		for (int j = n_inner; j < i->M1 - 1; j++)
		{
			x = (j - n_inner) * 1.0 / (i->M1 - 1 - n_inner);
			r = R111_ + pow(x, 1.05) * (R2 - R111_);
			//r = R1_ + (R11_ - R1_) * j / (n_inner);
			i->All_point[j]->x2 = r * cos(i->s);
			i->All_point[j]->y2 = r * sin(i->s);
		}



		R2 = i->Key_point[0]->y2;
		R3 = i->Key_point[1]->y2;

		for (int j = 0; j < i->M2 - 1; j++)
		{
			r = R2 + (R3 - R2) * (j + 1) / (i->M2);
			i->All_point[i->M1 + j]->x2 = i->Key_point[0]->x2;
			i->All_point[i->M1 + j]->y2 = r;
		}

		//i->All_point[i->M1 + i->M2 - 2]->x2 = i->Key_point[0]->x2;   // Устанавливаем точку до контакта 
		//i->All_point[i->M1 + i->M2 - 2]->y2 = R3 ;
		R3 = i->Key_point[1]->y2;

		//i->All_point[i->M1 + i->M2]->x2 = i->Key_point[0]->x2;   // Устанавливаем точку после контакта 
		//i->All_point[i->M1 + i->M2]->y2 = R3 ;

		R4 = i->Key_point[2]->y2;

		for (int j = 0; j < i->M3 - 1; j++)
		{
			//r = R3 + (R4 - R3) * (j) / (i->M3 - 1);
			//x = (j) / (1.0 * (i->M3 - 1));
			//r = R3 + (s_k * x + 3.0 * x * x - 3.0 * s_k * x * x - 2.0 * x * x * x + 2.0 * s_k * x * x * x) * (R4 - R3);
			r = R3 + (R4 - R3) * (j + 1) / (i->M3);
			i->All_point[i->M1 + i->M2 + j]->x2 = i->Key_point[0]->x2;
			i->All_point[i->M1 + i->M2 + j]->y2 = r;
		}

		//x = pow(R5_ / R4, 1.0 / (2.0 * this->M4));
		//double dd = 2.0 * (R5_ - R4) / (i->M4 * (i->M4 + 1.0));
		//r = R4;

		for (int j = 0; j < i->M4; j++)
		{
			//x = (j + 1.0) / (1.0 * i->M4);
			//r = R4 + (kv(x) * kv(x) + (1.0 - x) * (0.03 + 0.19 * jj / (this->B_Rails.size() - 1.0)) * x) * (R5_ - R4);
			//r = R4 * pow(kv(x), j + 1);  //R4 + (R5_ - R4) * (j + 1) / (i->M4);
			r = R4 + (R5_ - R4) * (j + 1) / (i->M4);
			i->All_point[i->M1 + i->M2 + i->M3 + j]->x2 = i->Key_point[0]->x2;
			i->All_point[i->M1 + i->M2 + i->M3 + j]->y2 = r;
		}

	}

	for (auto& i : this->C_Rails)
	{
		V = i->Key_point[0]->Vx * cos(i->s) + i->Key_point[0]->Vy * sin(i->s);
		i->Key_point[0]->x2 = i->Key_point[0]->x + dt * V * cos(i->s);
		i->Key_point[0]->y2 = i->Key_point[0]->y + dt * V * sin(i->s);

		i->Key_point[1]->x2 = Left_; // i->Key_point[1]->x;
		i->Key_point[1]->y2 = i->Key_point[1]->y;

		R2 = sqrt(kv(i->Key_point[0]->x2) + kv(i->Key_point[0]->y2));

		//double x;
		//x = pow(R11_ / R1_, 1.0 / n_inner);
		//for (int j = 0; j <= n_inner; j++)  // Передвинули точки до ударной волны
		//{
		//	r = R1_ * pow(x, j);
		//	//r = R1_ + (R11_ - R1_) * j / (n_inner);
		//	i->All_point[j]->x2 = r * cos(i->s);
		//	i->All_point[j]->y2 = r * sin(i->s);
		//}


		//for (int j = n_inner + 1; j <= m_inner; j++)  // Передвинули точки до ударной волны
		//{
		//	r = R11_ + (R111_ - R11_) * (j - n_inner) / (m_inner - n_inner);
		//	i->All_point[j]->x2 = r * cos(i->s);
		//	i->All_point[j]->y2 = r * sin(i->s);
		//}

		//for (int j = m_inner + 1; j < i->M1 - 1; j++)  // Передвинули точки до ударной волны
		//{
		//	//r = R11_ * pow(x, j - 7);
		//	r = R111_ + (R2 - R111_) * (j - m_inner) / (i->M1 - 1 - m_inner);
		//	i->All_point[j]->x2 = r * cos(i->s);
		//	i->All_point[j]->y2 = r * sin(i->s);
		//}

		double x;

		for (int j = 0; j < n_inner; j++)
		{
			x = j * 1.0 / (n_inner);
			r = R1_ + pow(x, 1.3) * (R111_ - R1_);
			//r = R1_ + (R11_ - R1_) * j / (n_inner);
			i->All_point[j]->x2 = r * cos(i->s);
			i->All_point[j]->y2 = r * sin(i->s);
		}

		for (int j = n_inner; j < i->M1 - 1; j++)
		{
			x = (j - n_inner) * 1.0 / (i->M1 - 1 - n_inner);
			r = R111_ + pow(x, 1.05) * (R2 - R111_);
			//r = R1_ + (R11_ - R1_) * j / (n_inner);
			i->All_point[j]->x2 = r * cos(i->s);
			i->All_point[j]->y2 = r * sin(i->s);
		}

		R2 = i->Key_point[0]->x2;

		R3 = i->Key_point[1]->x2;

		x = pow(R3 / (R2), 1.0 / (i->M2));

		for (int j = i->M1; j < i->All_point.size(); j++)
		{
			r = R2 * pow(x, j + 1 - i->M1);
			i->All_point[j]->x2 = r;
			i->All_point[j]->y2 = i->Key_point[0]->y2;
		}
	}

	for (int jj = 0; jj < this->D_Rails.size(); jj++)
	{
		double x;
		// Подвинем ключевые точки
		auto i = this->D_Rails[jj];
		V = i->Key_point[1]->Vy;
		i->Key_point[1]->x2 = i->Key_point[0]->x2;
		/*if (jj > 16)
		{
			i->Key_point[1]->y2 = y_Outer + (247.0 - y_Outer) * (jj - 16.0) / (this->D_Rails.size() - 1.0 - 16.0);
		}
		else
		{
			i->Key_point[1]->y2 = i->Key_point[1]->y + dt * V;
		}*/
		i->Key_point[1]->y2 = i->Key_point[1]->y + dt * V;
		/*if (jj > 22)
		{
			i->Key_point[1]->y2 = this->D_Rails[22]->Key_point[1]->y2;
		}*/

		V = i->Key_point[2]->Vy;
		i->Key_point[2]->x2 = i->Key_point[0]->x2;
		i->Key_point[2]->y2 = i->Key_point[2]->y + dt * V;

		i->Key_point[3]->x2 = i->Key_point[0]->x2;
		i->Key_point[3]->y2 = i->Key_point[3]->y;

		R2 = i->Key_point[0]->y2;

		R3 = i->Key_point[1]->y2;

		for (int j = 0; j < this->M2 - 1; j++)
		{
			r = R2 + (R3 - R2) * (j + 1) / (i->M2);
			i->All_point[j + 1]->x2 = i->Key_point[0]->x2;
			i->All_point[j + 1]->y2 = r;
		}

		//i->All_point[1 + i->M2 - 2]->x2 = i->Key_point[0]->x2;   // Устанавливаем точку до контакта 
		//i->All_point[1 + i->M2 - 2]->y2 = R3;
		R3 = i->Key_point[1]->y2;

		//i->All_point[1 + i->M2]->x2 = i->Key_point[0]->x2;   // Устанавливаем точку после контакта 
		//i->All_point[1 + i->M2]->y2 = R3;

		R4 = i->Key_point[2]->y2;

		for (int j = 0; j < this->M3 - 1; j++)
		{
			//r = R3 + (R4 - R3) * (j) / (i->M3 - 1);
			//x = (j) / (1.0 * (i->M3 - 1));
			//r = R3 + (s_k * x + 3.0 * x * x - 3.0 * s_k * x * x - 2.0 * x * x * x + 2.0 * s_k * x * x * x) * (R4 - R3);
			r = R3 + (R4 - R3) * (j + 1) / (i->M3);
			i->All_point[this->M2 + j + 1]->x2 = i->Key_point[0]->x2;
			i->All_point[this->M2 + j + 1]->y2 = r;
		}

		//double x = pow(R5_ / R4, 1.0 / (2.0 * this->M4));
		//double dd = 2.0 * (R5_ - R4) / (i->M4 * (i->M4 + 1.0));
		//r = R4;
		for (int j = 0; j < this->M4; j++)
		{
			//double x = (j + 1.0) / (1.0 * i->M4);
			//r = R4 + (kv(x) * kv(x) + (1.0 - x) * 0.22 * x) * (R5_ - R4);
			//r = R4 * pow(kv(x), j + 1);  //R4 + (R5_ - R4) * (j + 1) / (i->M4);
			r = R4 + (R5_ - R4) * (j + 1) / (i->M4);
			i->All_point[this->M2 + this->M3 + j + 1]->x2 = i->Key_point[0]->x2;
			i->All_point[this->M2 + this->M3 + j + 1]->y2 = r;
		}


	}
}

void Setka::Init_conditions(void)
{
	double x, y, dist;


	for (auto& i : this->All_Cells)
	{
		i->Get_Center(x, y);
		dist = sqrt(x * x + y * y);

		if (dist <= 0.1)
		{
			double ro = 1.0 / (kv(chi_real) * dist * dist);
			double P_E = ro * chi_real * chi_real / (ggg * 10.0 * 10.0);

			i->par[0].ro = ro;
			i->par[0].p = P_E;
			i->par[0].u = chi_real * x / dist;
			i->par[0].v = chi_real * y / dist;
			i->par[0].Q = ro;
			i->par[1] = i->par[0];
		}
	}
}

void Setka::Go_stationary(int step)
{
	cout << "Setka.cpp    " << "START " << step << endl;
	int now1 = 1;
	int now2 = 0;
	double T[2];
	T[0] = T[1] = 0.00000001;
	mutex mut;

	for (int st = 0; st < step; st ++)
	{
		if (st % 1000 == 0 && st > 0)
		{
			cout << "Setka.cpp    " << st << " " << T[now2] << endl;
		}
		now1 = (now1 + 1) % 2; // Какие параметры сейчас берём
		now2 = (now2 + 1) % 2; // Какие параметры сейчас меняем
		double time = 10000000;
		T[now2] = 100000000;

#pragma omp parallel for
		for (int num_cell = 0; num_cell < this->All_Cells.size(); num_cell++)
		{
			auto K = this->All_Cells[num_cell];
			K->Potok[0] = K->Potok[1] = K->Potok[2] = K->Potok[3] = K->Potok[4] = 0.0;
			double P[4];
			double PQ;
			Parametr par1 = K->par[now1];
			double n1, n2;
			double dist;
			double x, y;
			K->Get_Center(x, y);
			double rad = sqrt(kv(x) + kv(y));

			if (K->type == C_centr)
			{
				continue;
			}

			for (auto& i : K->Grans)
			{
				double x2, y2;
				double S = i->Get_square();
				i->Get_Center(x2, y2);
				dist = sqrt(kv(x - x2) + kv(y - y2));
				Parametr par2, par11; 
				par11 = par1;
				i->Get_par(par2, now1);
				i->Get_normal(n1, n2);
				double Vc, Vl, Vp;
				int met = 1;
				if (rad < 95 && i->type == Usualy)
				{
					double x3, y3;
					i->Sosed->Get_Center(x3, y3);
					double dd = sqrt(kv(x3) + kv(y3));
					par11.ro = par11.ro * kv(rad) / (kv(x2) + kv(y2));
					par2.ro = par2.ro * kv(dd) / (kv(x2) + kv(y2));
					par11.p = par11.p * pow(rad / sqrt(kv(x2) + kv(y2)), 2 * ggg);
					par2.p = par2.p * pow(dd / sqrt(kv(x2) + kv(y2)), 2 * ggg);
					polar_perenos(x, y, x2, y2, par11.u, par11.v);
					polar_perenos(x3, y3, x2, y2, par2.u, par2.v);
				}
				else if (i->type == Inner_sphere)
				{
					par11.ro = par11.ro * kv(rad) / (kv(x2) + kv(y2));
					par11.p = par11.p * pow(rad / sqrt(kv(x2) + kv(y2)), 2 * ggg);
				}
				else if (rad < 95 && i->type == Axis)
				{
					//par11.ro = par11.ro * kv(rad) / (kv(x2) + kv(y2));
					//par11.p = par11.p * pow(rad / sqrt(kv(x2) + kv(y2)), 2 * ggg);
					polar_perenos(x, y, x2, y2, par11.u, par11.v);
					par11.v = 0.0;
					par2 = par11;
					//par2.v = -par2.v;
				}
				else if (i->type == Axis)
				{
					par11.v = 0.0;
					par2 = par11;
				}
				/*if (num_cell == 1)
				{
					cout << "Setka.cpp    " << "-----" << endl;
					cout << "Setka.cpp    " << n1 << " " << n2 << endl;
					cout << "Setka.cpp    " << par11.ro << " " << par11.u << " " << par11.v << endl;
					cout << "Setka.cpp    " << par2.ro << " " << par2.u << " " << par2.v << endl;
				}*/
				/*if (step > 20000)
				{
					met = 1;
				}*/
				time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro, par11.Q, par11.p, par11.u, par11.v, par2.ro, par2.Q, //
							par2.p, par2.u, par2.v, 0.0, P, PQ, n1, n2, dist, met, Vl, Vc, Vp));
				for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
				{
					K->Potok[k] = K->Potok[k] + P[k] * S;
				}
				//cout << "Setka.cpp    " << "B5" << endl;
				K->Potok[4] = K->Potok[4] + PQ * S;
				//cout << "Setka.cpp    " << "B6" << endl;
			}

			if (time < T[now2])
			{
				mut.lock();
				T[now2] = time;
				mut.unlock();
			}

			double ro3, p3, u3, v3, Q33;
			double Volume = K->Get_Volume();

			ro3 = par1.ro - T[now1] * (K->Potok[0] / Volume + par1.ro * par1.v / y);
			Q33 = par1.Q - (T[now1] / Volume) * K->Potok[4] - T[now1] * par1.Q * par1.v / y;
			if (ro3 <= 0)
			{
				printf("Problemsssss  par1.ro < 0!   %lf,   %lf,   %lf\n", x, y, par1.ro);
				ro3 = 0.00001;
			}
			u3 = (par1.ro * par1.u - T[now1] * (K->Potok[1] / Volume + par1.ro * par1.v * par1.u / y)) / ro3;
			v3 = (par1.ro * par1.v - T[now1] * (K->Potok[2] / Volume + par1.ro * par1.v * par1.v / y)) / ro3;
			p3 = (((par1.p / (ggg - 1) + par1.ro * (par1.u * par1.u + par1.v * par1.v) * 0.5) - T[now1] * (K->Potok[3] / Volume +  //
				+ par1.v * (ggg * par1.p / (ggg - 1) + par1.ro * (par1.u * par1.u + par1.v * par1.v) * 0.5) / y)) - //
				0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
			if (p3 <= 0)
			{
				p3 = 0.000001;
			}

			K->par[now2].ro = ro3;
			K->par[now2].Q = Q33;
			K->par[now2].p = p3;
			K->par[now2].u = u3;
			K->par[now2].v = v3;
		}
	}
}

void Setka::Go_stationary_TVD(int step)
{
	cout << "Setka.cpp    " << "START " << step << endl;
	int now1 = 1;
	int now2 = 0;
	double T[2];
	T[0] = T[1] = 0.00000001;
	mutex mut;

	for (int st = 0; st < step; st++)
	{
		if (st % 1000 == 0 && st > 0)
		{
			cout << "Setka.cpp    " << st << " " << T[now2] << endl;
		}
		now1 = (now1 + 1) % 2; // Какие параметры сейчас берём
		now2 = (now2 + 1) % 2; // Какие параметры сейчас меняем
		double time = 10000000;
		T[now2] = 100000000;

#pragma omp parallel for
		for (int num_cell = 0; num_cell < this->All_Cells.size(); num_cell++)
		{
			auto K = this->All_Cells[num_cell];
			K->Potok[0] = K->Potok[1] = K->Potok[2] = K->Potok[3] = K->Potok[4] = 0.0;
			double P[4];
			double PQ;
			Parametr par1 = K->par[now1];
			double n1, n2;
			double dist;
			double x, y;
			K->Get_Center(x, y);

			if (K->type == C_centr)
			{
				continue;
			}

			for (auto& i : K->Grans)
			{
				double x2, y2;
				double S = i->Get_square();
				i->Get_Center(x2, y2);
				dist = sqrt(kv(x - x2) + kv(y - y2));
				Parametr par2;
				Parametr par3;
				if (i->type == Usualy)
				{
					i->Get_par_TVD(par3, now1);
					i->Gran_copy->Get_par_TVD(par2, now1);
				}
				else
				{
					par3 = par1;
					i->Get_par(par2, now1);
				}
				i->Get_normal(n1, n2);
				double Vc, Vl, Vp;
				int met = 0;
				/*if (step > 20000)
				{
					met = 1;
				}*/
				/*if (num_cell == 1)
				{
					cout << "Setka.cpp    " << par3.ro << " " << par2.ro << " " << x2 << " " << y2 << endl;
				}*/
				time = min(time, this->HLLC_2d_Korolkov_b_s(par3.ro, par3.Q, par3.p, par3.u, par3.v, par2.ro, par2.Q, //
					par2.p, par2.u, par2.v, 0.0, P, PQ, n1, n2, dist, met, Vl, Vc, Vp));
				for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
				{
					K->Potok[k] = K->Potok[k] + P[k] * S;
				}
				//cout << "Setka.cpp    " << "B5" << endl;
				K->Potok[4] = K->Potok[4] + PQ * S;
				//cout << "Setka.cpp    " << "B6" << endl;
			}

			if (time < T[now2])
			{
				mut.lock();
				T[now2] = time;
				mut.unlock();
			}

			double ro3, p3, u3, v3, Q33;
			double Volume = K->Get_Volume();

			ro3 = par1.ro - T[now1] * (K->Potok[0] / Volume + par1.ro * par1.v / y);
			Q33 = par1.Q - (T[now1] / Volume) * K->Potok[4] - T[now1] * par1.Q * par1.v / y;
			if (ro3 <= 0)
			{
				printf("Problemsssss  par1.ro < 0!   %lf,   %lf,   %lf\n", x, y, par1.ro);
				ro3 = 0.00001;
			}
			u3 = (par1.ro * par1.u - T[now1] * (K->Potok[1] / Volume + par1.ro * par1.v * par1.u / y)) / ro3;
			v3 = (par1.ro * par1.v - T[now1] * (K->Potok[2] / Volume + par1.ro * par1.v * par1.v / y)) / ro3;
			p3 = (((par1.p / (ggg - 1) + par1.ro * (par1.u * par1.u + par1.v * par1.v) * 0.5) - T[now1] * (K->Potok[3] / Volume +  //
				+par1.v * (ggg * par1.p / (ggg - 1) + par1.ro * (par1.u * par1.u + par1.v * par1.v) * 0.5) / y)) - //
				0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
			if (p3 <= 0)
			{
				p3 = 0.000001;
			}

			K->par[now2].ro = ro3;
			K->par[now2].Q = Q33;
			K->par[now2].p = p3;
			K->par[now2].u = u3;
			K->par[now2].v = v3;
		}
	}
}

void Setka::Go_stationary_5_komponent(int step)
{
	cout << "Setka.cpp    " << "START " << step << endl;
	int now1 = 1;
	int now2 = 0;
	double T[2];
	T[0] = T[1] = 0.00000001;
	mutex mut;
	double x_min, y_min;

	for (int st = 0; st < step; st++)
	{
		if (st % 1000 == 0 && st > 0)
		{
			cout << "Setka.cpp    " << st << " " << T[now2] << " " << x_min << " " << y_min << endl;
		}
		now1 = (now1 + 1) % 2; // Какие параметры сейчас берём
		now2 = (now2 + 1) % 2; // Какие параметры сейчас меняем
		double time = 10000000;
		T[now2] = 100000000;

		//omp_set_num_threads(4);
#pragma omp parallel for
		for (int num_cell = 0; num_cell < this->All_Cells.size(); num_cell++)
		{
			auto K = this->All_Cells[num_cell];
			K->Potok[0] = K->Potok[1] = K->Potok[2] = K->Potok[3] = K->Potok[4] = 0.0;
			K->Potok_H1[0] = K->Potok_H1[1] = K->Potok_H1[2] = K->Potok_H1[3] = 0.0;
			K->Potok_H2[0] = K->Potok_H2[1] = K->Potok_H2[2] = K->Potok_H2[3] = 0.0;
			K->Potok_H3[0] = K->Potok_H3[1] = K->Potok_H3[2] = K->Potok_H3[3] = 0.0;
			K->Potok_H4[0] = K->Potok_H4[1] = K->Potok_H4[2] = K->Potok_H4[3] = 0.0;

			double P[4];
			double PQ;
			Parametr par1 = K->par[now1];
			double n1, n2;
			double dist;
			double x, y;
			K->Get_Center(x, y);
			double radius = sqrt(kv(x) + kv(y));
			if (radius < R11_)
			{
				continue;
			}
			double Volume = K->Get_Volume();

			for (auto& i : K->Grans)
			{
				double x2, y2;
				double S = i->Get_square();
				i->Get_Center(x2, y2);
				dist = sqrt(kv(x - x2) + kv(y - y2));
				Parametr par2, par11;
				par11 = par1;
				i->Get_par(par2, now1);
				i->Get_normal(n1, n2);
				double Vc, Vl, Vp;
				int met = 1;
				double dis = sqrt(kv(x2) + kv(y2));
				if (radius < 85 && i->type == Usualy)
				{
					double x3, y3;
					i->Sosed->Get_Center(x3, y3);
					double dd = sqrt(kv(x3) + kv(y3));
					par11.ro = par11.ro * kv(radius) / kv(dis);
					par2.ro = par2.ro * kv(dd) / kv(dis);
					par11.ro_H1 = par11.ro_H1 * pow(radius / dis, H_pow);
					par2.ro_H1 = par2.ro_H1 * pow(dd / dis, H_pow);
					par11.Q = par11.Q * kv(radius) / kv(dis);
					par2.Q = par2.Q * kv(dd) / kv(dis);
					par11.p = par11.p * pow(radius / dis, 2 * ggg);
					par2.p = par2.p * pow(dd / dis, 2 * ggg);
					polar_perenos(x, y, x2, y2, par11.u, par11.v);
					polar_perenos(x3, y3, x2, y2, par2.u, par2.v);
					polar_perenos(x, y, x2, y2, par11.u_H1, par11.v_H1);
					polar_perenos(x3, y3, x2, y2, par2.u_H1, par2.v_H1);
				}
				else if (i->type == Inner_sphere)
				{
					par11.ro = par11.ro * kv(radius) / kv(dis);
					par11.ro_H1 = par11.ro_H1 * pow(radius / dis, H_pow);
					par11.Q = par11.Q * kv(radius) / kv(dis);
					par11.p = par11.p * pow(radius / dis, 2 * ggg);
				}
				else if (radius < 85 && i->type == Axis)
				{
					par11.ro = par11.ro * kv(radius) / kv(dis);
					par11.ro_H1 = par11.ro_H1 * pow(radius / dis, H_pow);
					par11.Q = par11.Q * kv(radius) / kv(dis);
					par11.p = par11.p * pow(radius / dis, 2 * ggg);
					polar_perenos(x, y, x2, y2, par11.u, par11.v);
					polar_perenos(x, y, x2, y2, par11.u_H1, par11.v_H1);
					//par11.v = 0.0;
					par2 = par11;
					par2.v = -par2.v;
					par2.v_H1 = -par2.v_H1;
					par2.v_H2 = -par2.v_H2;
					par2.v_H3 = -par2.v_H3;
					par2.v_H4 = -par2.v_H4;
				}

				if (K->type != C_centr)
				{
					time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro_H1, par11.Q, par11.p_H1, par11.u_H1, par11.v_H1, par2.ro_H1, par2.Q, //
						par2.p_H1, par2.u_H1, par2.v_H1, 0.0, P, PQ, n1, n2, dist, met, Vl, Vc, Vp));
					for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
					{
						K->Potok_H1[k] = K->Potok_H1[k] + P[k] * S;
					}
				}

				time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro_H2, par11.Q, par11.p_H2, par11.u_H2, par11.v_H2, par2.ro_H2, par2.Q, //
					par2.p_H2, par2.u_H2, par2.v_H2, 0.0, P, PQ, n1, n2, dist, met, Vl, Vc, Vp));
				for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
				{
					K->Potok_H2[k] = K->Potok_H2[k] + P[k] * S;
				}

				time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro_H3, par11.Q, par11.p_H3, par11.u_H3, par11.v_H3, par2.ro_H3, par2.Q, //
					par2.p_H3, par2.u_H3, par2.v_H3, 0.0, P, PQ, n1, n2, dist, met, Vl, Vc, Vp));
				for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
				{
					K->Potok_H3[k] = K->Potok_H3[k] + P[k] * S;
				}

				time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro_H4, par11.Q, par11.p_H4, par11.u_H4, par11.v_H4, par2.ro_H4, par2.Q, //
					par2.p_H4, par2.u_H4, par2.v_H4, 0.0, P, PQ, n1, n2, dist, met, Vl, Vc, Vp));
				for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
				{
					K->Potok_H4[k] = K->Potok_H4[k] + P[k] * S;
				}
				

				if (K->type != C_centr)
				{
					time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro, par11.Q, par11.p, par11.u, par11.v, par2.ro, par2.Q, //
						par2.p, par2.u, par2.v, 0.0, P, PQ, n1, n2, dist, met, Vl, Vc, Vp));
					for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
					{
						K->Potok[k] = K->Potok[k] + P[k] * S;
					}
					K->Potok[4] = K->Potok[4] + PQ * S;
				}
			}

			if (time < T[now2])
			{
				mut.lock();
				T[now2] = time;
				x_min = x;
				y_min = y;
				mut.unlock();
			}



			double U_M_H1, U_M_H2, U_M_H3, U_M_H4;
			double U_H1, U_H2, U_H3, U_H4;
			double sigma_H1, sigma_H2, sigma_H3, sigma_H4;
			double nu_H1, nu_H2, nu_H3, nu_H4;
			double q2_1, q2_2, q3;
			double u, v, ro, p, Q;
			if (par1.Q / par1.ro < 90)
			{
				u = par1.u * (chi_real / chi_);
				v = par1.v * (chi_real / chi_);
				ro = par1.ro / kv(chi_real / chi_);
				p = par1.p;
				Q = par1.Q / kv(chi_real / chi_);
			}
			else
			{
				u = par1.u;
				v = par1.v;
				ro = par1.ro;
				p = par1.p;
				Q = par1.Q;
			}
			double u_H1 = par1.u_H1, v_H1 = par1.v_H1, ro_H1 = par1.ro_H1, p_H1 = par1.p_H1;
			double u_H2 = par1.u_H2, v_H2 = par1.v_H2, ro_H2 = par1.ro_H2, p_H2 = par1.p_H2;
			double u_H3 = par1.u_H3, v_H3 = par1.v_H3, ro_H3 = par1.ro_H3, p_H3 = par1.p_H3;
			double u_H4 = par1.u_H4, v_H4 = par1.v_H4, ro_H4 = par1.ro_H4, p_H4 = par1.p_H4;

			U_M_H1 = sqrt(kv(u - u_H1) + kv(v - v_H1) + (64.0 / (9.0 * pi_)) //
				* (p / ro + 2.0 * p_H1 / ro_H1));
			U_M_H2 = sqrt(kv(u - u_H2) + kv(v - v_H2) + (64.0 / (9.0 * pi_)) //
				* (p / ro + 2.0 * p_H2 / ro_H2));
			U_M_H3 = sqrt(kv(u - u_H3) + kv(v - v_H3) + (64.0 / (9.0 * pi_)) //
				* (p / ro + 2.0 * p_H3 / ro_H3));
			U_M_H4 = sqrt(kv(u - u_H4) + kv(v - v_H4) + (64.0 / (9.0 * pi_)) //
				* (p / ro + 2.0 * p_H4 / ro_H4));

			U_H1 = sqrt(kv(u - u_H1) + kv(v - v_H1) + (4.0 / pi_) //
				* (p / ro + 2.0 * p_H1 / ro_H1));
			U_H2 = sqrt(kv(u - u_H2) + kv(v - v_H2) + (4.0 / pi_) //
				* (p / ro + 2.0 * p_H2 / ro_H2));
			U_H3 = sqrt(kv(u - u_H3) + kv(v - v_H3) + (4.0 / pi_) //
				* (p / ro + 2.0 * p_H3 / ro_H3));
			U_H4 = sqrt(kv(u - u_H4) + kv(v - v_H4) + (4.0 / pi_) //
				* (p / ro + 2.0 * p_H4 / ro_H4));

			sigma_H1 = kv(1.0 - a_2 * log(U_M_H1)); // 0.1243
			sigma_H2 = kv(1.0 - a_2 * log(U_M_H2));
			sigma_H3 = kv(1.0 - a_2 * log(U_M_H3)); // 0.1121     a_2
			sigma_H4 = kv(1.0 - a_2 * log(U_M_H4));

			nu_H1 = ro * ro_H1 * U_M_H1 * sigma_H1;
			nu_H2 = ro * ro_H2 * U_M_H2 * sigma_H2;
			nu_H3 = ro * ro_H3 * U_M_H3 * sigma_H3;
			nu_H4 = ro * ro_H4 * U_M_H4 * sigma_H4;

			q2_1 = (n_p_LISM_ / Kn_) * (nu_H1 * (u_H1 - u) + nu_H2 * (u_H2 - u) //
				+ nu_H3 * (u_H3 - u) + nu_H4 * (u_H4 - u));
			q2_2 = (n_p_LISM_ / Kn_) * (nu_H1 * (v_H1 - v) + nu_H2 * (v_H2 - v) //
				+ nu_H3 * (v_H3 - v) + nu_H4 * (v_H4 - v));
			q3 = (n_p_LISM_ / Kn_) * (nu_H1 * ((kv(u_H1) + kv(v_H1) - kv(u) - kv(v)) / 2.0 + //
				(U_H1 / U_M_H1) * (2.0 * p_H1 / ro_H1 - p / ro)) + //
				nu_H2 * ((kv(u_H2) + kv(v_H2) - kv(u) - kv(v)) / 2.0 + //
					(U_H2 / U_M_H2) * (2.0 * p_H2 / ro_H2 - p / ro)) + //
				nu_H3 * ((kv(u_H3) + kv(v_H3) - kv(u) - kv(v)) / 2.0 + //
					(U_H3 / U_M_H3) * (2.0 * p_H3 / ro_H3 - p / ro)) + //
				nu_H4 * ((kv(u_H4) + kv(v_H4) - kv(u) - kv(v)) / 2.0 + //
					(U_H4 / U_M_H4) * (2.0 * p_H4 / ro_H4 - p / ro))) / (chi_real / chi_);


			/*q2_1 = 0.0;
			q2_2 = 0.0;
			q3 = 0.0;*/

			double ro2_p, p2_p, V1_p, V2_p, QQ2;
			double ro3, p3, u3, v3, Q33;
			double kappa = Q / ro;

			u = par1.u;
			v = par1.v;
			ro = par1.ro;
			p = par1.p;
			Q = par1.Q;

			if (K->type != C_centr)
			{
				ro3 = ro - T[now1] * (K->Potok[0] / Volume + ro * v / y);
				Q33 = Q - (T[now1] / Volume) * K->Potok[4] - T[now1] * Q * v / y;
				if (ro3 <= 0)
				{
					printf("Problemsssss  ro < 0!   %lf,  %lf\n", x, y);
					ro3 = 0.00001;
				}
				u3 = (ro * u - T[now1] * (K->Potok[1] / Volume + ro * v * u / y - q2_1)) / ro3;
				v3 = (ro * v - T[now1] * (K->Potok[2] / Volume + ro * v * v / y - q2_2)) / ro3;
				p3 = (((p / (ggg - 1) + ro * (u * u + v * v) * 0.5) - T[now1] * (K->Potok[3] / Volume + //
					+v * (ggg * p / (ggg - 1) + ro * (u * u + v * v) * 0.5) / y - q3)) - //
					0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
				if (p3 <= 0)
				{
					p3 = 0.000001;
				}

				K->par[now2].ro = ro3;
				K->par[now2].Q = Q33;
				K->par[now2].p = p3;
				K->par[now2].u = u3;
				K->par[now2].v = v3;
			}
			else
			{
				K->par[now2].ro = K->par[now1].ro;
				K->par[now2].Q = K->par[now1].Q;
				K->par[now2].p = K->par[now1].p;
				K->par[now2].u = K->par[now1].u;
				K->par[now2].v = K->par[now1].v;
			}
			
			if (par1.Q / par1.ro < 90)
			{
				u = par1.u * (chi_real / chi_);
				v = par1.v * (chi_real / chi_);
				ro = par1.ro / kv(chi_real / chi_);
				p = par1.p;
				Q = par1.Q / kv(chi_real / chi_);
			}

			int k1 = 0, k2 = 0, k3 = 0, k4 = 0;

			if (kappa < 90)
			{
				k3 = 0;
				k4 = 0;
				if (((kv(u) + kv(v)) / (ggg * p / ro) > 1.0) && ((radius <= 400.0)))
				{
					k1 = 1;
					k2 = 0;
				}
				else
				{
					k1 = 0;
					k2 = 1;
				}
			}
			else
			{
				k1 = 0;
				k2 = 0;
				if (p / ro > 1.8)   // ( ((y < 8.0)&&(p / ro > (y * (-0.0238) + 0.36)))||( (y >= 8.0)&&(p / ro > 0.17) ) )
				{
					k4 = 0;
					k3 = 1;
				}
				else
				{
					k3 = 0;
					k4 = 1;
				}
			}

			double S1, S2;
			double q1_H1, q1_H2, q1_H3, q1_H4;
			double q21_H1, q21_H2, q21_H3, q21_H4;
			double q22_H1, q22_H2, q22_H3, q22_H4;
			double q3_H1, q3_H2, q3_H3, q3_H4;

			S1 = nu_H1 + nu_H2 + nu_H3 + nu_H4;
			S2 = nu_H1 * ((kv(u) + kv(v)) / 2.0 + (U_H1 / U_M_H1) * (p / ro))//
				+ nu_H2 * ((kv(u) + kv(v)) / 2.0 + (U_H2 / U_M_H2) * (p / ro))//
				+ nu_H3 * ((kv(u) + kv(v)) / 2.0 + (U_H3 / U_M_H3) * (p / ro))//
				+ nu_H4 * ((kv(u) + kv(v)) / 2.0 + (U_H4 / U_M_H4) * (p / ro));

			q1_H1 = (n_H_LISM_ / Kn_) * (k1 * S1 - nu_H1);
			q1_H2 = (n_H_LISM_ / Kn_) * (k2 * S1 - nu_H2);
			q1_H3 = (n_H_LISM_ / Kn_) * (k3 * S1 - nu_H3);
			q1_H4 = (n_H_LISM_ / Kn_) * (k4 * S1 - nu_H4);

			q21_H1 = (n_H_LISM_ / Kn_) * (k1 * S1 * u - nu_H1 * u_H1);
			q21_H2 = (n_H_LISM_ / Kn_) * (k2 * S1 * u - nu_H2 * u_H2);
			q21_H3 = (n_H_LISM_ / Kn_) * (k3 * S1 * u - nu_H3 * u_H3);
			q21_H4 = (n_H_LISM_ / Kn_) * (k4 * S1 * u - nu_H4 * u_H4);

			q22_H1 = (n_H_LISM_ / Kn_) * (k1 * S1 * v - nu_H1 * v_H1);
			q22_H2 = (n_H_LISM_ / Kn_) * (k2 * S1 * v - nu_H2 * v_H2);
			q22_H3 = (n_H_LISM_ / Kn_) * (k3 * S1 * v - nu_H3 * v_H3);
			q22_H4 = (n_H_LISM_ / Kn_) * (k4 * S1 * v - nu_H4 * v_H4);

			q3_H1 = (n_H_LISM_ / Kn_) * (k1 * S2 - nu_H1 * ((kv(u_H1) + kv(v_H1)) / 2.0 + (U_H1 / U_M_H1) * 2.0 * (p_H1 / ro_H1)));
			q3_H2 = (n_H_LISM_ / Kn_) * (k2 * S2 - nu_H2 * ((kv(u_H2) + kv(v_H2)) / 2.0 + (U_H2 / U_M_H2) * 2.0 * (p_H2 / ro_H2)));
			q3_H3 = (n_H_LISM_ / Kn_) * (k3 * S2 - nu_H3 * ((kv(u_H3) + kv(v_H3)) / 2.0 + (U_H3 / U_M_H3) * 2.0 * (p_H3 / ro_H3)));
			q3_H4 = (n_H_LISM_ / Kn_) * (k4 * S2 - nu_H4 * ((kv(u_H4) + kv(v_H4)) / 2.0 + (U_H4 / U_M_H4) * 2.0 * (p_H4 / ro_H4)));


			if (K->type != C_centr)
			{
				ro3 = ro_H1 - T[now1] * (K->Potok_H1[0] / Volume + ro_H1 * v_H1 / y - q1_H1);
				if (ro3 <= 0)
				{
					printf("Problemsssss  ro H1 < 0!\n");
					ro3 = 0.0000001;
				}
				u3 = (ro_H1 * u_H1 - T[now1] * (K->Potok_H1[1] / Volume + ro_H1 * v_H1 * u_H1 / y - q21_H1)) / ro3;
				v3 = (ro_H1 * v_H1 - T[now1] * (K->Potok_H1[2] / Volume + ro_H1 * v_H1 * v_H1 / y - q22_H1)) / ro3;
				p3 = (((p_H1 / (ggg - 1) + ro_H1 * (u_H1 * u_H1 + v_H1 * v_H1) * 0.5) - T[now1] * (K->Potok_H1[3] / Volume + //
					+v_H1 * (ggg * p_H1 / (ggg - 1) + ro_H1 * (u_H1 * u_H1 + v_H1 * v_H1) * 0.5) / y - q3_H1)) - //
					0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
				if (p3 <= 0)
				{
					p3 = 0.000001;
				}

				K->par[now2].ro_H1 = ro3;
				K->par[now2].p_H1 = p3;
				K->par[now2].u_H1 = u3;
				K->par[now2].v_H1 = v3;
			}
			else
			{
				K->par[now2].ro_H1 = K->par[now1].ro_H1;
				K->par[now2].p_H1 = K->par[now1].p_H1;
				K->par[now2].u_H1 = K->par[now1].u_H1;
				K->par[now2].v_H1 = K->par[now1].v_H1;
			}


			ro3 = ro_H2 - T[now1] * (K->Potok_H2[0] / Volume + ro_H2 * v_H2 / y - q1_H2);
			if (ro3 <= 0)
			{
				printf("Problemsssss  ro H2 < 0!  %lf,  %lf\n", x, y);
				ro3 = 0.000001;
			}
			u3 = (ro_H2 * u_H2 - T[now1] * (K->Potok_H2[1] / Volume + ro_H2 * v_H2 * u_H2 / y - q21_H2)) / ro3;
			v3 = (ro_H2 * v_H2 - T[now1] * (K->Potok_H2[2] / Volume + ro_H2 * v_H2 * v_H2 / y - q22_H2)) / ro3;
			p3 = (((p_H2 / (ggg - 1) + ro_H2 * (u_H2 * u_H2 + v_H2 * v_H2) * 0.5) - T[now1] * (K->Potok_H2[3] / Volume + //
				+v_H2 * (ggg * p_H2 / (ggg - 1) + ro_H2 * (u_H2 * u_H2 + v_H2 * v_H2) * 0.5) / y - q3_H2)) - //
				0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
			if (p3 <= 0)
			{
				p3 = 0.000001;
			}

			K->par[now2].ro_H2 = ro3;
			K->par[now2].p_H2 = p3;
			K->par[now2].u_H2 = u3;
			K->par[now2].v_H2 = v3;

			
			ro3 = ro_H3 - T[now1] * (K->Potok_H3[0] / Volume + ro_H3 * v_H3 / y - q1_H3);
			if (ro3 <= 0.0)
			{
				printf("Problemsssss  ro H3 < 0!  %lf,  %lf\n", x, y);
				ro3 = 0.000001;
			}
			u3 = (ro_H3 * u_H3 - T[now1] * (K->Potok_H3[1] / Volume + ro_H3 * v_H3 * u_H3 / y - q21_H3)) / ro3;
			v3 = (ro_H3 * v_H3 - T[now1] * (K->Potok_H3[2] / Volume + ro_H3 * v_H3 * v_H3 / y - q22_H3)) / ro3;
			p3 = (((p_H3 / (ggg - 1) + ro_H3 * (u_H3 * u_H3 + v_H3 * v_H3) * 0.5) - T[now1] * (K->Potok_H3[3] / Volume + //
				+v_H3 * (ggg * p_H3 / (ggg - 1) + ro_H3 * (u_H3 * u_H3 + v_H3 * v_H3) * 0.5) / y - q3_H3)) - //
				0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
			if (p3 <= 0)
			{
				p3 = 0.000001;
			}

			K->par[now2].ro_H3 = ro3;
			K->par[now2].p_H3 = p3;
			K->par[now2].u_H3 = u3;
			K->par[now2].v_H3 = v3;

			
			ro3 = ro_H4 - T[now1] * (K->Potok_H4[0] / Volume + ro_H4 * v_H4 / y - q1_H4);
			if (ro3 <= 0)
			{
				printf("Problemsssss  ro H4 < 0!  %lf,  %lf\n", x, y);
				ro3 = 0.000001;
			}
			u3 = (ro_H4 * u_H4 - T[now1] * (K->Potok_H4[1] / Volume + ro_H4 * v_H4 * u_H4 / y - q21_H4)) / ro3;
			v3 = (ro_H4 * v_H4 - T[now1] * (K->Potok_H4[2] / Volume + ro_H4 * v_H4 * v_H4 / y - q22_H4)) / ro3;
			p3 = (((p_H4 / (ggg - 1) + ro_H4 * (u_H4 * u_H4 + v_H4 * v_H4) * 0.5) - T[now1] * (K->Potok_H4[3] / Volume + //
				+v_H4 * (ggg * p_H4 / (ggg - 1) + ro_H4 * (u_H4 * u_H4 + v_H4 * v_H4) * 0.5) / y - q3_H4)) - //
				0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
			if (p3 <= 0)
			{
				p3 = 0.0000001;
			}

			K->par[now2].ro_H4 = ro3;
			K->par[now2].p_H4 = p3;
			K->par[now2].u_H4 = u3;
			K->par[now2].v_H4 = v3;
			

		}
	}
}

void Setka::Go_stationary_inner_infty(int step)
{
	// Стационарный решатель для внутренних ячеек в случае, когда Kn -> infty (особые источники)
	cout << "Setka.cpp    " << "START " << step << endl;
	int now1 = 1;
	int now2 = 0;
	double T[2];
	T[0] = T[1] = 0.00000001;
	mutex mut;
	double x_min = 0.0, y_min = 0.0;

	for (int st = 0; st < step; st++)
	{
		if (st % 10000 == 0 && st > 0)
		{
			cout << "Setka.cpp    " << st << " " << T[now2] << " " << x_min << " " << y_min << endl;
		}
		now1 = (now1 + 1) % 2; // Какие параметры сейчас берём
		now2 = (now2 + 1) % 2; // Какие параметры сейчас меняем
		double time = 10000000;
		T[now2] = 100000000;

		//omp_set_num_threads(4);
#pragma omp parallel for
		for (int num_cell = 0; num_cell < this->All_Cells_Inner.size(); num_cell++)
		{
			auto K = this->All_Cells_Inner[num_cell];
			if (K->type == C_centr)
			{
				continue;
			}
			K->Potok[0] = K->Potok[1] = K->Potok[2] = K->Potok[3] = K->Potok[4] = 0.0;

			double P[4];
			double PQ;
			Parametr par1 = K->par[now1];
			double n1, n2;
			double dist;
			double x, y;
			K->Get_Center(x, y);

			double radius = sqrt(kv(x) + kv(y));
			double Volume = K->Get_Volume();

			for (auto& i : K->Grans)
			{
				double x2, y2;
				double S = i->Get_square();
				i->Get_Center(x2, y2);
				dist = sqrt(kv(x - x2) + kv(y - y2));
				Parametr par2, par11, par3;
				par11 = par1;
				i->Get_par(par2, now1);
				i->Get_normal(n1, n2);

				double Vc, Vl, Vp;
				int met = 1;
				double dis = sqrt(kv(x2) + kv(y2));

				if (i->type == Usualy)
				{
					double x3, y3;
					i->Sosed->Get_Center(x3, y3);
					double dd = sqrt(kv(x3) + kv(y3));
					par11.ro = par11.ro * kv(radius) / kv(dis);
					par2.ro = par2.ro * kv(dd) / kv(dis);
					par11.Q = par11.Q * kv(radius) / kv(dis);
					par2.Q = par2.Q * kv(dd) / kv(dis);
					par11.p = par11.p * pow(radius / dis, 2 * ggg);
					par2.p = par2.p * pow(dd / dis, 2 * ggg);
					polar_perenos(x, y, x2, y2, par11.u, par11.v);
					polar_perenos(x3, y3, x2, y2, par2.u, par2.v);

					if (K->type == C_1 && i->Sosed->type == C_1)
					{
						Parametr pp1, pp2;
						i->Get_par_TVD_radial(pp1, now1);
						i->Gran_copy->Get_par_TVD_radial(pp2, now1);

						par11.ro = pp1.ro;
						par2.ro = pp2.ro;
						par11.p = pp1.p;
						par2.p = pp2.p;
					}
				}
				else if (i->type == Inner_sphere)
				{
					par11.ro = par2.ro;
					par11.Q = par2.Q;
					par11.p = par2.p;
					par11.u = par2.u;
					par11.v = par2.v;
				}
				else if (i->type == Axis)
				{
					par11.ro = par11.ro * kv(radius) / kv(dis);
					par11.Q = par11.Q * kv(radius) / kv(dis);
					par11.p = par11.p * pow(radius / dis, 2 * ggg);
					polar_perenos(x, y, x2, y2, par11.u, par11.v);
					par2 = par11;
					par2.v = -par2.v;
				}


				if (K->type != C_centr)
				{
					time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro, par11.Q, par11.p, par11.u, par11.v, par2.ro, par2.Q, //
						par2.p, par2.u, par2.v, 0.0, P, PQ, n1, n2, dist, met, Vl, Vc, Vp));
					for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
					{
						K->Potok[k] = K->Potok[k] + P[k] * S;
					}
					K->Potok[4] = K->Potok[4] + PQ * S;
				}
			}

			if (time < T[now2])
			{
				mut.lock();
				T[now2] = time;
				x_min = x;
				y_min = y;
				mut.unlock();
			}



			double U_M_H4;
			double U_H4;
			double sigma_H4;
			double nu_H4;
			double q2_1, q2_2, q3;
			double u, v, ro, p, Q;

			u = par1.u * (chi_real / chi_);
			v = par1.v *(chi_real / chi_);
			ro = par1.ro / kv(chi_real / chi_);
			p = par1.p;
			Q = par1.Q / kv(chi_real / chi_);


			double u_H4 = Velosity_inf, v_H4 = 0.0, ro_H4 = 1.0, p_H4 = 0.5;

			U_M_H4 = sqrt(kv(u - u_H4) + kv(v - v_H4) + (64.0 / (9.0 * pi_)) //
				* (p / ro + 2.0 * p_H4 / ro_H4));

			U_H4 = sqrt(kv(u - u_H4) + kv(v - v_H4) + (4.0 / pi_) //
				* (p / ro + 2.0 * p_H4 / ro_H4));

			sigma_H4 = kv(1.0 - a_2 * log(U_M_H4));

			nu_H4 = ro * ro_H4 * U_M_H4 * sigma_H4;

			q2_1 = (n_p_LISM_ / Kn_) * (nu_H4 * (u_H4 - u));
			q2_2 = (n_p_LISM_ / Kn_) * (nu_H4 * (v_H4 - v));
			q3 = (n_p_LISM_ / Kn_) * (nu_H4 * ((kv(u_H4) + kv(v_H4) - kv(u) - kv(v)) / 2.0 + //
					(U_H4 / U_M_H4) * (2.0 * p_H4 / ro_H4 - p / ro))) / (chi_real / chi_);


			double ro2_p, p2_p, V1_p, V2_p, QQ2;
			double ro3, p3, u3, v3, Q33;
			double kappa = Q / ro;

			u = par1.u;
			v = par1.v;
			ro = par1.ro;
			p = par1.p;
			Q = par1.Q;

			if (K->type != C_centr)
			{
				ro3 = ro - T[now1] * (K->Potok[0] / Volume + ro * v / y);
				Q33 = Q - (T[now1] / Volume) * K->Potok[4] - T[now1] * Q * v / y;
				if (ro3 <= 0)
				{
					printf("Problemsssss  ro < 0!   %lf,  %lf\n", x, y);
					ro3 = 0.00001;
				}
				u3 = (ro * u - T[now1] * (K->Potok[1] / Volume + ro * v * u / y - q2_1)) / ro3;
				v3 = (ro * v - T[now1] * (K->Potok[2] / Volume + ro * v * v / y - q2_2)) / ro3;
				p3 = (((p / (ggg - 1) + ro * (u * u + v * v) * 0.5) - T[now1] * (K->Potok[3] / Volume + //
					+v * (ggg * p / (ggg - 1) + ro * (u * u + v * v) * 0.5) / y - q3)) - //
					0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
				if (p3 <= 0)
				{
					p3 = 0.000001;
				}

				K->par[now2].ro = ro3;
				K->par[now2].Q = Q33;
				K->par[now2].p = p3;
				K->par[now2].u = u3;
				K->par[now2].v = v3;
			}
			else
			{
				K->par[now2].ro = K->par[now1].ro;
				K->par[now2].Q = K->par[now1].Q;
				K->par[now2].p = K->par[now1].p;
				K->par[now2].u = K->par[now1].u;
				K->par[now2].v = K->par[now1].v;
			}
		}
	}
}

void Setka::Go_5_komponent_infty(int step, bool dvig)
{
	// Нестационарный решатель в случае, когда Kn -> infty (особые источники)

	//cout << "Setka.cpp    " << "START " << step << endl;
	int now1 = 1;
	int now2 = 0;
	double T[2];
	T[0] = T[1] = 0.00000001;
	mutex mut;

	for (int st = 0; st < step; st++)
	{
		//cout << "Setka.cpp    " << st << endl;
		if (st % 10000 == 0 && st > 0)
		{
			cout << "Setka.cpp    " << st << " " << T[now2] << endl;
		}
		now1 = (now1 + 1) % 2; // Какие параметры сейчас берём
		now2 = (now2 + 1) % 2; // Какие параметры сейчас меняем
		double time = 10000000;
		T[now2] = 100000000;

		if (dvig)
		{
			this->Move_surface(now1, T[now1]);
			this->Move_Setka_Calculate_2(T[now1]);
		}

		//omp_set_num_threads(4);
#pragma omp parallel for
		for (int num_cell = 0; num_cell < this->All_Cells.size(); num_cell++)
		{
			auto K = this->All_Cells[num_cell];
			K->Potok[0] = K->Potok[1] = K->Potok[2] = K->Potok[3] = K->Potok[4] = 0.0;

			double P[4];
			double PQ;
			Parametr par1 = K->par[now1];
			double n1, n2;
			double dist;
			double x, y;
			K->Get_Center(x, y);

			//y = K->Get_Volume_rotate(5.0) / (2.0 * pi_ * K->Get_Volume() / 72.0);

			double radius = sqrt(kv(x) + kv(y));
			if (radius < R11_)
			{
				continue;
			}
			double Volume = K->Get_Volume(); // K->Get_Volume_rotate(alpha_rot);
			//cout << "Setka.cpp    " << "Vol = " << Volume << endl;
			double Volume2 = K->Get_Volume_posle(); // K->Get_Volume_posle_rotate(alpha_rot);
			//cout << "Setka.cpp    " << "Vol2 = " << Volume2 << endl;
			double W = 0.0;
			bool np = true;

			for (auto& i : K->Grans)
			{
				double x2, y2, x4, y4;
				//double S = i->Get_square_rotate(alpha_rot); 
				double S = i->Get_square();
				//cout << "Setka.cpp    " << "S = " << S << endl;
				i->Get_Center(x2, y2);
				i->Get_Center_posle(x4, y4);
				dist = sqrt(kv(x - x2) + kv(y - y2));
				Parametr par2, par11;
				par11 = par1;
				i->Get_par(par2, now1);
				i->Get_normal(n1, n2);
				double Vc, Vl, Vp;
				W = ((x4 - x2) * n1 + (y4 - y2) * n2) / T[now1];
				int met = 1;  // 1
				double dis = sqrt(kv(x2) + kv(y2));

				if (i->type == Usualy)
				{
					double x3, y3;
					i->Sosed->Get_Center(x3, y3);
					double dd = sqrt(kv(x3) + kv(y3));

					if (K->type == C_1)
					{
						par11.ro = par11.ro * kv(radius) / kv(dis);
						par11.Q = par11.Q * kv(radius) / kv(dis);
						par11.p = par11.p * pow(radius / dis, 2 * ggg);
						polar_perenos(x, y, x2, y2, par11.u, par11.v);

						Parametr pp1;
						i->Get_par_TVD_radial(pp1, now1);
						par11.ro = pp1.ro;
						par11.p = pp1.p;
					}
					else
					{
						i->Get_par_TVD(par11, now1);
					}

					if (i->Sosed->type == C_1)
					{
						par2.ro = par2.ro * kv(dd) / kv(dis);
						par2.Q = par2.Q * kv(dd) / kv(dis);
						par2.p = par2.p * pow(dd / dis, 2 * ggg);
						polar_perenos(x3, y3, x2, y2, par2.u, par2.v);

						Parametr pp2;
						i->Gran_copy->Get_par_TVD_radial(pp2, now1);
						par2.ro = pp2.ro;
						par2.p = pp2.p;
					}
					else
					{
						i->Gran_copy->Get_par_TVD(par2, now1);
					}
				}
				else if (i->type == Inner_sphere)
				{
					par11.ro = par11.ro * kv(radius) / kv(dis);
					par11.Q = par11.Q * kv(radius) / kv(dis);
					par11.p = par11.p * pow(radius / dis, 2 * ggg);
				}
				else if (K->type == C_1 && i->type == Axis)
				{
					par11.ro = par11.ro * kv(radius) / kv(dis);
					par11.Q = par11.Q * kv(radius) / kv(dis);
					par11.p = par11.p * pow(radius / dis, 2 * ggg);
					polar_perenos(x, y, x2, y2, par11.u, par11.v);
					//par11.v = 0.0;
					par2 = par11;
					par2.v = -par2.v;
				}
				else if (K->type != C_1 && i->type == Axis)
				{
					i->Get_par_TVD(par11, now1);
					par2 = par11;
					par2.v = -par2.v;
				}


				/*if ((i->A->type == P_Contact && i->B->type != P_Contact)||(i->A->type != P_Contact && i->B->type == P_Contact))
				{
					met = 0;
				}*/

				//np = false;
				/*bool god = false;

				if ((i->A->type == P_Contact && i->B->type == P_Contact))
				{
					god = true;
				}
				else if ((i->A->type == P_Inner_shock && i->B->type == P_Inner_shock))
				{
					god = true;
				}
				else if ((i->A->type == P_Outer_shock && i->B->type == P_Outer_shock))
				{
					god = true;
				}


				if ((i->A->type == P_Contact && i->B->type != P_Contact))
				{
					np = true;
				}
				else if ((i->A->type != P_Contact && i->B->type == P_Contact))
				{
					np = true;
				}
				else if ((i->A->type == P_Inner_shock && i->B->type != P_Inner_shock))
				{
					np = true;
				}
				else if ((i->A->type != P_Inner_shock && i->B->type == P_Inner_shock))
				{
					np = true;
				}
				else if ((i->A->type == P_Outer_shock && i->B->type != P_Outer_shock))
				{
					np = true;
				}
				else if ((i->A->type != P_Outer_shock && i->B->type == P_Outer_shock))
				{
					np = true;
				}*/

				np = true;  // Сглаживание
				if ((i->A->type == P_Contact && i->B->type == P_Contact))
				{
					np = false;
					met = 1;
				}

				//np = false;

				/*if (y < 30 && x > 120 && x < 160)
				{
					np = true;
				}*/


				if (false)//(god == true)
				{
					double Vl, Vp, VV;
					vector<double> qqq1(5);
					vector<double> qqq2(5);
					vector<double> qqq(5);
					vector<double> n(3);
					n[0] = n1;
					n[1] = n2;
					n[2] = 0.0;

					qqq1[0] = par11.ro;
					qqq1[1] = par11.p;
					qqq1[2] = par11.u;
					qqq1[3] = par11.v;
					qqq1[4] = 0.0;

					qqq2[0] = par2.ro;
					qqq2[1] = par2.p;
					qqq2[2] = par2.u;
					qqq2[3] = par2.v;
					qqq2[4] = 0.0;

					auto SS = Solvers();
					time = min(time, SS.Godunov_Solver_Alexashov(qqq1, qqq2, n, qqq, Vl, Vp, VV, W, dist));

					K->Potok[0] = K->Potok[0] + qqq[0] * S;
					K->Potok[1] = K->Potok[1] + qqq[1] * S;
					K->Potok[2] = K->Potok[2] + qqq[2] * S;
					K->Potok[3] = K->Potok[3] + qqq[4] * S;

					time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro, par11.Q, par11.p, par11.u, par11.v, par2.ro, par2.Q, //
						par2.p, par2.u, par2.v, W, P, PQ, n1, n2, dist, met, Vl, Vc, Vp, np));

					K->Potok[4] = K->Potok[4] + PQ * S;
				}
				else
				{
					time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro, par11.Q, par11.p, par11.u, par11.v, par2.ro, par2.Q, //
						par2.p, par2.u, par2.v, W, P, PQ, n1, n2, dist, met, Vl, Vc, Vp, np));
					for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
					{
						K->Potok[k] = K->Potok[k] + P[k] * S;
					}
					K->Potok[4] = K->Potok[4] + PQ * S;
				}

			}

			if (time < T[now2])
			{
				mut.lock();
				T[now2] = time;
				mut.unlock();
			}



			double U_M_H4;
			double U_H4;
			double sigma_H4;
			double nu_H4;
			double q2_1, q2_2, q3;
			double u, v, ro, p, Q;
			if (par1.Q / par1.ro < 90)//(K->type == C_centr || K->type == C_1 || K->type == C_2 || K->type == C_3)//(par1.Q / par1.ro < 90)//
			{
				u = par1.u * (chi_real / chi_);
				v = par1.v * (chi_real / chi_);
				ro = par1.ro / kv(chi_real / chi_);
				p = par1.p;
				Q = par1.Q / kv(chi_real / chi_);
			}
			else
			{
				u = par1.u;
				v = par1.v;
				ro = par1.ro;
				p = par1.p;
				Q = par1.Q;
			}

			double u_H4 = par1.u_H4, v_H4 = par1.v_H4, ro_H4 = par1.ro_H4, p_H4 = par1.p_H4;

			U_M_H4 = sqrt(kv(u - u_H4) + kv(v - v_H4) + (64.0 / (9.0 * pi)) //
				* (p / ro + 2.0 * p_H4 / ro_H4));

			U_H4 = sqrt(kv(u - u_H4) + kv(v - v_H4) + (4.0 / pi) //
				* (p / ro + 2.0 * p_H4 / ro_H4));

			sigma_H4 = kv(1.0 - a_2 * log(U_M_H4));

			nu_H4 = ro * ro_H4 * U_M_H4 * sigma_H4;

			q2_1 = (n_p_LISM_ / Kn_) * (nu_H4 * (u_H4 - u));
			q2_2 = (n_p_LISM_ / Kn_) * (nu_H4 * (v_H4 - v));

			if (par1.Q / par1.ro < 90)//(K->type == C_centr || K->type == C_1 || K->type == C_2 || K->type == C_3)//(par1.Q / par1.ro < 90)//
			{
				q3 = (n_p_LISM_ / Kn_) * (nu_H4 * ((kv(u_H4) + kv(v_H4) - kv(u) - kv(v)) / 2.0 + //
						(U_H4 / U_M_H4) * (2.0 * p_H4 / ro_H4 - p / ro))) / (chi_real / chi_);
			}
			else
			{
				q3 = (n_p_LISM_ / Kn_) * (nu_H4 * ((kv(u_H4) + kv(v_H4) - kv(u) - kv(v)) / 2.0 + //
						(U_H4 / U_M_H4) * (2.0 * p_H4 / ro_H4 - p / ro)));
			}


			/*q2_1 = 0.0;
			q2_2 = 0.0;
			q3 = 0.0;*/

			double ro2_p, p2_p, V1_p, V2_p, QQ2;
			double ro3, p3, u3, v3, Q33;
			double kappa = Q / ro;


			u = par1.u;
			v = par1.v;
			ro = par1.ro;
			p = par1.p;
			Q = par1.Q;

			if (K->type != C_centr)
			{
				//ro3 = ro * Volume / Volume2 - T[now1] * (K->Potok[0] / Volume2);  // В декартовых
				ro3 = ro * Volume / Volume2 - T[now1] * (K->Potok[0] / Volume2 + ro * v / y);  // В цилиндрических

				//Q33 = Q * Volume / Volume2 - (T[now1] / Volume2) * K->Potok[4];
				Q33 = Q * Volume / Volume2 - (T[now1] / Volume2) * K->Potok[4] - T[now1] * Q * v / y;
				if (ro3 <= 0.0)
				{
					printf("Problemsssss  ro < 0!    %lf, %lf\n", x, y);
					ro3 = 0.00001;
				}
				//u3 = (ro * u * Volume / Volume2 - T[now1] * (K->Potok[1] / Volume2  - q2_1)) / ro3;
				//v3 = (ro * v * Volume / Volume2 - T[now1] * ( (K->Potok[2] - 2.0 * p * sin(0.5 * alpha_rot * pi_ / 180.0) * K->Get_Volume()) / Volume2  - q2_2)) / ro3;
				//p3 = (((p / (ggg - 1) + ro * (u * u + v * v) * 0.5) * Volume / Volume2 - T[now1] * (K->Potok[3] / Volume2 - q3)) - //
				//	0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
				u3 = (ro * u * Volume / Volume2 - T[now1] * (K->Potok[1] / Volume2 + ro * v * u / y - q2_1)) / ro3;
				v3 = (ro * v * Volume / Volume2 - T[now1] * (K->Potok[2] / Volume2 + ro * v * v / y - q2_2)) / ro3;
				p3 = (((p / (ggg - 1) + ro * (u * u + v * v) * 0.5) * Volume / Volume2 - T[now1] * (K->Potok[3] / Volume2 + //
					+v * (ggg * p / (ggg - 1) + ro * (u * u + v * v) * 0.5) / y - q3)) - //
					0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
				if (p3 <= 0)
				{
					p3 = 0.000001;
				}

				K->par[now2].ro = ro3;
				K->par[now2].Q = Q33;
				K->par[now2].p = p3;
				K->par[now2].u = u3;
				K->par[now2].v = v3;
			}
			else
			{
				K->par[now2].ro = K->par[now1].ro;
				K->par[now2].Q = K->par[now1].Q;
				K->par[now2].p = K->par[now1].p;
				K->par[now2].u = K->par[now1].u;
				K->par[now2].v = K->par[now1].v;
			}


		}


		for (auto& i : this->All_Points)
		{
			i->x = i->x2;
			i->y = i->y2;
			i->Vx = 0.0;
			i->Vy = 0.0;
			i->count = 0;
		}
	}
}

void Setka::Go_stationary_5_komponent_inner(int step)
{
	cout << "Setka.cpp    " << "START " << step << endl;
	int now1 = 1;
	int now2 = 0;
	double T[2];
	T[0] = T[1] = 0.00000001;
	mutex mut;
	double x_min = 0.0, y_min = 0.0;

	for (int st = 0; st < step; st++)
	{
		if (st % 10000 == 0 && st > 0)
		{
			cout << "Setka.cpp    " << st << " " << T[now2] << " " << x_min << " " << y_min << endl;
		}
		now1 = (now1 + 1) % 2; // Какие параметры сейчас берём
		now2 = (now2 + 1) % 2; // Какие параметры сейчас меняем
		double time = 10000000;
		T[now2] = 100000000;

		//omp_set_num_threads(4);
#pragma omp parallel for
		for (int num_cell = 0; num_cell < this->All_Cells_Inner.size(); num_cell++)
		{
			auto K = this->All_Cells_Inner[num_cell];
			if (K->type == C_centr)
			{
				continue;
			}
			K->Potok[0] = K->Potok[1] = K->Potok[2] = K->Potok[3] = K->Potok[4] = 0.0;
			K->Potok_H1[0] = K->Potok_H1[1] = K->Potok_H1[2] = K->Potok_H1[3] = 0.0;
			K->Potok_H2[0] = K->Potok_H2[1] = K->Potok_H2[2] = K->Potok_H2[3] = 0.0;
			K->Potok_H3[0] = K->Potok_H3[1] = K->Potok_H3[2] = K->Potok_H3[3] = 0.0;
			K->Potok_H4[0] = K->Potok_H4[1] = K->Potok_H4[2] = K->Potok_H4[3] = 0.0;

			double P[4];
			double PQ;
			Parametr par1 = K->par[now1];
			double n1, n2;
			double dist;
			double x, y;
			K->Get_Center(x, y);

			//cout << "Setka.cpp    " << "Do = " << y << endl;
			//y = K->Get_Volume_rotate(360.0) / (2.0 * pi_ * K->Get_Volume());
			//y = K->Get_Volume_rotate(5.0) / (2.0 * pi_ * K->Get_Volume() / 72.0);
			//cout << "Setka.cpp    " << "Posle = " << y << endl;

			double radius = sqrt(kv(x) + kv(y));
			double Volume = K->Get_Volume();

			for (auto& i : K->Grans)
			{
				double x2, y2;
				double S = i->Get_square();
				i->Get_Center(x2, y2);
				dist = sqrt(kv(x - x2) + kv(y - y2));
				Parametr par2, par11, par3;
				par11 = par1;
				i->Get_par(par2, now1);
				i->Get_normal(n1, n2);

				double Vc, Vl, Vp;
				int met = 1;
				double dis = sqrt(kv(x2) + kv(y2));

				if (i->type == Usualy)
				{
					double x3, y3;
					i->Sosed->Get_Center(x3, y3);
					double dd = sqrt(kv(x3) + kv(y3));
					par11.ro = par11.ro * kv(radius) / kv(dis);
					par2.ro = par2.ro * kv(dd) / kv(dis);
					par11.ro_H1 = par11.ro_H1 * pow(radius / dis, H_pow);
					par2.ro_H1 = par2.ro_H1 * pow(dd / dis, H_pow);
					//par11.p_H1 = par11.p_H1 * pow(radius / dis, 2 * ggg);
					//par2.p_H1 = par2.p_H1 * pow(dd / dis, 2 * ggg);
					par11.Q = par11.Q * kv(radius) / kv(dis);
					par2.Q = par2.Q * kv(dd) / kv(dis);
					par11.p = par11.p * pow(radius / dis, 2 * ggg);
					par2.p = par2.p * pow(dd / dis, 2 * ggg);
					polar_perenos(x, y, x2, y2, par11.u, par11.v);
					polar_perenos(x3, y3, x2, y2, par2.u, par2.v);
					polar_perenos(x, y, x2, y2, par11.u_H1, par11.v_H1);
					polar_perenos(x3, y3, x2, y2, par2.u_H1, par2.v_H1);

					if (K->type == C_1 && i->Sosed->type == C_1)
					{
						Parametr pp1, pp2;
						i->Get_par_TVD_radial(pp1, now1);
						i->Gran_copy->Get_par_TVD_radial(pp2, now1);
						//par11.u_H1 = pp1.u_H1;
						//par11.v_H1 = pp1.v_H1;
						//par2.u_H1 = pp2.u_H1;
						//par2.v_H1 = pp2.v_H1;
						par11.p_H1 = pp1.p_H1;
						par2.p_H1 = pp2.p_H1;
						par11.ro_H1 = pp1.ro_H1;
						par2.ro_H1 = pp2.ro_H1;

						par11.ro = pp1.ro;
						par2.ro = pp2.ro;
						par11.p = pp1.p;
						par2.p = pp2.p;
					}
				}
				else if (i->type == Inner_sphere)
				{
					/*par11.ro = par11.ro * kv(radius) / kv(dis);
					par11.ro_H1 = par11.ro_H1 * pow(radius / dis, H_pow);
					par11.Q = par11.Q * kv(radius) / kv(dis);
					par11.p = par11.p * pow(radius / dis, 2 * ggg);*/

					//cout << "Setka.cpp    " << "inner " << dis << endl;
					/*par2.ro_H3 = par1.ro_H3;
					par2.p_H3 = par1.p_H3;
					par2.u_H3 = par1.u_H3;
					par2.v_H3 = par1.v_H3;
					par2.ro_H4 = par1.ro_H4;
					par2.p_H4 = par1.p_H4;
					par2.u_H4 = par1.u_H4;
					par2.v_H4 = par1.v_H4;*/
					par11.ro = par2.ro;
					par11.Q = par2.Q;
					par11.p = par2.p;
					par11.u = par2.u;
					par11.v = par2.v;
					par11.ro_H1 = par2.ro_H1;
					par11.p_H1 = par2.p_H1;
					par11.u_H1 = par2.u_H1;
					par11.v_H1 = par2.v_H1;
				}
				else if (i->type == Axis)
				{
					par11.ro = par11.ro * kv(radius) / kv(dis);
					par11.Q = par11.Q * kv(radius) / kv(dis);
					par11.p = par11.p * pow(radius / dis, 2 * ggg);
					par11.ro_H1 = par11.ro_H1 * pow(radius / dis, H_pow);
					par11.p_H1 = par11.p_H1 * pow(radius / dis, 2 * ggg);
					polar_perenos(x, y, x2, y2, par11.u, par11.v);
					polar_perenos(x, y, x2, y2, par11.u_H1, par11.v_H1);
					par2 = par11;
					par2.v = -par2.v;
					par2.v_H1 = -par2.v_H1;
					par2.v_H2 = -par2.v_H2;
					par2.v_H3 = -par2.v_H3;
					par2.v_H4 = -par2.v_H4;
				}

				if (K->type != C_centr)
				{
					time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro_H1, par11.Q, par11.p_H1, par11.u_H1, par11.v_H1, par2.ro_H1, par2.Q, //
						par2.p_H1, par2.u_H1, par2.v_H1, 0.0, P, PQ, n1, n2, dist, met, Vl, Vc, Vp));
					for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
					{
						K->Potok_H1[k] = K->Potok_H1[k] + P[k] * S;
					}
				}

				time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro_H2, par11.Q, par11.p_H2, par11.u_H2, par11.v_H2, par2.ro_H2, par2.Q, //
					par2.p_H2, par2.u_H2, par2.v_H2, 0.0, P, PQ, n1, n2, dist, met, Vl, Vc, Vp));
				for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
				{
					K->Potok_H2[k] = K->Potok_H2[k] + P[k] * S;
				}

				time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro_H3, par11.Q, par11.p_H3, par11.u_H3, par11.v_H3, par2.ro_H3, par2.Q, //
					par2.p_H3, par2.u_H3, par2.v_H3, 0.0, P, PQ, n1, n2, dist, met, Vl, Vc, Vp));
				for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
				{
					K->Potok_H3[k] = K->Potok_H3[k] + P[k] * S;
				}

				time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro_H4, par11.Q, par11.p_H4, par11.u_H4, par11.v_H4, par2.ro_H4, par2.Q, //
					par2.p_H4, par2.u_H4, par2.v_H4, 0.0, P, PQ, n1, n2, dist, met, Vl, Vc, Vp));
				for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
				{
					K->Potok_H4[k] = K->Potok_H4[k] + P[k] * S;
				}


				if (K->type != C_centr)
				{
					/*double Vl, Vp, VV;
					vector<double> qqq1(5);
					vector<double> qqq2(5);
					vector<double> qqq(5);
					vector<double> n(3);
					n[0] = n1;
					n[1] = n2;
					n[2] = 0.0;

					qqq1[0] = par11.ro;
					qqq1[1] = par11.p;
					qqq1[2] = par11.u;
					qqq1[3] = par11.v;
					qqq1[4] = 0.0;

					qqq2[0] = par2.ro;
					qqq2[1] = par2.p;
					qqq2[2] = par2.u;
					qqq2[3] = par2.v;
					qqq2[4] = 0.0;

					auto SS = Solvers();
					time = min(time, SS.Godunov_Solver_Alexashov(qqq1, qqq2, n, qqq, Vl, Vp, VV, 0.0, dist));

					K->Potok[0] = K->Potok[0] + qqq[0] * S;
					K->Potok[1] = K->Potok[1] + qqq[1] * S;
					K->Potok[2] = K->Potok[2] + qqq[2] * S;
					K->Potok[3] = K->Potok[3] + qqq[4] * S;*/


					time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro, par11.Q, par11.p, par11.u, par11.v, par2.ro, par2.Q, //
						par2.p, par2.u, par2.v, 0.0, P, PQ, n1, n2, dist, met, Vl, Vc, Vp));
					for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
					{
						K->Potok[k] = K->Potok[k] + P[k] * S;
					}
					K->Potok[4] = K->Potok[4] + PQ * S;
				}
			}

			if (time < T[now2])
			{
				mut.lock();
				T[now2] = time;
				x_min = x;
				y_min = y;
				mut.unlock();
			}



			double U_M_H1, U_M_H2, U_M_H3, U_M_H4;
			double U_H1, U_H2, U_H3, U_H4;
			double sigma_H1, sigma_H2, sigma_H3, sigma_H4;
			double nu_H1, nu_H2, nu_H3, nu_H4;
			double q2_1, q2_2, q3;
			double u, v, ro, p, Q;

			u = par1.u;// *(chi_real / chi_);
			v = par1.v;// *(chi_real / chi_);
			ro = par1.ro;// / kv(chi_real / chi_);
			p = par1.p;
			Q = par1.Q;// / kv(chi_real / chi_);

			double u_H1 = par1.u_H1, v_H1 = par1.v_H1, ro_H1 = par1.ro_H1, p_H1 = par1.p_H1;
			double u_H2 = par1.u_H2, v_H2 = par1.v_H2, ro_H2 = par1.ro_H2, p_H2 = par1.p_H2;
			double u_H3 = par1.u_H3, v_H3 = par1.v_H3, ro_H3 = par1.ro_H3, p_H3 = par1.p_H3;
			double u_H4 = par1.u_H4, v_H4 = par1.v_H4, ro_H4 = par1.ro_H4, p_H4 = par1.p_H4;

			U_M_H1 = sqrt(kv(u - u_H1) + kv(v - v_H1) + (64.0 / (9.0 * pi_)) //
				* (p / ro + 2.0 * p_H1 / ro_H1));
			U_M_H2 = sqrt(kv(u - u_H2) + kv(v - v_H2) + (64.0 / (9.0 * pi_)) //
				* (p / ro + 2.0 * p_H2 / ro_H2));
			U_M_H3 = sqrt(kv(u - u_H3) + kv(v - v_H3) + (64.0 / (9.0 * pi_)) //
				* (p / ro + 2.0 * p_H3 / ro_H3));
			U_M_H4 = sqrt(kv(u - u_H4) + kv(v - v_H4) + (64.0 / (9.0 * pi_)) //
				* (p / ro + 2.0 * p_H4 / ro_H4));

			U_H1 = sqrt(kv(u - u_H1) + kv(v - v_H1) + (4.0 / pi_) //
				* (p / ro + 2.0 * p_H1 / ro_H1));
			U_H2 = sqrt(kv(u - u_H2) + kv(v - v_H2) + (4.0 / pi_) //
				* (p / ro + 2.0 * p_H2 / ro_H2));
			U_H3 = sqrt(kv(u - u_H3) + kv(v - v_H3) + (4.0 / pi_) //
				* (p / ro + 2.0 * p_H3 / ro_H3));
			U_H4 = sqrt(kv(u - u_H4) + kv(v - v_H4) + (4.0 / pi_) //
				* (p / ro + 2.0 * p_H4 / ro_H4));

			sigma_H1 = kv(1.0 - a_2 * log(U_M_H1)); // 0.1243
			sigma_H2 = kv(1.0 - a_2 * log(U_M_H2));
			sigma_H3 = kv(1.0 - a_2 * log(U_M_H3)); // 0.1121     a_2
			sigma_H4 = kv(1.0 - a_2 * log(U_M_H4));

			nu_H1 = ro * ro_H1 * U_M_H1 * sigma_H1;
			nu_H2 = ro * ro_H2 * U_M_H2 * sigma_H2;
			nu_H3 = ro * ro_H3 * U_M_H3 * sigma_H3;
			nu_H4 = ro * ro_H4 * U_M_H4 * sigma_H4;

			q2_1 = (n_p_LISM_ / Kn_) * (nu_H1 * (u_H1 - u) + nu_H2 * (u_H2 - u) //
				+ nu_H3 * (u_H3 - u) + nu_H4 * (u_H4 - u));
			q2_2 = (n_p_LISM_ / Kn_) * (nu_H1 * (v_H1 - v) + nu_H2 * (v_H2 - v) //
				+ nu_H3 * (v_H3 - v) + nu_H4 * (v_H4 - v));
			q3 = (n_p_LISM_ / Kn_) * (nu_H1 * ((kv(u_H1) + kv(v_H1) - kv(u) - kv(v)) / 2.0 + //
				(U_H1 / U_M_H1) * (2.0 * p_H1 / ro_H1 - p / ro)) + //
				nu_H2 * ((kv(u_H2) + kv(v_H2) - kv(u) - kv(v)) / 2.0 + //
					(U_H2 / U_M_H2) * (2.0 * p_H2 / ro_H2 - p / ro)) + //
				nu_H3 * ((kv(u_H3) + kv(v_H3) - kv(u) - kv(v)) / 2.0 + //
					(U_H3 / U_M_H3) * (2.0 * p_H3 / ro_H3 - p / ro)) + //
				nu_H4 * ((kv(u_H4) + kv(v_H4) - kv(u) - kv(v)) / 2.0 + //
					(U_H4 / U_M_H4) * (2.0 * p_H4 / ro_H4 - p / ro)));// / (chi_real / chi_);


			/*q2_1 = 0.0;
			q2_2 = 0.0;
			q3 = 0.0;*/

			double ro2_p, p2_p, V1_p, V2_p, QQ2;
			double ro3, p3, u3, v3, Q33;
			double kappa = Q / ro;

			u = par1.u;
			v = par1.v;
			ro = par1.ro;
			p = par1.p;
			Q = par1.Q;

			if (K->type != C_centr)
			{
				ro3 = ro - T[now1] * (K->Potok[0] / Volume + ro * v / y);
				Q33 = Q - (T[now1] / Volume) * K->Potok[4] - T[now1] * Q * v / y;
				if (ro3 <= 0)
				{
					printf("Problemsssss  ro < 0!   %lf,  %lf\n", x, y);
					ro3 = 0.00001;
				}
				u3 = (ro * u - T[now1] * (K->Potok[1] / Volume + ro * v * u / y - q2_1)) / ro3;
				v3 = (ro * v - T[now1] * (K->Potok[2] / Volume + ro * v * v / y - q2_2)) / ro3;
				p3 = (((p / (ggg - 1) + ro * (u * u + v * v) * 0.5) - T[now1] * (K->Potok[3] / Volume + //
					+v * (ggg * p / (ggg - 1) + ro * (u * u + v * v) * 0.5) / y - q3)) - //
					0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
				if (p3 <= 0)
				{
					p3 = 0.000001;
				}

				K->par[now2].ro = ro3;
				K->par[now2].Q = Q33;
				K->par[now2].p = p3;
				K->par[now2].u = u3;
				K->par[now2].v = v3;
			}
			else
			{
				K->par[now2].ro = K->par[now1].ro;
				K->par[now2].Q = K->par[now1].Q;
				K->par[now2].p = K->par[now1].p;
				K->par[now2].u = K->par[now1].u;
				K->par[now2].v = K->par[now1].v;
			}

			u = par1.u;// *(chi_real / chi_);
			v = par1.v;// *(chi_real / chi_);
			ro = par1.ro;// / kv(chi_real / chi_);
			p = par1.p;
			Q = par1.Q;// / kv(chi_real / chi_);
			

			int k1 = 1, k2 = 0, k3 = 0, k4 = 0;


			double S1, S2;
			double q1_H1, q1_H2, q1_H3, q1_H4;
			double q21_H1, q21_H2, q21_H3, q21_H4;
			double q22_H1, q22_H2, q22_H3, q22_H4;
			double q3_H1, q3_H2, q3_H3, q3_H4;

			S1 = nu_H1 + nu_H2 + nu_H3 + nu_H4;
			S2 = nu_H1 * ((kv(u) + kv(v)) / 2.0 + (U_H1 / U_M_H1) * (p / ro))//
				+ nu_H2 * ((kv(u) + kv(v)) / 2.0 + (U_H2 / U_M_H2) * (p / ro))//
				+ nu_H3 * ((kv(u) + kv(v)) / 2.0 + (U_H3 / U_M_H3) * (p / ro))//
				+ nu_H4 * ((kv(u) + kv(v)) / 2.0 + (U_H4 / U_M_H4) * (p / ro));

			q1_H1 = (n_H_LISM_ / Kn_) * (k1 * S1 - nu_H1);
			q1_H2 = (n_H_LISM_ / Kn_) * (k2 * S1 - nu_H2);
			q1_H3 = (n_H_LISM_ / Kn_) * (k3 * S1 - nu_H3);
			q1_H4 = (n_H_LISM_ / Kn_) * (k4 * S1 - nu_H4);

			q21_H1 = (n_H_LISM_ / Kn_) * (k1 * S1 * u - nu_H1 * u_H1);
			q21_H2 = (n_H_LISM_ / Kn_) * (k2 * S1 * u - nu_H2 * u_H2);
			q21_H3 = (n_H_LISM_ / Kn_) * (k3 * S1 * u - nu_H3 * u_H3);
			q21_H4 = (n_H_LISM_ / Kn_) * (k4 * S1 * u - nu_H4 * u_H4);

			q22_H1 = (n_H_LISM_ / Kn_) * (k1 * S1 * v - nu_H1 * v_H1);
			q22_H2 = (n_H_LISM_ / Kn_) * (k2 * S1 * v - nu_H2 * v_H2);
			q22_H3 = (n_H_LISM_ / Kn_) * (k3 * S1 * v - nu_H3 * v_H3);
			q22_H4 = (n_H_LISM_ / Kn_) * (k4 * S1 * v - nu_H4 * v_H4);

			q3_H1 = (n_H_LISM_ / Kn_) * (k1 * S2 - nu_H1 * ((kv(u_H1) + kv(v_H1)) / 2.0 + (U_H1 / U_M_H1) * 2.0 * (p_H1 / ro_H1)));
			q3_H2 = (n_H_LISM_ / Kn_) * (k2 * S2 - nu_H2 * ((kv(u_H2) + kv(v_H2)) / 2.0 + (U_H2 / U_M_H2) * 2.0 * (p_H2 / ro_H2)));
			q3_H3 = (n_H_LISM_ / Kn_) * (k3 * S2 - nu_H3 * ((kv(u_H3) + kv(v_H3)) / 2.0 + (U_H3 / U_M_H3) * 2.0 * (p_H3 / ro_H3)));
			q3_H4 = (n_H_LISM_ / Kn_) * (k4 * S2 - nu_H4 * ((kv(u_H4) + kv(v_H4)) / 2.0 + (U_H4 / U_M_H4) * 2.0 * (p_H4 / ro_H4)));

			if (K->type == C_centr)
			{
				q1_H1 = 0.0;
				q1_H2 = 0.0;
				q1_H3 = 0.0;
				q1_H4 = 0.0;

				q21_H1 = 0.0;
				q21_H2 = 0.0;
				q21_H3 = 0.0;
				q21_H4 = 0.0;

				q22_H1 = 0.0;
				q22_H2 = 0.0;
				q22_H3 = 0.0;
				q22_H4 = 0.0;

				q3_H1 = 0.0;
				q3_H2 = 0.0;
				q3_H3 = 0.0;
				q3_H4 = 0.0;
			}

			//cout << "Setka.cpp    " << x << " " << y << " " << q21_H1 << " " << q22_H1 << " " << q3_H1 << " " << S2 << " " << nu_H1 * ((kv(u_H1) + kv(v_H1)) / 2.0 + (U_H1 / U_M_H1) * 2.0 * (p_H1 / ro_H1)) << endl;
			//cout << "Setka.cpp    " << x << " " << y << " " << par1.ro_H3 << " " << par1.ro_H4 << endl;
			//cout << "Setka.cpp    " << nu_H1 * ((kv(u) + kv(v)) / 2.0 + (U_H1 / U_M_H1) * (p / ro)) << " " << //
			//	nu_H2 * ((kv(u) + kv(v)) / 2.0 + (U_H2 / U_M_H2) * (p / ro)) << " " << //
			//	nu_H3 * ((kv(u) + kv(v)) / 2.0 + (U_H3 / U_M_H3) * (p / ro)) << " " << //
			//	nu_H4 * ((kv(u) + kv(v)) / 2.0 + (U_H4 / U_M_H4) * (p / ro)) << endl;
			//cout << "Setka.cpp    " << nu_H1 << " " << nu_H2 << " " << nu_H3 << " " << nu_H4 << endl;
			//cout << "Setka.cpp    " << U_H1 / U_M_H1 << " " << U_H2 / U_M_H2 << " " << U_H3 / U_M_H3 << " " << U_H4 / U_M_H4 << endl;
			/*if (radius < 100)
			{
				 q3_H1 = q22_H1 = q21_H1 = q1_H1 = 0.0;
			}*/


			if (K->type != C_centr)
			{
				ro3 = ro_H1 - T[now1] * (K->Potok_H1[0] / Volume + ro_H1 * v_H1 / y - q1_H1);
				if (ro3 <= 0)
				{
					printf("Problemsssss  ro H1 < 0!  %lf,  %lf\n", x, y);
					ro3 = 0.0003;
				}
				u3 = (ro_H1 * u_H1 - T[now1] * (K->Potok_H1[1] / Volume + ro_H1 * v_H1 * u_H1 / y - q21_H1)) / ro3;
				v3 = (ro_H1 * v_H1 - T[now1] * (K->Potok_H1[2] / Volume + ro_H1 * v_H1 * v_H1 / y - q22_H1)) / ro3;
				p3 = (((p_H1 / (ggg - 1) + ro_H1 * (u_H1 * u_H1 + v_H1 * v_H1) * 0.5) - T[now1] * (K->Potok_H1[3] / Volume + //
					+v_H1 * (ggg * p_H1 / (ggg - 1) + ro_H1 * (u_H1 * u_H1 + v_H1 * v_H1) * 0.5) / y - q3_H1)) - //
					0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
				if (p3 <= 0)
				{
					p3 = 0.000001;
				}

				//if (num_cell == 0)
				//{
				//	cout << "Setka.cpp    " << ro3 << " " << u3 << " " << v3 << " " << q1_H1 << endl;
				//	//ro3 = 0.0003;
				//}

				K->par[now2].ro_H1 = ro3;
				K->par[now2].p_H1 = p3;
				K->par[now2].u_H1 = u3;
				K->par[now2].v_H1 = v3;
			}
			else
			{
				K->par[now2].ro_H1 = K->par[now1].ro_H1;
				K->par[now2].p_H1 = K->par[now1].p_H1;
				K->par[now2].u_H1 = K->par[now1].u_H1;
				K->par[now2].v_H1 = K->par[now1].v_H1;
			}


			ro3 = ro_H2 - T[now1] * (K->Potok_H2[0] / Volume + ro_H2 * v_H2 / y - q1_H2);
			if (ro3 <= 0)
			{
				printf("Problemsssss  ro H2 < 0!  %lf,  %lf\n", x, y);
				ro3 = 0.000001;
			}
			u3 = (ro_H2 * u_H2 - T[now1] * (K->Potok_H2[1] / Volume + ro_H2 * v_H2 * u_H2 / y - q21_H2)) / ro3;
			v3 = (ro_H2 * v_H2 - T[now1] * (K->Potok_H2[2] / Volume + ro_H2 * v_H2 * v_H2 / y - q22_H2)) / ro3;
			p3 = (((p_H2 / (ggg - 1) + ro_H2 * (u_H2 * u_H2 + v_H2 * v_H2) * 0.5) - T[now1] * (K->Potok_H2[3] / Volume + //
				+v_H2 * (ggg * p_H2 / (ggg - 1) + ro_H2 * (u_H2 * u_H2 + v_H2 * v_H2) * 0.5) / y - q3_H2)) - //
				0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
			if (p3 <= 0)
			{
				p3 = 0.000001;
			}

			K->par[now2].ro_H2 = ro3;
			K->par[now2].p_H2 = p3;
			K->par[now2].u_H2 = u3;
			K->par[now2].v_H2 = v3;


			ro3 = ro_H3 - T[now1] * (K->Potok_H3[0] / Volume + ro_H3 * v_H3 / y - q1_H3);
			if (ro3 <= 0.0)
			{
				printf("Problemsssss  ro H3 < 0!  %lf,  %lf\n", x, y);
				ro3 = 0.000001;
			}
			u3 = (ro_H3 * u_H3 - T[now1] * (K->Potok_H3[1] / Volume + ro_H3 * v_H3 * u_H3 / y - q21_H3)) / ro3;
			v3 = (ro_H3 * v_H3 - T[now1] * (K->Potok_H3[2] / Volume + ro_H3 * v_H3 * v_H3 / y - q22_H3)) / ro3;
			p3 = (((p_H3 / (ggg - 1) + ro_H3 * (u_H3 * u_H3 + v_H3 * v_H3) * 0.5) - T[now1] * (K->Potok_H3[3] / Volume + //
				+v_H3 * (ggg * p_H3 / (ggg - 1) + ro_H3 * (u_H3 * u_H3 + v_H3 * v_H3) * 0.5) / y - q3_H3)) - //
				0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
			if (p3 <= 0)
			{
				p3 = 0.000001;
			}

			K->par[now2].ro_H3 = ro3;
			K->par[now2].p_H3 = p3;
			K->par[now2].u_H3 = u3;
			K->par[now2].v_H3 = v3;


			ro3 = ro_H4 - T[now1] * (K->Potok_H4[0] / Volume + ro_H4 * v_H4 / y - q1_H4);
			if (ro3 <= 0)
			{
				printf("Problemsssss  ro H4 < 0!  %lf,  %lf\n", x, y);
				ro3 = 0.000001;
			}
			u3 = (ro_H4 * u_H4 - T[now1] * (K->Potok_H4[1] / Volume + ro_H4 * v_H4 * u_H4 / y - q21_H4)) / ro3;
			v3 = (ro_H4 * v_H4 - T[now1] * (K->Potok_H4[2] / Volume + ro_H4 * v_H4 * v_H4 / y - q22_H4)) / ro3;
			p3 = (((p_H4 / (ggg - 1) + ro_H4 * (u_H4 * u_H4 + v_H4 * v_H4) * 0.5) - T[now1] * (K->Potok_H4[3] / Volume + //
				+v_H4 * (ggg * p_H4 / (ggg - 1) + ro_H4 * (u_H4 * u_H4 + v_H4 * v_H4) * 0.5) / y - q3_H4)) - //
				0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
			if (p3 <= 0)
			{
				p3 = 0.0000001;
			}

			K->par[now2].ro_H4 = ro3;
			K->par[now2].p_H4 = p3;
			K->par[now2].u_H4 = u3;
			K->par[now2].v_H4 = v3;


		}
	}
}

void Setka::Go_stationary_5_komponent_inner_2(int step)
{
	// Здесь объёмы - это поверхности вращения (считаем как-бы на декартовой сетке трёхмерной)
	cout << "Setka.cpp    " << "START " << step << endl;
	int now1 = 1;
	int now2 = 0;
	double T[2];
	T[0] = T[1] = 0.00000001;
	mutex mut;
	double x_min = 0.0, y_min = 0.0;

	for (int st = 0; st < step; st++)
	{
		if (st % 10000 == 0 && st > 0)
		{
			cout << "Setka.cpp    " << st << " " << T[now2] << " " << x_min << " " << y_min << endl;
		}
		now1 = (now1 + 1) % 2; // Какие параметры сейчас берём
		now2 = (now2 + 1) % 2; // Какие параметры сейчас меняем
		double time = 10000000;
		T[now2] = 100000000;

		//omp_set_num_threads(4);
#pragma omp parallel for
		for (int num_cell = 0; num_cell < this->All_Cells_Inner.size(); num_cell++)
		{
			auto K = this->All_Cells_Inner[num_cell];
			if (K->type == C_centr)
			{
				continue;
			}
			K->Potok[0] = K->Potok[1] = K->Potok[2] = K->Potok[3] = K->Potok[4] = 0.0;
			K->Potok_H1[0] = K->Potok_H1[1] = K->Potok_H1[2] = K->Potok_H1[3] = 0.0;
			K->Potok_H2[0] = K->Potok_H2[1] = K->Potok_H2[2] = K->Potok_H2[3] = 0.0;
			K->Potok_H3[0] = K->Potok_H3[1] = K->Potok_H3[2] = K->Potok_H3[3] = 0.0;
			K->Potok_H4[0] = K->Potok_H4[1] = K->Potok_H4[2] = K->Potok_H4[3] = 0.0;

			double P[4];
			double PQ;
			Parametr par1 = K->par[now1];
			double n1, n2;
			double dist;
			double x, y;
			K->Get_Center(x, y);

			double vs1 = K->Get_Volume();
			double vs2 = K->Get_Volume_posle();

			//y = K->Get_Volume_rotate(360.0) / (2.0 * pi_ * vs1);
			//double y2 = K->Get_Volume_posle_rotate(360.0) / (2.0 * pi_ * vs2);

			double radius = sqrt(kv(x) + kv(y));
			double Volume = K->Get_Volume_rotate(360.0) / (2.0 * pi_);
			double Volume2 = K->Get_Volume_posle_rotate(360.0) / (2.0 * pi_);

			for (auto& i : K->Grans)
			{
				double x2, y2;
				double S = i->Get_square() * (i->A->y + i->B->y) * 0.5;
				i->Get_Center(x2, y2);
				dist = sqrt(kv(x - x2) + kv(y - y2));
				Parametr par2, par11, par3;
				par11 = par1;
				i->Get_par(par2, now1);
				i->Get_normal(n1, n2);

				double Vc, Vl, Vp;
				int met = 1;
				double dis = sqrt(kv(x2) + kv(y2));

				if (i->type == Usualy)
				{
					double x3, y3;
					i->Sosed->Get_Center(x3, y3);
					double dd = sqrt(kv(x3) + kv(y3));
					par11.ro = par11.ro * kv(radius) / kv(dis);
					par2.ro = par2.ro * kv(dd) / kv(dis);
					par11.ro_H1 = par11.ro_H1 * pow(radius / dis, H_pow);
					par2.ro_H1 = par2.ro_H1 * pow(dd / dis, H_pow);
					//par11.p_H1 = par11.p_H1 * pow(radius / dis, 2 * ggg);
					//par2.p_H1 = par2.p_H1 * pow(dd / dis, 2 * ggg);
					par11.Q = par11.Q * kv(radius) / kv(dis);
					par2.Q = par2.Q * kv(dd) / kv(dis);
					par11.p = par11.p * pow(radius / dis, 2 * ggg);
					par2.p = par2.p * pow(dd / dis, 2 * ggg);
					polar_perenos(x, y, x2, y2, par11.u, par11.v);
					polar_perenos(x3, y3, x2, y2, par2.u, par2.v);
					polar_perenos(x, y, x2, y2, par11.u_H1, par11.v_H1);
					polar_perenos(x3, y3, x2, y2, par2.u_H1, par2.v_H1);

					if (K->type == C_1 && i->Sosed->type == C_1)
					{
						Parametr pp1, pp2;
						i->Get_par_TVD_radial(pp1, now1);
						i->Gran_copy->Get_par_TVD_radial(pp2, now1);
						//par11.u_H1 = pp1.u_H1;
						//par11.v_H1 = pp1.v_H1;
						//par2.u_H1 = pp2.u_H1;
						//par2.v_H1 = pp2.v_H1;
						par11.p_H1 = pp1.p_H1;
						par2.p_H1 = pp2.p_H1;
						par11.ro_H1 = pp1.ro_H1;
						par2.ro_H1 = pp2.ro_H1;

						par11.ro = pp1.ro;
						par2.ro = pp2.ro;
						par11.p = pp1.p;
						par2.p = pp2.p;
					}
				}
				else if (i->type == Inner_sphere)
				{
					/*par11.ro = par11.ro * kv(radius) / kv(dis);
					par11.ro_H1 = par11.ro_H1 * pow(radius / dis, H_pow);
					par11.Q = par11.Q * kv(radius) / kv(dis);
					par11.p = par11.p * pow(radius / dis, 2 * ggg);*/

					//cout << "Setka.cpp    " << "inner " << dis << endl;
					/*par2.ro_H3 = par1.ro_H3;
					par2.p_H3 = par1.p_H3;
					par2.u_H3 = par1.u_H3;
					par2.v_H3 = par1.v_H3;
					par2.ro_H4 = par1.ro_H4;
					par2.p_H4 = par1.p_H4;
					par2.u_H4 = par1.u_H4;
					par2.v_H4 = par1.v_H4;*/
					par11.ro = par2.ro;
					par11.Q = par2.Q;
					par11.p = par2.p;
					par11.u = par2.u;
					par11.v = par2.v;
					par11.ro_H1 = par2.ro_H1;
					par11.p_H1 = par2.p_H1;
					par11.u_H1 = par2.u_H1;
					par11.v_H1 = par2.v_H1;
				}
				else if (i->type == Axis)
				{
					par11.ro = par11.ro * kv(radius) / kv(dis);
					par11.Q = par11.Q * kv(radius) / kv(dis);
					par11.p = par11.p * pow(radius / dis, 2 * ggg);
					par11.ro_H1 = par11.ro_H1 * pow(radius / dis, H_pow);
					par11.p_H1 = par11.p_H1 * pow(radius / dis, 2 * ggg);
					polar_perenos(x, y, x2, y2, par11.u, par11.v);
					polar_perenos(x, y, x2, y2, par11.u_H1, par11.v_H1);
					par2 = par11;
					par2.v = -par2.v;
					par2.v_H1 = -par2.v_H1;
					par2.v_H2 = -par2.v_H2;
					par2.v_H3 = -par2.v_H3;
					par2.v_H4 = -par2.v_H4;
				}

				if (K->type != C_centr)
				{
					time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro_H1, par11.Q, par11.p_H1, par11.u_H1, par11.v_H1, par2.ro_H1, par2.Q, //
						par2.p_H1, par2.u_H1, par2.v_H1, 0.0, P, PQ, n1, n2, dist, met, Vl, Vc, Vp));
					for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
					{
						K->Potok_H1[k] = K->Potok_H1[k] + P[k] * S;
					}
				}

				time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro_H2, par11.Q, par11.p_H2, par11.u_H2, par11.v_H2, par2.ro_H2, par2.Q, //
					par2.p_H2, par2.u_H2, par2.v_H2, 0.0, P, PQ, n1, n2, dist, met, Vl, Vc, Vp));
				for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
				{
					K->Potok_H2[k] = K->Potok_H2[k] + P[k] * S;
				}

				time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro_H3, par11.Q, par11.p_H3, par11.u_H3, par11.v_H3, par2.ro_H3, par2.Q, //
					par2.p_H3, par2.u_H3, par2.v_H3, 0.0, P, PQ, n1, n2, dist, met, Vl, Vc, Vp));
				for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
				{
					K->Potok_H3[k] = K->Potok_H3[k] + P[k] * S;
				}

				time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro_H4, par11.Q, par11.p_H4, par11.u_H4, par11.v_H4, par2.ro_H4, par2.Q, //
					par2.p_H4, par2.u_H4, par2.v_H4, 0.0, P, PQ, n1, n2, dist, met, Vl, Vc, Vp));
				for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
				{
					K->Potok_H4[k] = K->Potok_H4[k] + P[k] * S;
				}


				if (K->type != C_centr)
				{
					/*double Vl, Vp, VV;
					vector<double> qqq1(5);
					vector<double> qqq2(5);
					vector<double> qqq(5);
					vector<double> n(3);
					n[0] = n1;
					n[1] = n2;
					n[2] = 0.0;

					qqq1[0] = par11.ro;
					qqq1[1] = par11.p;
					qqq1[2] = par11.u;
					qqq1[3] = par11.v;
					qqq1[4] = 0.0;

					qqq2[0] = par2.ro;
					qqq2[1] = par2.p;
					qqq2[2] = par2.u;
					qqq2[3] = par2.v;
					qqq2[4] = 0.0;

					auto SS = Solvers();
					time = min(time, SS.Godunov_Solver_Alexashov(qqq1, qqq2, n, qqq, Vl, Vp, VV, 0.0, dist));

					K->Potok[0] = K->Potok[0] + qqq[0] * S;
					K->Potok[1] = K->Potok[1] + qqq[1] * S;
					K->Potok[2] = K->Potok[2] + qqq[2] * S;
					K->Potok[3] = K->Potok[3] + qqq[4] * S;*/


					time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro, par11.Q, par11.p, par11.u, par11.v, par2.ro, par2.Q, //
						par2.p, par2.u, par2.v, 0.0, P, PQ, n1, n2, dist, met, Vl, Vc, Vp));
					for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
					{
						K->Potok[k] = K->Potok[k] + P[k] * S;
					}
					K->Potok[4] = K->Potok[4] + PQ * S;
				}
			}

			if (time < T[now2])
			{
				mut.lock();
				T[now2] = time;
				x_min = x;
				y_min = y;
				mut.unlock();
			}



			double U_M_H1, U_M_H2, U_M_H3, U_M_H4;
			double U_H1, U_H2, U_H3, U_H4;
			double sigma_H1, sigma_H2, sigma_H3, sigma_H4;
			double nu_H1, nu_H2, nu_H3, nu_H4;
			double q2_1, q2_2, q3;
			double u, v, ro, p, Q;

			u = par1.u;// *(chi_real / chi_);
			v = par1.v;// *(chi_real / chi_);
			ro = par1.ro;// / kv(chi_real / chi_);
			p = par1.p;
			Q = par1.Q;// / kv(chi_real / chi_);

			double u_H1 = par1.u_H1, v_H1 = par1.v_H1, ro_H1 = par1.ro_H1, p_H1 = par1.p_H1;
			double u_H2 = par1.u_H2, v_H2 = par1.v_H2, ro_H2 = par1.ro_H2, p_H2 = par1.p_H2;
			double u_H3 = par1.u_H3, v_H3 = par1.v_H3, ro_H3 = par1.ro_H3, p_H3 = par1.p_H3;
			double u_H4 = par1.u_H4, v_H4 = par1.v_H4, ro_H4 = par1.ro_H4, p_H4 = par1.p_H4;

			U_M_H1 = sqrt(kv(u - u_H1) + kv(v - v_H1) + (64.0 / (9.0 * pi_)) //
				* (p / ro + 2.0 * p_H1 / ro_H1));
			U_M_H2 = sqrt(kv(u - u_H2) + kv(v - v_H2) + (64.0 / (9.0 * pi_)) //
				* (p / ro + 2.0 * p_H2 / ro_H2));
			U_M_H3 = sqrt(kv(u - u_H3) + kv(v - v_H3) + (64.0 / (9.0 * pi_)) //
				* (p / ro + 2.0 * p_H3 / ro_H3));
			U_M_H4 = sqrt(kv(u - u_H4) + kv(v - v_H4) + (64.0 / (9.0 * pi_)) //
				* (p / ro + 2.0 * p_H4 / ro_H4));

			U_H1 = sqrt(kv(u - u_H1) + kv(v - v_H1) + (4.0 / pi_) //
				* (p / ro + 2.0 * p_H1 / ro_H1));
			U_H2 = sqrt(kv(u - u_H2) + kv(v - v_H2) + (4.0 / pi_) //
				* (p / ro + 2.0 * p_H2 / ro_H2));
			U_H3 = sqrt(kv(u - u_H3) + kv(v - v_H3) + (4.0 / pi_) //
				* (p / ro + 2.0 * p_H3 / ro_H3));
			U_H4 = sqrt(kv(u - u_H4) + kv(v - v_H4) + (4.0 / pi_) //
				* (p / ro + 2.0 * p_H4 / ro_H4));

			sigma_H1 = kv(1.0 - a_2 * log(U_M_H1)); // 0.1243
			sigma_H2 = kv(1.0 - a_2 * log(U_M_H2));
			sigma_H3 = kv(1.0 - a_2 * log(U_M_H3)); // 0.1121     a_2
			sigma_H4 = kv(1.0 - a_2 * log(U_M_H4));

			nu_H1 = ro * ro_H1 * U_M_H1 * sigma_H1;
			nu_H2 = ro * ro_H2 * U_M_H2 * sigma_H2;
			nu_H3 = ro * ro_H3 * U_M_H3 * sigma_H3;
			nu_H4 = ro * ro_H4 * U_M_H4 * sigma_H4;

			q2_1 = (n_p_LISM_ / Kn_) * (nu_H1 * (u_H1 - u) + nu_H2 * (u_H2 - u) //
				+ nu_H3 * (u_H3 - u) + nu_H4 * (u_H4 - u));
			q2_2 = (n_p_LISM_ / Kn_) * (nu_H1 * (v_H1 - v) + nu_H2 * (v_H2 - v) //
				+ nu_H3 * (v_H3 - v) + nu_H4 * (v_H4 - v));
			q3 = (n_p_LISM_ / Kn_) * (nu_H1 * ((kv(u_H1) + kv(v_H1) - kv(u) - kv(v)) / 2.0 + //
				(U_H1 / U_M_H1) * (2.0 * p_H1 / ro_H1 - p / ro)) + //
				nu_H2 * ((kv(u_H2) + kv(v_H2) - kv(u) - kv(v)) / 2.0 + //
					(U_H2 / U_M_H2) * (2.0 * p_H2 / ro_H2 - p / ro)) + //
				nu_H3 * ((kv(u_H3) + kv(v_H3) - kv(u) - kv(v)) / 2.0 + //
					(U_H3 / U_M_H3) * (2.0 * p_H3 / ro_H3 - p / ro)) + //
				nu_H4 * ((kv(u_H4) + kv(v_H4) - kv(u) - kv(v)) / 2.0 + //
					(U_H4 / U_M_H4) * (2.0 * p_H4 / ro_H4 - p / ro)));// / (chi_real / chi_);


			/*q2_1 = 0.0;
			q2_2 = 0.0;
			q3 = 0.0;*/

			double ro2_p, p2_p, V1_p, V2_p, QQ2;
			double ro3, p3, u3, v3, Q33;
			double kappa = Q / ro;

			u = par1.u;
			v = par1.v;
			ro = par1.ro;
			p = par1.p;
			Q = par1.Q;

			if (K->type != C_centr)
			{
				/*ro3 = ro - T[now1] * (K->Potok[0] / Volume + ro * v / y);
				Q33 = Q - (T[now1] / Volume) * K->Potok[4] - T[now1] * Q * v / y;*/
				ro3 = ro - T[now1] * (K->Potok[0] / Volume);
				Q33 = Q - (T[now1] / Volume) * K->Potok[4] - T[now1];
				if (ro3 <= 0)
				{
					printf("Problemsssss  ro < 0!   %lf,  %lf\n", x, y);
					ro3 = 0.00001;
				}
				//u3 = (ro * u - T[now1] * (K->Potok[1] / Volume + ro * v * u / y - q2_1)) / ro3;
				//v3 = (ro * v - T[now1] * (K->Potok[2] / Volume + ro * v * v / y - q2_2)) / ro3;
				//p3 = (((p / (ggg - 1) + ro * (u * u + v * v) * 0.5) - T[now1] * (K->Potok[3] / Volume + //
				//	+v * (ggg * p / (ggg - 1) + ro * (u * u + v * v) * 0.5) / y - q3)) - //
				//	0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
				u3 = (ro * u - T[now1] * (K->Potok[1] / Volume - q2_1)) / ro3;
				v3 = (ro * v - T[now1] * ((K->Potok[2] - p * 0.5 * (vs1 + vs2)) / Volume - q2_2)) / ro3;
				p3 = (((p / (ggg - 1) + ro * (u * u + v * v) * 0.5) - T[now1] * (K->Potok[3] / Volume  - q3)) - //
					0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
				if (p3 <= 0)
				{
					p3 = 0.000001;
				}

				K->par[now2].ro = ro3;
				K->par[now2].Q = Q33;
				K->par[now2].p = p3;
				K->par[now2].u = u3;
				K->par[now2].v = v3;
			}
			else
			{
				K->par[now2].ro = K->par[now1].ro;
				K->par[now2].Q = K->par[now1].Q;
				K->par[now2].p = K->par[now1].p;
				K->par[now2].u = K->par[now1].u;
				K->par[now2].v = K->par[now1].v;
			}

			u = par1.u;// *(chi_real / chi_);
			v = par1.v;// *(chi_real / chi_);
			ro = par1.ro;// / kv(chi_real / chi_);
			p = par1.p;
			Q = par1.Q;// / kv(chi_real / chi_);


			int k1 = 1, k2 = 0, k3 = 0, k4 = 0;


			double S1, S2;
			double q1_H1, q1_H2, q1_H3, q1_H4;
			double q21_H1, q21_H2, q21_H3, q21_H4;
			double q22_H1, q22_H2, q22_H3, q22_H4;
			double q3_H1, q3_H2, q3_H3, q3_H4;

			S1 = nu_H1 + nu_H2 + nu_H3 + nu_H4;
			S2 = nu_H1 * ((kv(u) + kv(v)) / 2.0 + (U_H1 / U_M_H1) * (p / ro))//
				+ nu_H2 * ((kv(u) + kv(v)) / 2.0 + (U_H2 / U_M_H2) * (p / ro))//
				+ nu_H3 * ((kv(u) + kv(v)) / 2.0 + (U_H3 / U_M_H3) * (p / ro))//
				+ nu_H4 * ((kv(u) + kv(v)) / 2.0 + (U_H4 / U_M_H4) * (p / ro));

			q1_H1 = (n_H_LISM_ / Kn_) * (k1 * S1 - nu_H1);
			q1_H2 = (n_H_LISM_ / Kn_) * (k2 * S1 - nu_H2);
			q1_H3 = (n_H_LISM_ / Kn_) * (k3 * S1 - nu_H3);
			q1_H4 = (n_H_LISM_ / Kn_) * (k4 * S1 - nu_H4);

			q21_H1 = (n_H_LISM_ / Kn_) * (k1 * S1 * u - nu_H1 * u_H1);
			q21_H2 = (n_H_LISM_ / Kn_) * (k2 * S1 * u - nu_H2 * u_H2);
			q21_H3 = (n_H_LISM_ / Kn_) * (k3 * S1 * u - nu_H3 * u_H3);
			q21_H4 = (n_H_LISM_ / Kn_) * (k4 * S1 * u - nu_H4 * u_H4);

			q22_H1 = (n_H_LISM_ / Kn_) * (k1 * S1 * v - nu_H1 * v_H1);
			q22_H2 = (n_H_LISM_ / Kn_) * (k2 * S1 * v - nu_H2 * v_H2);
			q22_H3 = (n_H_LISM_ / Kn_) * (k3 * S1 * v - nu_H3 * v_H3);
			q22_H4 = (n_H_LISM_ / Kn_) * (k4 * S1 * v - nu_H4 * v_H4);

			q3_H1 = (n_H_LISM_ / Kn_) * (k1 * S2 - nu_H1 * ((kv(u_H1) + kv(v_H1)) / 2.0 + (U_H1 / U_M_H1) * 2.0 * (p_H1 / ro_H1)));
			q3_H2 = (n_H_LISM_ / Kn_) * (k2 * S2 - nu_H2 * ((kv(u_H2) + kv(v_H2)) / 2.0 + (U_H2 / U_M_H2) * 2.0 * (p_H2 / ro_H2)));
			q3_H3 = (n_H_LISM_ / Kn_) * (k3 * S2 - nu_H3 * ((kv(u_H3) + kv(v_H3)) / 2.0 + (U_H3 / U_M_H3) * 2.0 * (p_H3 / ro_H3)));
			q3_H4 = (n_H_LISM_ / Kn_) * (k4 * S2 - nu_H4 * ((kv(u_H4) + kv(v_H4)) / 2.0 + (U_H4 / U_M_H4) * 2.0 * (p_H4 / ro_H4)));

			if (K->type == C_centr)
			{
				q1_H1 = 0.0;
				q1_H2 = 0.0;
				q1_H3 = 0.0;
				q1_H4 = 0.0;

				q21_H1 = 0.0;
				q21_H2 = 0.0;
				q21_H3 = 0.0;
				q21_H4 = 0.0;

				q22_H1 = 0.0;
				q22_H2 = 0.0;
				q22_H3 = 0.0;
				q22_H4 = 0.0;

				q3_H1 = 0.0;
				q3_H2 = 0.0;
				q3_H3 = 0.0;
				q3_H4 = 0.0;
			}

			//cout << "Setka.cpp    " << x << " " << y << " " << q21_H1 << " " << q22_H1 << " " << q3_H1 << " " << S2 << " " << nu_H1 * ((kv(u_H1) + kv(v_H1)) / 2.0 + (U_H1 / U_M_H1) * 2.0 * (p_H1 / ro_H1)) << endl;
			//cout << "Setka.cpp    " << x << " " << y << " " << par1.ro_H3 << " " << par1.ro_H4 << endl;
			//cout << "Setka.cpp    " << nu_H1 * ((kv(u) + kv(v)) / 2.0 + (U_H1 / U_M_H1) * (p / ro)) << " " << //
			//	nu_H2 * ((kv(u) + kv(v)) / 2.0 + (U_H2 / U_M_H2) * (p / ro)) << " " << //
			//	nu_H3 * ((kv(u) + kv(v)) / 2.0 + (U_H3 / U_M_H3) * (p / ro)) << " " << //
			//	nu_H4 * ((kv(u) + kv(v)) / 2.0 + (U_H4 / U_M_H4) * (p / ro)) << endl;
			//cout << "Setka.cpp    " << nu_H1 << " " << nu_H2 << " " << nu_H3 << " " << nu_H4 << endl;
			//cout << "Setka.cpp    " << U_H1 / U_M_H1 << " " << U_H2 / U_M_H2 << " " << U_H3 / U_M_H3 << " " << U_H4 / U_M_H4 << endl;
			/*if (radius < 100)
			{
				 q3_H1 = q22_H1 = q21_H1 = q1_H1 = 0.0;
			}*/


			if (K->type != C_centr)
			{
				//ro3 = ro_H1 - T[now1] * (K->Potok_H1[0] / Volume + ro_H1 * v_H1 / y - q1_H1);
				ro3 = ro_H1 - T[now1] * (K->Potok_H1[0] / Volume - q1_H1);
				if (ro3 <= 0)
				{
					printf("Problemsssss  ro H1 < 0!  %lf,  %lf\n", x, y);
					ro3 = 0.0003;
				}
				//u3 = (ro_H1 * u_H1 - T[now1] * (K->Potok_H1[1] / Volume + ro_H1 * v_H1 * u_H1 / y - q21_H1)) / ro3;
				//v3 = (ro_H1 * v_H1 - T[now1] * (K->Potok_H1[2] / Volume + ro_H1 * v_H1 * v_H1 / y - q22_H1)) / ro3;
				//p3 = (((p_H1 / (ggg - 1) + ro_H1 * (u_H1 * u_H1 + v_H1 * v_H1) * 0.5) - T[now1] * (K->Potok_H1[3] / Volume + //
				//	+v_H1 * (ggg * p_H1 / (ggg - 1) + ro_H1 * (u_H1 * u_H1 + v_H1 * v_H1) * 0.5) / y - q3_H1)) - //
				//	0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
				u3 = (ro_H1 * u_H1 - T[now1] * (K->Potok_H1[1] / Volume - q21_H1)) / ro3;
				v3 = (ro_H1 * v_H1 - T[now1] * ((K->Potok_H1[2] - p_H1 * 0.5 * (vs1 + vs2)) / Volume - q22_H1)) / ro3;
				p3 = (((p_H1 / (ggg - 1) + ro_H1 * (u_H1 * u_H1 + v_H1 * v_H1) * 0.5) - T[now1] * (K->Potok_H1[3] / Volume - q3_H1)) - //
					0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);

				if (p3 <= 0)
				{
					p3 = 0.000001;
				}

				//if (num_cell == 0)
				//{
				//	cout << "Setka.cpp    " << ro3 << " " << u3 << " " << v3 << " " << q1_H1 << endl;
				//	//ro3 = 0.0003;
				//}

				K->par[now2].ro_H1 = ro3;
				K->par[now2].p_H1 = p3;
				K->par[now2].u_H1 = u3;
				K->par[now2].v_H1 = v3;
			}
			else
			{
				K->par[now2].ro_H1 = K->par[now1].ro_H1;
				K->par[now2].p_H1 = K->par[now1].p_H1;
				K->par[now2].u_H1 = K->par[now1].u_H1;
				K->par[now2].v_H1 = K->par[now1].v_H1;
			}


			//ro3 = ro_H2 - T[now1] * (K->Potok_H2[0] / Volume + ro_H2 * v_H2 / y - q1_H2);
			ro3 = ro_H2 - T[now1] * (K->Potok_H2[0] / Volume - q1_H2);
			if (ro3 <= 0)
			{
				printf("Problemsssss  ro H2 < 0!  %lf,  %lf\n", x, y);
				ro3 = 0.000001;
			}
			//u3 = (ro_H2 * u_H2 - T[now1] * (K->Potok_H2[1] / Volume + ro_H2 * v_H2 * u_H2 / y - q21_H2)) / ro3;
			//v3 = (ro_H2 * v_H2 - T[now1] * (K->Potok_H2[2] / Volume + ro_H2 * v_H2 * v_H2 / y - q22_H2)) / ro3;
			//p3 = (((p_H2 / (ggg - 1) + ro_H2 * (u_H2 * u_H2 + v_H2 * v_H2) * 0.5) - T[now1] * (K->Potok_H2[3] / Volume + //
			//	+v_H2 * (ggg * p_H2 / (ggg - 1) + ro_H2 * (u_H2 * u_H2 + v_H2 * v_H2) * 0.5) / y - q3_H2)) - //
			//	0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
			u3 = (ro_H2 * u_H2 - T[now1] * (K->Potok_H2[1] / Volume - q21_H2)) / ro3;
			v3 = (ro_H2 * v_H2 - T[now1] * ((K->Potok_H2[2] - p_H2 * 0.5 * (vs1 + vs2)) / Volume - q22_H2)) / ro3;
			p3 = (((p_H2 / (ggg - 1) + ro_H2 * (u_H2 * u_H2 + v_H2 * v_H2) * 0.5) - T[now1] * (K->Potok_H2[3] / Volume - q3_H2)) - //
				0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);

			if (p3 <= 0)
			{
				p3 = 0.000001;
			}

			K->par[now2].ro_H2 = ro3;
			K->par[now2].p_H2 = p3;
			K->par[now2].u_H2 = u3;
			K->par[now2].v_H2 = v3;


			//ro3 = ro_H3 - T[now1] * (K->Potok_H3[0] / Volume + ro_H3 * v_H3 / y - q1_H3);
			ro3 = ro_H3 - T[now1] * (K->Potok_H3[0] / Volume - q1_H3);
			if (ro3 <= 0.0)
			{
				printf("Problemsssss  ro H3 < 0!  %lf,  %lf\n", x, y);
				ro3 = 0.000001;
			}
			//u3 = (ro_H3 * u_H3 - T[now1] * (K->Potok_H3[1] / Volume + ro_H3 * v_H3 * u_H3 / y - q21_H3)) / ro3;
			//v3 = (ro_H3 * v_H3 - T[now1] * (K->Potok_H3[2] / Volume + ro_H3 * v_H3 * v_H3 / y - q22_H3)) / ro3;
			//p3 = (((p_H3 / (ggg - 1) + ro_H3 * (u_H3 * u_H3 + v_H3 * v_H3) * 0.5) - T[now1] * (K->Potok_H3[3] / Volume + //
			//	+v_H3 * (ggg * p_H3 / (ggg - 1) + ro_H3 * (u_H3 * u_H3 + v_H3 * v_H3) * 0.5) / y - q3_H3)) - //
			//	0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
			u3 = (ro_H3 * u_H3 - T[now1] * (K->Potok_H3[1] / Volume - q21_H3)) / ro3;
			v3 = (ro_H3 * v_H3 - T[now1] * ((K->Potok_H3[2] - p_H3 * 0.5 * (vs1 + vs2)) / Volume - q22_H3)) / ro3;
			p3 = (((p_H3 / (ggg - 1) + ro_H3 * (u_H3 * u_H3 + v_H3 * v_H3) * 0.5) - T[now1] * (K->Potok_H3[3] / Volume - q3_H3)) - //
				0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
			if (p3 <= 0)
			{
				p3 = 0.000001;
			}

			K->par[now2].ro_H3 = ro3;
			K->par[now2].p_H3 = p3;
			K->par[now2].u_H3 = u3;
			K->par[now2].v_H3 = v3;


			//ro3 = ro_H4 - T[now1] * (K->Potok_H4[0] / Volume + ro_H4 * v_H4 / y - q1_H4);
			ro3 = ro_H4 - T[now1] * (K->Potok_H4[0] / Volume - q1_H4);
			if (ro3 <= 0)
			{
				printf("Problemsssss  ro H4 < 0!  %lf,  %lf\n", x, y);
				ro3 = 0.000001;
			}
			//u3 = (ro_H4 * u_H4 - T[now1] * (K->Potok_H4[1] / Volume + ro_H4 * v_H4 * u_H4 / y - q21_H4)) / ro3;
			//v3 = (ro_H4 * v_H4 - T[now1] * (K->Potok_H4[2] / Volume + ro_H4 * v_H4 * v_H4 / y - q22_H4)) / ro3;
			//p3 = (((p_H4 / (ggg - 1) + ro_H4 * (u_H4 * u_H4 + v_H4 * v_H4) * 0.5) - T[now1] * (K->Potok_H4[3] / Volume + //
			//	+v_H4 * (ggg * p_H4 / (ggg - 1) + ro_H4 * (u_H4 * u_H4 + v_H4 * v_H4) * 0.5) / y - q3_H4)) - //
			//	0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
			u3 = (ro_H4 * u_H4 - T[now1] * (K->Potok_H4[1] / Volume - q21_H4)) / ro3;
			v3 = (ro_H4 * v_H4 - T[now1] * ((K->Potok_H4[2] - p_H4 * 0.5 * (vs1 + vs2)) / Volume - q22_H4)) / ro3;
			p3 = (((p_H4 / (ggg - 1) + ro_H4 * (u_H4 * u_H4 + v_H4 * v_H4) * 0.5) - T[now1] * (K->Potok_H4[3] / Volume - q3_H4)) - //
				0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
			if (p3 <= 0)
			{
				p3 = 0.0000001;
			}

			K->par[now2].ro_H4 = ro3;
			K->par[now2].p_H4 = p3;
			K->par[now2].u_H4 = u3;
			K->par[now2].v_H4 = v3;


		}
	}
}

void Setka::Go_stationary_5_komponent_inner_MK2(int step)
{
	cout << "Setka.cpp    " << "START " << step << endl;
	int now1 = 1;
	int now2 = 0;
	double T[2];
	T[0] = T[1] = 0.00000001;
	mutex mut;
	double x_min = 0.0, y_min = 0.0;

	for (int st = 0; st < step; st++)
	{
		if (st % 1000 == 0 && st > 0)
		{
			cout << "Setka.cpp    " << st << " " << T[now2] << " " << x_min << " " << y_min << endl;
		}
		now1 = (now1 + 1) % 2; // Какие параметры сейчас берём
		now2 = (now2 + 1) % 2; // Какие параметры сейчас меняем
		double time = 10000000;
		T[now2] = 100000000;

		//omp_set_num_threads(4);
#pragma omp parallel for
		for (int num_cell = 0; num_cell < this->All_Cells_Inner.size(); num_cell++)
		{
			auto K = this->All_Cells_Inner[num_cell];
			if (K->type == C_centr)
			{
				continue;
			}
			K->Potok[0] = K->Potok[1] = K->Potok[2] = K->Potok[3] = K->Potok[4] = 0.0;

			double P[4];
			double PQ;
			Parametr par1 = K->par[now1];
			double n1, n2;
			double dist;
			double x, y;
			K->Get_Center(x, y);

			double vs1 = K->Get_Volume();
			double vs2 = K->Get_Volume_posle();

			//y = K->Get_Volume_rotate(360.0) / (2.0 * pi_ * vs1);
			//double y2 = K->Get_Volume_posle_rotate(360.0) / (2.0 * pi_ * vs2);

			double radius = sqrt(kv(x) + kv(y));
			double Volume = K->Get_Volume_rotate(360.0) / (2.0 * pi_);
			double Volume2 = K->Get_Volume_posle_rotate(360.0) / (2.0 * pi_);

			for (auto& i : K->Grans)
			{
				double x2, y2;
				double S = i->Get_square() * (i->A->y + i->B->y) * 0.5;
				i->Get_Center(x2, y2);
				dist = sqrt(kv(x - x2) + kv(y - y2));
				Parametr par2, par11, par3;
				par11 = par1;
				i->Get_par(par2, now1);
				i->Get_normal(n1, n2);

				double Vc, Vl, Vp;
				int met = 1;
				double dis = sqrt(kv(x2) + kv(y2));

				if (i->type == Usualy)
				{
					double x3, y3;
					i->Sosed->Get_Center(x3, y3);
					double dd = sqrt(kv(x3) + kv(y3));
					par11.ro = par11.ro * kv(radius) / kv(dis);
					par2.ro = par2.ro * kv(dd) / kv(dis);
					//par11.p_H1 = par11.p_H1 * pow(radius / dis, 2 * ggg);
					//par2.p_H1 = par2.p_H1 * pow(dd / dis, 2 * ggg);
					par11.Q = par11.Q * kv(radius) / kv(dis);
					par2.Q = par2.Q * kv(dd) / kv(dis);
					par11.p = par11.p * pow(radius / dis, 2.0 * ggg);
					par2.p = par2.p * pow(dd / dis, 2.0 * ggg);
					polar_perenos(x, y, x2, y2, par11.u, par11.v);
					polar_perenos(x3, y3, x2, y2, par2.u, par2.v);

					if (K->type == C_1 && i->Sosed->type == C_1)
					{
						Parametr pp1, pp2;
						i->Get_par_TVD_radial(pp1, now1);
						i->Gran_copy->Get_par_TVD_radial(pp2, now1);

						par11.ro = pp1.ro;
						par2.ro = pp2.ro;
						par11.p = pp1.p;
						par2.p = pp2.p;
					}
				}
				else if (i->type == Inner_sphere)
				{
					par11.ro = par2.ro;
					par11.Q = par2.Q;
					par11.p = par2.p;
					par11.u = par2.u;
					par11.v = par2.v;
				}
				else if (i->type == Axis)
				{
					// Сюда не попадаем, т.к. вращательная симметрия

					/*par11.ro = par11.ro * kv(radius) / kv(dis);
					par11.Q = par11.Q * kv(radius) / kv(dis);
					par11.p = par11.p * pow(radius / dis, 2.0 * ggg);
					polar_perenos(x, y, x2, y2, par11.u, par11.v);
					par2 = par11;
					par2.v = -par2.v;
					met = 0;*/

					i->Get_par_TVD(par2, now1);
					par11 = par2;
					polar_perenos(x, y, x2, y2, par11.u, par11.v);
					par2 = par11;
					par2.v = -par2.v;

				}

				if (K->type != C_centr)
				{
					time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro, par11.Q, par11.p, par11.u, par11.v, par2.ro, par2.Q, //
						par2.p, par2.u, par2.v, 0.0, P, PQ, n1, n2, dist, met, Vl, Vc, Vp));
					for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
					{
						K->Potok[k] = K->Potok[k] + P[k] * S;
					}
					K->Potok[4] = K->Potok[4] + PQ * S;
				}
			}

			if (time < T[now2])
			{
				mut.lock();
				T[now2] = time;
				x_min = x;
				y_min = y;
				mut.unlock();
			}



			double u = par1.u;
			double v = par1.v;
			double ro = par1.ro;
			double p = par1.p;
			double Q = par1.Q;

			double a1, a2, a3;
			double b1, b2, b3;
			double c1, c2, c3;
			double q2_1, q2_2, q3;

			if (polusum == true)
			{
				K->Get_Sourse_MK1(a1, a2, a3, u, v, ro, p);
				K->Get_Sourse_MK2(b1, b2, b3, u, v, ro, p);
				K->Get_Sourse_MK3(c1, c2, c3, u, v, ro, p);
				q2_1 = (a1 + b1 + c1) / 3.0;
				q2_2 = (a2 + b2 + c2) / 3.0;
				q3 = (a3 + b3 + c3) / 3.0;
			}
			else
			{
				K->Get_Sourse_MK1(a1, a2, a3, u, v, ro, p);
				q2_1 = a1;
				q2_2 = a2;
				q3 = a3;
			}

			double ro3, Q33, u3, v3, p3;

			if (K->type != C_centr)
			{
				/*ro3 = ro - T[now1] * (K->Potok[0] / Volume + ro * v / y);
				Q33 = Q - (T[now1] / Volume) * K->Potok[4] - T[now1] * Q * v / y;*/
				ro3 = ro - T[now1] * (K->Potok[0] / Volume);
				Q33 = Q - (T[now1] / Volume) * K->Potok[4] - T[now1];
				if (ro3 <= 0)
				{
					printf("Problemsssss  ro < 0!   %lf,  %lf\n", x, y);
					ro3 = 0.00001;
				}
				//u3 = (ro * u - T[now1] * (K->Potok[1] / Volume + ro * v * u / y - q2_1)) / ro3;
				//v3 = (ro * v - T[now1] * (K->Potok[2] / Volume + ro * v * v / y - q2_2)) / ro3;
				//p3 = (((p / (ggg - 1) + ro * (u * u + v * v) * 0.5) - T[now1] * (K->Potok[3] / Volume + //
				//	+v * (ggg * p / (ggg - 1) + ro * (u * u + v * v) * 0.5) / y - q3)) - //
				//	0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
				u3 = (ro * u - T[now1] * (K->Potok[1] / Volume - q2_1)) / ro3;
				v3 = (ro * v - T[now1] * ((K->Potok[2] - p * 0.5 * (vs1 + vs2)) / Volume - q2_2)) / ro3;
				p3 = (((p / (ggg - 1) + ro * (u * u + v * v) * 0.5) - T[now1] * (K->Potok[3] / Volume - q3)) - //
					0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
				if (p3 <= 0)
				{
					p3 = 0.000001;
				}

				K->par[now2].ro = ro3;
				K->par[now2].Q = Q33;
				K->par[now2].p = p3;
				K->par[now2].u = u3;
				K->par[now2].v = v3;
			}
			else
			{
				K->par[now2].ro = K->par[now1].ro;
				K->par[now2].Q = K->par[now1].Q;
				K->par[now2].p = K->par[now1].p;
				K->par[now2].u = K->par[now1].u;
				K->par[now2].v = K->par[now1].v;
			}


		}
	}
}

void Setka::Go_stationary_5_komponent_inner_MK(int step)
{
	//cout << "Setka.cpp    " << "START " << step << endl;
	int now1 = 1;
	int now2 = 0;
	double T[2];
	T[0] = T[1] = 0.00000001;
	mutex mut;
	double x_min = 0.0, y_min = 0.0;

	for (int st = 0; st < step; st++)
	{
		if (st % 10000 == 0 && st > 0)
		{
			cout << "Setka.cpp    " << st << " " << T[now2] << " " << x_min << " " << y_min << endl;
		}
		now1 = (now1 + 1) % 2; // Какие параметры сейчас берём
		now2 = (now2 + 1) % 2; // Какие параметры сейчас меняем
		double time = 10000000;
		T[now2] = 100000000;

		//omp_set_num_threads(4);
#pragma omp parallel for
		for (int num_cell = 0; num_cell < this->All_Cells_Inner.size(); num_cell++)
		{
			auto K = this->All_Cells_Inner[num_cell];
			K->Potok[0] = K->Potok[1] = K->Potok[2] = K->Potok[3] = K->Potok[4] = 0.0;

			double P[4];
			double PQ;
			Parametr par1 = K->par[now1];
			double n1, n2;
			double dist;
			double x, y;
			K->Get_Center(x, y);
			double radius = sqrt(kv(x) + kv(y));
			double Volume = K->Get_Volume();

			for (auto& i : K->Grans)
			{
				double x2, y2;
				double S = i->Get_square();
				i->Get_Center(x2, y2);
				dist = sqrt(kv(x - x2) + kv(y - y2));
				Parametr par2, par11, par3;
				par11 = par1;
				i->Get_par(par2, now1);
				i->Get_normal(n1, n2);

				double Vc, Vl, Vp;
				int met = 1;
				double dis = sqrt(kv(x2) + kv(y2));

				if (i->type == Usualy)
				{
					double x3, y3;
					i->Sosed->Get_Center(x3, y3);
					double dd = sqrt(kv(x3) + kv(y3));
					par11.ro = par11.ro * kv(radius) / kv(dis);
					par2.ro = par2.ro * kv(dd) / kv(dis);
					par11.Q = par11.Q * kv(radius) / kv(dis);
					par2.Q = par2.Q * kv(dd) / kv(dis);
					par11.p = par11.p * pow(radius / dis, 2 * ggg);
					par2.p = par2.p * pow(dd / dis, 2 * ggg);
					polar_perenos(x, y, x2, y2, par11.u, par11.v);
					polar_perenos(x3, y3, x2, y2, par2.u, par2.v);

					if (K->type == C_1 && i->Sosed->type == C_1)
					{
						Parametr pp1, pp2;
						i->Get_par_TVD_radial(pp1, now1);
						i->Gran_copy->Get_par_TVD_radial(pp2, now1);

						par11.ro = pp1.ro;
						par2.ro = pp2.ro;
						par11.p = pp1.p;
						par2.p = pp2.p;
					}
				}
				else if (i->type == Inner_sphere)
				{
					par11.ro = par2.ro;
					par11.Q = par2.Q;
					par11.p = par2.p;
					par11.u = par2.u;
					par11.v = par2.v;
				}
				else if (i->type == Axis)
				{
					par11.ro = par11.ro * kv(radius) / kv(dis);
					par11.Q = par11.Q * kv(radius) / kv(dis);
					par11.p = par11.p * pow(radius / dis, 2 * ggg);
					polar_perenos(x, y, x2, y2, par11.u, par11.v);
					par2 = par11;
					par2.v = -par2.v;

					/*i->Get_par_TVD(par2, now1);
					par11 = par2;
					polar_perenos(x, y, x2, y2, par11.u, par11.v);
					par2 = par11;
					par2.v = -par2.v;*/
				}

				if (K->type != C_centr)
				{
					time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro, par11.Q, par11.p, par11.u, par11.v, par2.ro, par2.Q, //
						par2.p, par2.u, par2.v, 0.0, P, PQ, n1, n2, dist, met, Vl, Vc, Vp));
					for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
					{
						K->Potok[k] = K->Potok[k] + P[k] * S;
					}
					K->Potok[4] = K->Potok[4] + PQ * S;
				}
			}

			if (time < T[now2])
			{
				mut.lock();
				T[now2] = time;
				x_min = x;
				y_min = y;
				mut.unlock();
			}

			/////////////////////////////////////////////////////////////////////////
			
			double a1, a2, a3;
			double q2_1, q2_2, q3;
			double u, v, ro, p, Q;

			/*u = par1.u * (chi_real / chi_);
			v = par1.v * (chi_real / chi_);
			ro = par1.ro / kv(chi_real / chi_);
			p = par1.p;
			Q = par1.Q / kv(chi_real / chi_);*/
			u = par1.u;
			v = par1.v;
			ro = par1.ro;
			p = par1.p;
			Q = par1.Q;

			K->Get_Sourse_MK1(a1, a2, a3, u, v, ro, p, false);
			q2_1 = a1;
			q2_2 = a2;
			q3 = a3; // / (chi_real / chi_);

			//q2_1 = q2_2 = q3 = 0.0;

			//cout << q2_1 << " " << q2_2 << " " << q3 << endl;

			u = par1.u;
			v = par1.v;
			ro = par1.ro;
			p = par1.p;
			Q = par1.Q;

			double ro3, Q33, u3, v3, p3;

			if (K->type != C_centr)
			{
				ro3 = ro - T[now1] * (K->Potok[0] / Volume + ro * v / y);
				Q33 = Q - (T[now1] / Volume) * K->Potok[4] - T[now1] * Q * v / y;
				if (ro3 <= 0)
				{
					printf("Problemsssss  ro < 0!   %lf,  %lf,  %lf,  %lf,  %lf,  %lf,  %lf, %lf,\n", x, y, ro3, ro, p, q2_1, q2_2, q3);
					ro3 = 0.0000001;
					exit(-1);
				}
				u3 = (ro * u - T[now1] * (K->Potok[1] / Volume + ro * v * u / y - q2_1)) / ro3;
				v3 = (ro * v - T[now1] * (K->Potok[2] / Volume + ro * v * v / y - q2_2)) / ro3;
				p3 = (((p / (ggg - 1) + ro * (u * u + v * v) * 0.5) - T[now1] * (K->Potok[3] / Volume + //
					+v * (ggg * p / (ggg - 1) + ro * (u * u + v * v) * 0.5) / y - q3)) - //
					0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
				if (p3 <= 0)
				{
					p3 = 0.000001;
				}

				K->par[now2].ro = ro3;
				K->par[now2].Q = Q33;
				K->par[now2].p = p3;
				K->par[now2].u = u3;
				K->par[now2].v = v3;
			}
			else
			{
				K->par[now2].ro = K->par[now1].ro;
				K->par[now2].Q = K->par[now1].Q;
				K->par[now2].p = K->par[now1].p;
				K->par[now2].u = K->par[now1].u;
				K->par[now2].v = K->par[now1].v;
			}


		}
	}
}

void Setka::Go(int step)
{
	cout << "Setka.cpp    " << "START " << step << endl;
	int now1 = 1;
	int now2 = 0;
	double T[2];
	T[0] = T[1] = 0.00000001;
	mutex mut;
	double xx, yy;

	vector <double> q1(5);
	vector <double> q2(5);
	vector <double> q(5);
	vector <double> n(3);
	double dsl, dsp;

	Solvers Sol;

	for (int st = 0; st < step; st++)
	{
		if (st % 1000 == 0 && st > 0)
		{
			cout << "Setka.cpp    " << st << " " << T[now2] << " " << xx << " " << yy << endl;
		}
		now1 = (now1 + 1) % 2; // Какие параметры сейчас берём
		now2 = (now2 + 1) % 2; // Какие параметры сейчас меняем
		double time = 10000000;
		T[now2] = 100000000;

		this->Move_surface(now1);
		this->Move_Setka_Calculate_2(T[now1]);

#pragma omp parallel for
		for (int num_cell = 0; num_cell < this->All_Cells.size(); num_cell++)
		{
			//cout << "Setka.cpp    " << "A1" << endl;
			auto K = this->All_Cells[num_cell];
			K->Potok[0] = K->Potok[1] = K->Potok[2] = K->Potok[3] = K->Potok[4] = 0.0;
			double P[4];
			double PQ;
			//cout << "Setka.cpp    " << "A2" << endl;
			Parametr par1 = K->par[now1];
			double n1, n2;
			double dist;
			double x, y;
			K->Get_Center(x, y);
			double W;
			//cout << "Setka.cpp    " << "A3" << endl;

			if (sqrt(kv(x) + kv(y)) <= 70.0)
			{
				K->par[1].u = K->par[0].u = 36.12 / (chi_real / chi_) * x / sqrt(x * x + y * y);       // Перенормировка
				K->par[1].v = K->par[0].v = 36.12 / (chi_real / chi_) * y / sqrt(x * x + y * y);
				K->par[1].ro = K->par[0].ro = 116.667 * kv(chi_real / chi_) / (x * x + y * y);
				K->par[1].p = K->par[0].p = kv(36.12 / (chi_real / chi_)) * (116.667 * kv(chi_real / chi_)) / (ggg * kv(10.0)) * pow(1.0 / sqrt(x * x + y * y), 2.0 * ggg);
				continue;
			}

			for (auto& i : K->Grans)
			{
				double x2, y2;
				double x3, y3;
				double S = i->Get_square();
				i->Get_Center(x2, y2);
				i->Get_Center_posle(x3, y3);
				dist = sqrt(kv(x - x2) + kv(y - y2));
				Parametr par2;
				i->Get_par(par2, now1);
				i->Get_normal(n1, n2);
				W = ((x3 - x2) * n1 + (y3 - y2) * n2) / T[now1];
				double Vc, Vl, Vp;


				if (i->A->type == P_Contact && i->B->type == P_Contact)
				{
					time = min(time, this->HLLC_2d_Korolkov_b_s(par1.ro, par1.Q, par1.p, par1.u, par1.v, par2.ro, par2.Q, //
						par2.p, par2.u, par2.v, W, P, PQ, n1, n2, dist, 1, Vl, Vc, Vp));
				}
				else
				{
					time = min(time, this->HLLC_2d_Korolkov_b_s(par1.ro, par1.Q, par1.p, par1.u, par1.v, par2.ro, par2.Q, //
						par2.p, par2.u, par2.v, W, P, PQ, n1, n2, dist, 0, Vl, Vc, Vp));
				}

				for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
				{
					K->Potok[k] = K->Potok[k] + P[k] * S;
				}
				//cout << "Setka.cpp    " << "B5" << endl;
				K->Potok[4] = K->Potok[4] + PQ * S;
				//cout << "Setka.cpp    " << "B6" << endl;
			}

			if (time < T[now2])
			{
				mut.lock();
				T[now2] = time;
				xx = x;
				yy = y;
				mut.unlock();
			}

			double ro3, p3, u3, v3, Q33;
			double Volume = K->Get_Volume();
			double Volume2 = K->Get_Volume_posle();

			ro3 = par1.ro * Volume/Volume2 - T[now1] * (K->Potok[0] / Volume2 + par1.ro * par1.v / y);
			Q33 = par1.Q * Volume / Volume2 - (T[now1] / Volume2) * K->Potok[4] - T[now1] * par1.Q * par1.v / y;
			if (ro3 <= 0)
			{
				printf("Problemsssss  par1.ro < 0!   %lf,  %lf\n", x, y);
				ro3 = 0.00001;
			}
			u3 = (par1.ro * par1.u * Volume / Volume2 - T[now1] * (K->Potok[1] / Volume2 + par1.ro * par1.v * par1.u / y)) / ro3;
			v3 = (par1.ro * par1.v * Volume / Volume2 - T[now1] * (K->Potok[2] / Volume2 + par1.ro * par1.v * par1.v / y)) / ro3;
			p3 = (((par1.p / (ggg - 1) + par1.ro * (par1.u * par1.u + par1.v * par1.v) * 0.5) * Volume / Volume2 - T[now1] * (K->Potok[3] / Volume2 +  //
				+par1.v * (ggg * par1.p / (ggg - 1) + par1.ro * (par1.u * par1.u + par1.v * par1.v) * 0.5) / y)) - //
				0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
			if (p3 <= 0)
			{
				p3 = 0.000001;
			}

			K->par[now2].ro = ro3;
			K->par[now2].Q = Q33;
			K->par[now2].p = p3;
			K->par[now2].u = u3;
			K->par[now2].v = v3;
		}

		for (auto& i : this->All_Points)
		{
			i->x = i->x2;
			i->y = i->y2;
			i->Vx = 0.0;
			i->Vy = 0.0;
			i->count = 0;
		}
	}
}

void Setka::Go_2(int step)
{
	//cout << "Setka.cpp    " << "START " << step << endl;
	int now1 = 1;
	int now2 = 0;
	double T[2];
	T[0] = T[1] = 0.00000001;
	mutex mut;

	for (int st = 0; st < step; st++)
	{
		//cout << "Setka.cpp    " << st << endl;
		if (st % 10000 == 0 && st > 0)
		{
			cout << "Setka.cpp    " << st << " " << T[now2] << endl;
		}
		now1 = (now1 + 1) % 2; // Какие параметры сейчас берём
		now2 = (now2 + 1) % 2; // Какие параметры сейчас меняем
		double time = 10000000;
		T[now2] = 100000000;

		this->Move_surface(now1, T[now1]);
		this->Move_Setka_Calculate_2(T[now1]);

		//omp_set_num_threads(4);
#pragma omp parallel for
		for (int num_cell = 0; num_cell < this->All_Cells.size(); num_cell++)
		{
			auto K = this->All_Cells[num_cell];
			K->Potok[0] = K->Potok[1] = K->Potok[2] = K->Potok[3] = K->Potok[4] = 0.0;

			double P[4];
			double PQ;
			Parametr par1 = K->par[now1];
			double n1, n2;
			double dist;
			double x, y;
			K->Get_Center(x, y);

			y = K->Get_Volume_rotate(360.0) / (2.0 * pi_ * K->Get_Volume());

			double radius = sqrt(kv(x) + kv(y));
			if (K->type == C_1 || K->type == C_centr)//(sqrt(kv(x) + kv(y)) <= 70.0)
			{
				K->par[1].u = K->par[0].u = 36.12 / (chi_real / chi_) * x / sqrt(x * x + y * y);       // Перенормировка
				K->par[1].v = K->par[0].v = 36.12 / (chi_real / chi_) * y / sqrt(x * x + y * y);
				K->par[1].ro = K->par[0].ro = 116.667 * kv(chi_real / chi_) / (x * x + y * y);
				K->par[1].p = K->par[0].p = kv(36.12 / (chi_real / chi_)) * (116.667 * kv(chi_real / chi_)) / (ggg * kv(10.0)) * pow(1.0 / sqrt(x * x + y * y), 2.0 * ggg);
				continue;
			}
			double Volume = K->Get_Volume(); // K->Get_Volume_rotate(alpha_rot);
			//cout << "Setka.cpp    " << "Vol = " << Volume << endl;
			double Volume2 = K->Get_Volume_posle(); // K->Get_Volume_posle_rotate(alpha_rot);
			//cout << "Setka.cpp    " << "Vol2 = " << Volume2 << endl;
			double W = 0.0;
			bool np = true;

			for (auto& i : K->Grans)
			{
				double x2, y2, x4, y4;
				//double S = i->Get_square_rotate(alpha_rot); 
				double S = i->Get_square();
				//cout << "Setka.cpp    " << "S = " << S << endl;
				i->Get_Center(x2, y2);
				i->Get_Center_posle(x4, y4);
				dist = sqrt(kv(x - x2) + kv(y - y2));
				Parametr par2, par11;
				par11 = par1;
				i->Get_par(par2, now1);
				i->Get_normal(n1, n2);
				double Vc, Vl, Vp;
				W = ((x4 - x2) * n1 + (y4 - y2) * n2) / T[now1];
				int met = 1;  // 1
				double dis = sqrt(kv(x2) + kv(y2));

				if (i->type == Usualy)
				{
					double x3, y3;
					i->Sosed->Get_Center(x3, y3);
					double dd = sqrt(kv(x3) + kv(y3));

					if (K->type == C_1)
					{
						par11.ro = par11.ro * kv(radius) / kv(dis);
						par11.Q = par11.Q * kv(radius) / kv(dis);
						par11.p = par11.p * pow(radius / dis, 2 * ggg);
						polar_perenos(x, y, x2, y2, par11.u, par11.v);

						Parametr pp1;
						i->Get_par_TVD_radial(pp1, now1);
						par11.ro = pp1.ro;
						par11.p = pp1.p;
					}
					else
					{
						i->Get_par_TVD(par11, now1);
					}

					if (i->Sosed->type == C_1)
					{
						par2.ro = par2.ro * kv(dd) / kv(dis);
						par2.Q = par2.Q * kv(dd) / kv(dis);
						par2.p = par2.p * pow(dd / dis, 2 * ggg);
						polar_perenos(x3, y3, x2, y2, par2.u, par2.v);

						Parametr pp2;
						i->Gran_copy->Get_par_TVD_radial(pp2, now1);
						par2.ro = pp2.ro;
						par2.p = pp2.p;
					}
					else
					{
						i->Gran_copy->Get_par_TVD(par2, now1);
					}
				}
				else if (i->type == Inner_sphere)
				{
					par11.ro = par11.ro * kv(radius) / kv(dis);
					par11.Q = par11.Q * kv(radius) / kv(dis);
					par11.p = par11.p * pow(radius / dis, 2 * ggg);
				}
				else if (K->type == C_1 && i->type == Axis)
				{
					par11.ro = par11.ro * kv(radius) / kv(dis);
					par11.Q = par11.Q * kv(radius) / kv(dis);
					par11.p = par11.p * pow(radius / dis, 2 * ggg);
					polar_perenos(x, y, x2, y2, par11.u, par11.v);
					//par11.v = 0.0;
					par2 = par11;
					par2.v = -par2.v;
				}
				else if (K->type != C_1 && i->type == Axis)
				{
					i->Get_par_TVD(par11, now1);
					par2 = par11;
					par11.v = 0.0;
					par2.v = 0.0;
					//par2.v = -par2.v;
				}


				/*if ((i->A->type == P_Contact && i->B->type != P_Contact)||(i->A->type != P_Contact && i->B->type == P_Contact))
				{
					met = 0;
				}*/

				np = false;
				bool god = false;

				if ((i->A->type == P_Contact && i->B->type == P_Contact))
				{
					god = true;
				}
				else if ((i->A->type == P_Inner_shock && i->B->type == P_Inner_shock))
				{
					god = true;
				}
				else if ((i->A->type == P_Outer_shock && i->B->type == P_Outer_shock))
				{
					god = true;
				}


				if ((i->A->type == P_Contact && i->B->type != P_Contact))
				{
					np = true;
				}
				else if ((i->A->type != P_Contact && i->B->type == P_Contact))
				{
					np = true;
				}
				else if ((i->A->type == P_Inner_shock && i->B->type != P_Inner_shock))
				{
					np = true;
				}
				else if ((i->A->type != P_Inner_shock && i->B->type == P_Inner_shock))
				{
					np = true;
				}
				else if ((i->A->type == P_Outer_shock && i->B->type != P_Outer_shock))
				{
					np = true;
				}
				else if ((i->A->type != P_Outer_shock && i->B->type == P_Outer_shock))
				{
					np = true;
				}

				//np = false;

				/*if (y < 30 && x > 120 && x < 160)
				{
					np = true;
				}*/

				if (god == true)
				{
					double Vl, Vp, VV;
					vector<double> qqq1(5);
					vector<double> qqq2(5);
					vector<double> qqq(5);
					vector<double> n(3);
					n[0] = n1;
					n[1] = n2;
					n[2] = 0.0;

					qqq1[0] = par11.ro;
					qqq1[1] = par11.p;
					qqq1[2] = par11.u;
					qqq1[3] = par11.v;
					qqq1[4] = 0.0;

					qqq2[0] = par2.ro;
					qqq2[1] = par2.p;
					qqq2[2] = par2.u;
					qqq2[3] = par2.v;
					qqq2[4] = 0.0;

					auto SS = Solvers();
					time = min(time, SS.Godunov_Solver_Alexashov(qqq1, qqq2, n, qqq, Vl, Vp, VV, W, dist));

					K->Potok[0] = K->Potok[0] + qqq[0] * S;
					K->Potok[1] = K->Potok[1] + qqq[1] * S;
					K->Potok[2] = K->Potok[2] + qqq[2] * S;
					K->Potok[3] = K->Potok[3] + qqq[4] * S;

					time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro, par11.Q, par11.p, par11.u, par11.v, par2.ro, par2.Q, //
						par2.p, par2.u, par2.v, W, P, PQ, n1, n2, dist, met, Vl, Vc, Vp, np));

					K->Potok[4] = K->Potok[4] + PQ * S;
				}
				else
				{
					time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro, par11.Q, par11.p, par11.u, par11.v, par2.ro, par2.Q, //
						par2.p, par2.u, par2.v, W, P, PQ, n1, n2, dist, met, Vl, Vc, Vp, np));
					for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
					{
						K->Potok[k] = K->Potok[k] + P[k] * S;
					}
					K->Potok[4] = K->Potok[4] + PQ * S;
				}

			}

			if (time < T[now2])
			{
				mut.lock();
				T[now2] = time;
				mut.unlock();
			}

			double ro3, Q33, u3, v3, p3;
			double u, v, ro, p, Q;

			u = par1.u;
			v = par1.v;
			ro = par1.ro;
			p = par1.p;
			Q = par1.Q;
			

			//ro3 = ro * Volume / Volume2 - T[now1] * (K->Potok[0] / Volume2);  // В декартовых
			ro3 = ro * Volume / Volume2 - T[now1] * (K->Potok[0] / Volume2 + ro * v / y);  // В цилиндрических

			//Q33 = Q * Volume / Volume2 - (T[now1] / Volume2) * K->Potok[4];
			Q33 = Q * Volume / Volume2 - (T[now1] / Volume2) * K->Potok[4] - T[now1] * Q * v / y;
			if (ro3 <= 0.0)
			{
				printf("Problemsssss  ro < 0!    %lf, %lf\n", x, y);
				ro3 = 0.00001;
			}
			//u3 = (ro * u * Volume / Volume2 - T[now1] * (K->Potok[1] / Volume2  - q2_1)) / ro3;
			//v3 = (ro * v * Volume / Volume2 - T[now1] * ( (K->Potok[2] - 2.0 * p * sin(0.5 * alpha_rot * pi_ / 180.0) * K->Get_Volume()) / Volume2  - q2_2)) / ro3;
			//p3 = (((p / (ggg - 1) + ro * (u * u + v * v) * 0.5) * Volume / Volume2 - T[now1] * (K->Potok[3] / Volume2 - q3)) - //
			//	0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
			u3 = (ro * u * Volume / Volume2 - T[now1] * (K->Potok[1] / Volume2 + ro * v * u / y )) / ro3;
			v3 = (ro * v * Volume / Volume2 - T[now1] * (K->Potok[2] / Volume2 + ro * v * v / y )) / ro3;
			p3 = (((p / (ggg - 1) + ro * (u * u + v * v) * 0.5) * Volume / Volume2 - T[now1] * (K->Potok[3] / Volume2 + //
				+v * (ggg * p / (ggg - 1) + ro * (u * u + v * v) * 0.5) / y)) - //
				0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
			if (p3 <= 0)
			{
				p3 = 0.000001;
			}

			K->par[now2].ro = ro3;
			K->par[now2].Q = Q33;
			K->par[now2].p = p3;
			K->par[now2].u = u3;
			K->par[now2].v = v3;

		}


		for (auto& i : this->All_Points)
		{
			i->x = i->x2;
			i->y = i->y2;
			i->Vx = 0.0;
			i->Vy = 0.0;
			i->count = 0;
		}
	}
}

void Setka::Go_3(int step)
// Идея в том, чтобы использовать расчёты в угле (3д задача, сведённая к 2д)
{
	//cout << "Setka.cpp    " << "START " << step << endl;
	int now1 = 1;
	int now2 = 0;
	double T[2];
	T[0] = T[1] = 0.00000001;
	mutex mut;

	for (int st = 0; st < step; st++)
	{
		//cout << "Setka.cpp    " << st << endl;
		if (st % 10000 == 0 && st > 0)
		{
			cout << "Setka.cpp    " << st << " " << T[now2] << endl;
		}
		now1 = (now1 + 1) % 2; // Какие параметры сейчас берём
		now2 = (now2 + 1) % 2; // Какие параметры сейчас меняем
		double time = 10000000;
		T[now2] = 100000000;

		this->Move_surface(now1, T[now1]);
		this->Move_Setka_Calculate_2(T[now1]);

		//omp_set_num_threads(4);
#pragma omp parallel for
		for (int num_cell = 0; num_cell < this->All_Cells.size(); num_cell++)
		{
			auto K = this->All_Cells[num_cell];
			K->Potok[0] = K->Potok[1] = K->Potok[2] = K->Potok[3] = K->Potok[4] = 0.0;
			K->Potok_H1[0] = K->Potok_H1[1] = K->Potok_H1[2] = K->Potok_H1[3] = 0.0;
			K->Potok_H2[0] = K->Potok_H2[1] = K->Potok_H2[2] = K->Potok_H2[3] = 0.0;
			K->Potok_H3[0] = K->Potok_H3[1] = K->Potok_H3[2] = K->Potok_H3[3] = 0.0;
			K->Potok_H4[0] = K->Potok_H4[1] = K->Potok_H4[2] = K->Potok_H4[3] = 0.0;

			double P[4];
			double PP[8];
			double PQ;
			Parametr par1 = K->par[now1];
			double n1, n2;
			double dist;
			double x, y;
			K->Get_Center(x, y);

			double vs1 = K->Get_Volume();
			double vs2 = K->Get_Volume_posle();

			y = K->Get_Volume_rotate(360.0) / (2.0 * pi_ * vs1);
			double y2 = K->Get_Volume_posle_rotate(360.0) / (2.0 * pi_ * vs2);

			double radius = sqrt(kv(x) + kv(y));
			if ((K->type == C_1 || K->type == C_centr)&& x > -350.0)//(sqrt(kv(x) + kv(y)) <= 70.0)
			{
				K->par[1].u = K->par[0].u = 36.12 / (chi_real / chi_) * x / sqrt(x * x + y * y);       // Перенормировка
				K->par[1].v = K->par[0].v = 36.12 / (chi_real / chi_) * y / sqrt(x * x + y * y);
				K->par[1].ro = K->par[0].ro = 116.667 * kv(chi_real / chi_) / (x * x + y * y);
				K->par[1].p = K->par[0].p = kv(36.12 / (chi_real / chi_)) * (116.667 * kv(chi_real / chi_)) / (ggg * kv(10.0)) * pow(1.0 / sqrt(x * x + y * y), 2.0 * ggg);
				continue;
			}



			double Volume = vs1 * y;
			//cout << "Setka.cpp    " << "Vol = " << Volume << endl;
			double Volume2 = vs2 * y2;
			//cout << "Setka.cpp    " << "Vol2 = " << Volume2 << endl;
			double W = 0.0;
			bool np = true;

			for (auto& i : K->Grans)
			{
				double x2, y2, x4, y4;
				double S = i->Get_square() * (i->A->y + i->B->y) * 0.5;
				//double S = i->Get_square();
				//cout << "Setka.cpp    " << "S = " << S << endl;
				i->Get_Center(x2, y2);
				i->Get_Center_posle(x4, y4);
				dist = sqrt(kv(x - x2) + kv(y - y2));
				Parametr par2, par11;
				par11 = par1;
				i->Get_par(par2, now1);
				i->Get_normal(n1, n2);
				double Vc, Vl, Vp;
				W = ((x4 - x2) * n1 + (y4 - y2) * n2) / T[now1];
				int met = 1;  // 1
				double dis = sqrt(kv(x2) + kv(y2));

				if (i->type == Usualy)
				{
					double x3, y3;
					i->Sosed->Get_Center(x3, y3);
					double dd = sqrt(kv(x3) + kv(y3));

					if (K->type == C_1)
					{
						par11.ro = par11.ro * kv(radius) / kv(dis);
						par11.ro_H1 = par11.ro_H1 * pow(radius / dis, H_pow);
						par11.Q = par11.Q * kv(radius) / kv(dis);
						par11.p = par11.p * pow(radius / dis, 2 * ggg);
						polar_perenos(x, y, x2, y2, par11.u, par11.v);
						polar_perenos(x, y, x2, y2, par11.u_H1, par11.v_H1);

						Parametr pp1;
						i->Get_par_TVD_radial(pp1, now1);
						par11.ro = pp1.ro;
						par11.p = pp1.p;
						//par11.u_H1 = pp1.u_H1;
						//par11.v_H1 = pp1.v_H1;
						par11.p_H1 = pp1.p_H1;
						par11.ro_H1 = pp1.ro_H1;
					}
					else
					{
						i->Get_par_TVD(par11, now1);
					}

					if (i->Sosed->type == C_1)
					{
						par2.ro = par2.ro * kv(dd) / kv(dis);
						par2.ro_H1 = par2.ro_H1 * pow(dd / dis, H_pow);
						par2.Q = par2.Q * kv(dd) / kv(dis);
						par2.p = par2.p * pow(dd / dis, 2 * ggg);
						polar_perenos(x3, y3, x2, y2, par2.u, par2.v);
						polar_perenos(x3, y3, x2, y2, par2.u_H1, par2.v_H1);

						Parametr pp2;
						i->Gran_copy->Get_par_TVD_radial(pp2, now1);
						par2.ro = pp2.ro;
						par2.p = pp2.p;
						//par2.u_H1 = pp2.u_H1;
						//par2.v_H1 = pp2.v_H1;
						par2.p_H1 = pp2.p_H1;
						par2.ro_H1 = pp2.ro_H1;
					}
					else
					{
						i->Gran_copy->Get_par_TVD(par2, now1);
					}
				}
				else if (i->type == Inner_sphere)
				{
					par11.ro = par11.ro * kv(radius) / kv(dis);
					par11.ro_H1 = par11.ro_H1 * pow(radius / dis, H_pow);
					par11.Q = par11.Q * kv(radius) / kv(dis);
					par11.p = par11.p * pow(radius / dis, 2 * ggg);
				}
				else if (K->type == C_1 && i->type == Axis)
				{
					par11.ro = par11.ro * kv(radius) / kv(dis);
					par11.ro_H1 = par11.ro_H1 * pow(radius / dis, H_pow);
					par11.Q = par11.Q * kv(radius) / kv(dis);
					par11.p = par11.p * pow(radius / dis, 2 * ggg);
					polar_perenos(x, y, x2, y2, par11.u, par11.v);
					polar_perenos(x, y, x2, y2, par11.u_H1, par11.v_H1);
					//par11.v = 0.0;
					par2 = par11;
					par2.v = 0.0;
					par11.v = 0.0;
				}
				else if (K->type != C_1 && i->type == Axis)
				{
					i->Get_par_TVD(par11, now1);
					par2 = par11;
					par2.v = 0.0;
					par11.v = 0.0;
				}


				/*if ((i->A->type == P_Contact && i->B->type != P_Contact)||(i->A->type != P_Contact && i->B->type == P_Contact))
				{
					met = 0;
				}*/

				np = false;
				bool god = false;

				if ((i->A->type == P_Contact && i->B->type == P_Contact))
				{
					god = true;
				}
				else if ((i->A->type == P_Inner_shock && i->B->type == P_Inner_shock))
				{
					god = true;
				}
				else if ((i->A->type == P_Outer_shock && i->B->type == P_Outer_shock))
				{
					god = true;
				}


				if ((i->A->type == P_Contact && i->B->type != P_Contact))
				{
					np = true;
				}
				else if ((i->A->type != P_Contact && i->B->type == P_Contact))
				{
					np = true;
				}
				else if ((i->A->type == P_Inner_shock && i->B->type != P_Inner_shock))
				{
					np = true;
				}
				else if ((i->A->type != P_Inner_shock && i->B->type == P_Inner_shock))
				{
					np = true;
				}
				else if ((i->A->type == P_Outer_shock && i->B->type != P_Outer_shock))
				{
					np = true;
				}
				else if ((i->A->type != P_Outer_shock && i->B->type == P_Outer_shock))
				{
					np = true;
				}

				//np = false;

				/*if (y < 30 && x > 120 && x < 160)
				{
					np = true;
				}*/

				if (god == true)
				{
					double Vl, Vp, VV;
					vector<double> qqq1(5);
					vector<double> qqq2(5);
					vector<double> qqq(5);
					vector<double> n(3);
					n[0] = n1;
					n[1] = n2;
					n[2] = 0.0;

					qqq1[0] = par11.ro;
					qqq1[1] = par11.p;
					qqq1[2] = par11.u;
					qqq1[3] = par11.v;
					qqq1[4] = 0.0;

					qqq2[0] = par2.ro;
					qqq2[1] = par2.p;
					qqq2[2] = par2.u;
					qqq2[3] = par2.v;
					qqq2[4] = 0.0;

					auto SS = Solvers();
					time = min(time, SS.Godunov_Solver_Alexashov(qqq1, qqq2, n, qqq, Vl, Vp, VV, W, dist));

					K->Potok[0] = K->Potok[0] + qqq[0] * S;
					K->Potok[1] = K->Potok[1] + qqq[1] * S;
					K->Potok[2] = K->Potok[2] + qqq[2] * S;
					K->Potok[3] = K->Potok[3] + qqq[4] * S;

					time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro, par11.Q, par11.p, par11.u, par11.v, par2.ro, par2.Q, //
						par2.p, par2.u, par2.v, W, P, PQ, n1, n2, dist, met, Vl, Vc, Vp, np));

					K->Potok[4] = K->Potok[4] + PQ * S;
				}
				else
				{
					time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro, par11.Q, par11.p, par11.u, par11.v, par2.ro, par2.Q, //
						par2.p, par2.u, par2.v, W, P, PQ, n1, n2, dist, met, Vl, Vc, Vp, np));
					for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
					{
						K->Potok[k] = K->Potok[k] + P[k] * S;
					}
					K->Potok[4] = K->Potok[4] + PQ * S;
				}

			}



			// Блок законов сохранения в ячейке

			if (time < T[now2])
			{
				mut.lock();
				T[now2] = time;
				mut.unlock();
			}


			double u, v, ro, p, Q;



			u = par1.u;
			v = par1.v;
			ro = par1.ro;
			p = par1.p;
			Q = par1.Q;

			double ro3, Q33, u3, v3, p3;

			ro3 = ro * Volume / Volume2 - T[now1] * (K->Potok[0] / Volume2);  // В декартовых
			Q33 = Q * Volume / Volume2 - (T[now1] / Volume2) * K->Potok[4];
			if (ro3 <= 0.0)
			{
				printf("Problemsssss  ro < 0!    %lf, %lf\n", x, y);
				ro3 = 0.00001;
			}
			u3 = (ro * u * Volume / Volume2 - T[now1] * (K->Potok[1] / Volume2 )) / ro3;
			//v3 = (ro * v * Volume / Volume2 - T[now1] * ((K->Potok[2] - 2.0 * p * sin(0.5 * alpha_rot * pi_ / 180.0) * K->Get_Volume()) / Volume2)) / ro3;
			v3 = (ro * v * Volume / Volume2 - T[now1] * ((K->Potok[2] - p * 0.5 * (vs1 + vs2)) / Volume2)) / ro3;
			p3 = (((p / (ggg - 1) + ro * (u * u + v * v) * 0.5) * Volume / Volume2 - T[now1] * (K->Potok[3] / Volume2 )) - //
				0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
			if (p3 <= 0)
			{
				p3 = 0.000001;
			}

			K->par[now2].ro = ro3;
			K->par[now2].Q = Q33;
			K->par[now2].p = p3;
			K->par[now2].u = u3;
			K->par[now2].v = v3;

		}


		for (auto& i : this->All_Points)
		{
			i->x = i->x2;
			i->y = i->y2;
			i->Vx = 0.0;
			i->Vy = 0.0;
			i->count = 0;
		}
	}
}

void Setka::Go_5_komponent(int step)
{
	//cout << "Setka.cpp    " << "START " << step << endl;
	int now1 = 1;
	int now2 = 0;
	double T[2];
	T[0] = T[1] = 0.00000001;
	mutex mut;

	for (int st = 0; st < step; st++)
	{
		//cout << "Setka.cpp    " << st << endl;
		if (st % 10000 == 0 && st > 0)
		{
			cout << "Setka.cpp    " << st << " " << T[now2] << endl;
		}
		now1 = (now1 + 1) % 2; // Какие параметры сейчас берём
		now2 = (now2 + 1) % 2; // Какие параметры сейчас меняем
		double time = 10000000;
		T[now2] = 100000000;

		this->Move_surface(now1, T[now1]);
		this->Move_Setka_Calculate(T[now1]);
		
		//omp_set_num_threads(4);
#pragma omp parallel for
		for (int num_cell = 0; num_cell < this->All_Cells.size(); num_cell++)
		{
			auto K = this->All_Cells[num_cell];
			K->Potok[0] = K->Potok[1] = K->Potok[2] = K->Potok[3] = K->Potok[4] = 0.0;
			K->Potok_H1[0] = K->Potok_H1[1] = K->Potok_H1[2] = K->Potok_H1[3] = 0.0;
			K->Potok_H2[0] = K->Potok_H2[1] = K->Potok_H2[2] = K->Potok_H2[3] = 0.0;
			K->Potok_H3[0] = K->Potok_H3[1] = K->Potok_H3[2] = K->Potok_H3[3] = 0.0;
			K->Potok_H4[0] = K->Potok_H4[1] = K->Potok_H4[2] = K->Potok_H4[3] = 0.0;

			double P[4];
			double PQ;
			Parametr par1 = K->par[now1];
			double n1, n2;
			double dist;
			double x, y;
			K->Get_Center(x, y);

			//y = K->Get_Volume_rotate(5.0) / (2.0 * pi_ * K->Get_Volume() / 72.0);

			double radius = sqrt(kv(x) + kv(y));
			if (radius < R11_)
			{
				continue;
			}
			double Volume = K->Get_Volume(); // K->Get_Volume_rotate(alpha_rot);
			//cout << "Setka.cpp    " << "Vol = " << Volume << endl;
			double Volume2 = K->Get_Volume_posle(); // K->Get_Volume_posle_rotate(alpha_rot);
			//cout << "Setka.cpp    " << "Vol2 = " << Volume2 << endl;
			double W = 0.0;
			bool np = true;

			for (auto& i : K->Grans)
			{
				double x2, y2, x4, y4;
				//double S = i->Get_square_rotate(alpha_rot); 
				double S = i->Get_square();
				//cout << "Setka.cpp    " << "S = " << S << endl;
				i->Get_Center(x2, y2);
				i->Get_Center_posle(x4, y4);
				dist = sqrt(kv(x - x2) + kv(y - y2));
				Parametr par2, par11;
				par11 = par1;
				i->Get_par(par2, now1);
				i->Get_normal(n1, n2);
				double Vc, Vl, Vp;
				W = ((x4 - x2) * n1 + (y4 - y2) * n2) / T[now1];
				int met = 1;  // 1
				double dis = sqrt(kv(x2) + kv(y2));

				if (i->type == Usualy)
				{
					double x3, y3;
					i->Sosed->Get_Center(x3, y3);
					double dd = sqrt(kv(x3) + kv(y3));

					if (K->type == C_1)
					{
						par11.ro = par11.ro * kv(radius) / kv(dis);
						par11.ro_H1 = par11.ro_H1 * pow(radius / dis, H_pow);
						par11.Q = par11.Q * kv(radius) / kv(dis);
						par11.p = par11.p * pow(radius / dis, 2 * ggg);
						polar_perenos(x, y, x2, y2, par11.u, par11.v);
						polar_perenos(x, y, x2, y2, par11.u_H1, par11.v_H1);

						Parametr pp1;
						i->Get_par_TVD_radial(pp1, now1);
						par11.ro = pp1.ro;
						par11.p = pp1.p;
						//par11.u_H1 = pp1.u_H1;
						//par11.v_H1 = pp1.v_H1;
						par11.p_H1 = pp1.p_H1;
						par11.ro_H1 = pp1.ro_H1;
					}
					else
					{
						i->Get_par_TVD(par11, now1);
					}

					if (i->Sosed->type == C_1)
					{
						par2.ro = par2.ro * kv(dd) / kv(dis);
						par2.ro_H1 = par2.ro_H1 * pow(dd / dis, H_pow);
						par2.Q = par2.Q * kv(dd) / kv(dis);
						par2.p = par2.p * pow(dd / dis, 2 * ggg);
						polar_perenos(x3, y3, x2, y2, par2.u, par2.v);
						polar_perenos(x3, y3, x2, y2, par2.u_H1, par2.v_H1);

						Parametr pp2;
						i->Gran_copy->Get_par_TVD_radial(pp2, now1);
						par2.ro = pp2.ro;
						par2.p = pp2.p;
						//par2.u_H1 = pp2.u_H1;
						//par2.v_H1 = pp2.v_H1;
						par2.p_H1 = pp2.p_H1;
						par2.ro_H1 = pp2.ro_H1;
					}
					else
					{
						i->Gran_copy->Get_par_TVD(par2, now1);
					}
				}
				else if (i->type == Inner_sphere)
				{
					par11.ro = par11.ro * kv(radius) / kv(dis);
					par11.ro_H1 = par11.ro_H1 * pow(radius / dis, H_pow);
					par11.Q = par11.Q * kv(radius) / kv(dis);
					par11.p = par11.p * pow(radius / dis, 2 * ggg);
				}
				else if (K->type == C_1 && i->type == Axis)
				{
					par11.ro = par11.ro * kv(radius) / kv(dis);
					par11.ro_H1 = par11.ro_H1 * pow(radius / dis, H_pow);
					par11.Q = par11.Q * kv(radius) / kv(dis);
					par11.p = par11.p * pow(radius / dis, 2 * ggg);
					polar_perenos(x, y, x2, y2, par11.u, par11.v);
					polar_perenos(x, y, x2, y2, par11.u_H1, par11.v_H1);
					//par11.v = 0.0;
					par2 = par11;
					par2.v = -par2.v;
					par2.v_H1 = -par2.v_H1;
					par2.v_H2 = -par2.v_H2;
					par2.v_H3 = -par2.v_H3;
					par2.v_H4 = -par2.v_H4;
				}
				else if (K->type != C_1 && i->type == Axis)
				{
					i->Get_par_TVD(par11, now1);
					par2 = par11;
					par2.v = -par2.v;
					par2.v_H1 = -par2.v_H1;
					par2.v_H2 = -par2.v_H2;
					par2.v_H3 = -par2.v_H3;
					par2.v_H4 = -par2.v_H4;
				}
				

				/*if ((i->A->type == P_Contact && i->B->type != P_Contact)||(i->A->type != P_Contact && i->B->type == P_Contact))
				{
					met = 0;
				}*/

				np = false;
				bool god = false;

				if ((i->A->type == P_Contact && i->B->type == P_Contact))
				{
					god = true;
				}
				else if ((i->A->type == P_Inner_shock && i->B->type == P_Inner_shock))
				{
					god = true;
				}
				else if ((i->A->type == P_Outer_shock && i->B->type == P_Outer_shock))
				{
					god = true;
				}
				
				
				if ((i->A->type == P_Contact && i->B->type != P_Contact))
				{
					np = true;
				}
				else if ((i->A->type != P_Contact && i->B->type == P_Contact))
				{
					np = true;
				}
				else if ((i->A->type == P_Inner_shock && i->B->type != P_Inner_shock))
				{
					np = true;
				}
				else if ((i->A->type != P_Inner_shock && i->B->type == P_Inner_shock))
				{
					np = true;
				}
				else if ((i->A->type == P_Outer_shock && i->B->type != P_Outer_shock))
				{
					np = true;
				}
				else if ((i->A->type != P_Outer_shock && i->B->type == P_Outer_shock))
				{
					np = true;
				}

				//np = false;

				/*if (y < 30 && x > 120 && x < 160)
				{
					np = true;
				}*/

				if (K->type != C_centr)
				{
					time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro_H1, par11.Q, par11.p_H1, par11.u_H1, par11.v_H1, par2.ro_H1, par2.Q, //
						par2.p_H1, par2.u_H1, par2.v_H1, W, P, PQ, n1, n2, dist, met, Vl, Vc, Vp));
					for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
					{
						K->Potok_H1[k] = K->Potok_H1[k] + P[k] * S;
					}
				}

				time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro_H2, par11.Q, par11.p_H2, par11.u_H2, par11.v_H2, par2.ro_H2, par2.Q, //
					par2.p_H2, par2.u_H2, par2.v_H2, W, P, PQ, n1, n2, dist, met, Vl, Vc, Vp));
				for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
				{
					K->Potok_H2[k] = K->Potok_H2[k] + P[k] * S;
				}

				time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro_H3, par11.Q, par11.p_H3, par11.u_H3, par11.v_H3, par2.ro_H3, par2.Q, //
					par2.p_H3, par2.u_H3, par2.v_H3, W, P, PQ, n1, n2, dist, met, Vl, Vc, Vp));
				for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
				{
					K->Potok_H3[k] = K->Potok_H3[k] + P[k] * S;
				}

				time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro_H4, par11.Q, par11.p_H4, par11.u_H4, par11.v_H4, par2.ro_H4, par2.Q, //
					par2.p_H4, par2.u_H4, par2.v_H4, W, P, PQ, n1, n2, dist, met, Vl, Vc, Vp));
				for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
				{
					K->Potok_H4[k] = K->Potok_H4[k] + P[k] * S;
				}

				if (god == true)
				{
					double Vl, Vp, VV;
					vector<double> qqq1(5);
					vector<double> qqq2(5);
					vector<double> qqq(5);
					vector<double> n(3);
					n[0] = n1;
					n[1] = n2;
					n[2] = 0.0;

					qqq1[0] = par11.ro;
					qqq1[1] = par11.p;
					qqq1[2] = par11.u;
					qqq1[3] = par11.v;
					qqq1[4] = 0.0;

					qqq2[0] = par2.ro;
					qqq2[1] = par2.p;
					qqq2[2] = par2.u;
					qqq2[3] = par2.v;
					qqq2[4] = 0.0;

					auto SS = Solvers();
					time = min(time, SS.Godunov_Solver_Alexashov(qqq1, qqq2, n, qqq, Vl, Vp, VV, W, dist));

					K->Potok[0] = K->Potok[0] + qqq[0] * S;
					K->Potok[1] = K->Potok[1] + qqq[1] * S;
					K->Potok[2] = K->Potok[2] + qqq[2] * S;
					K->Potok[3] = K->Potok[3] + qqq[4] * S;

					time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro, par11.Q, par11.p, par11.u, par11.v, par2.ro, par2.Q, //
						par2.p, par2.u, par2.v, W, P, PQ, n1, n2, dist, met, Vl, Vc, Vp, np));

					K->Potok[4] = K->Potok[4] + PQ * S;
				}
				else
				{
					time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro, par11.Q, par11.p, par11.u, par11.v, par2.ro, par2.Q, //
						par2.p, par2.u, par2.v, W, P, PQ, n1, n2, dist, met, Vl, Vc, Vp, np));
					for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
					{
						K->Potok[k] = K->Potok[k] + P[k] * S;
					}
					K->Potok[4] = K->Potok[4] + PQ * S;
				}
				
			}

			if (time < T[now2])
			{
				mut.lock();
				T[now2] = time;
				mut.unlock();
			}



			double U_M_H1, U_M_H2, U_M_H3, U_M_H4;
			double U_H1, U_H2, U_H3, U_H4;
			double sigma_H1, sigma_H2, sigma_H3, sigma_H4;
			double nu_H1, nu_H2, nu_H3, nu_H4;
			double q2_1, q2_2, q3;
			double u, v, ro, p, Q;
			if (false)//(K->type == C_centr || K->type == C_1 || K->type == C_2 || K->type == C_3)//(par1.Q / par1.ro < 90)//
			{
				u = par1.u * (chi_real/chi_);
				v = par1.v * (chi_real / chi_);
				ro = par1.ro / kv(chi_real / chi_);
				p = par1.p;
				Q = par1.Q / kv(chi_real / chi_);
			}
			else
			{
				u = par1.u;
				v = par1.v;
				ro = par1.ro;
				p = par1.p;
				Q = par1.Q;
			}
		
			
			double u_H1 = par1.u_H1, v_H1 = par1.v_H1, ro_H1 = par1.ro_H1, p_H1 = par1.p_H1;
			double u_H2 = par1.u_H2, v_H2 = par1.v_H2, ro_H2 = par1.ro_H2, p_H2 = par1.p_H2;
			double u_H3 = par1.u_H3, v_H3 = par1.v_H3, ro_H3 = par1.ro_H3, p_H3 = par1.p_H3;
			double u_H4 = par1.u_H4, v_H4 = par1.v_H4, ro_H4 = par1.ro_H4, p_H4 = par1.p_H4;

			U_M_H1 = sqrt(kv(u - u_H1) + kv(v - v_H1) + (64.0 / (9.0 * pi)) //
				* (p / ro + 2.0 * p_H1 / ro_H1));
			U_M_H2 = sqrt(kv(u - u_H2) + kv(v - v_H2) + (64.0 / (9.0 * pi)) //
				* (p / ro + 2.0 * p_H2 / ro_H2));
			U_M_H3 = sqrt(kv(u - u_H3) + kv(v - v_H3) + (64.0 / (9.0 * pi)) //
				* (p / ro + 2.0 * p_H3 / ro_H3));
			U_M_H4 = sqrt(kv(u - u_H4) + kv(v - v_H4) + (64.0 / (9.0 * pi)) //
				* (p / ro + 2.0 * p_H4 / ro_H4));

			U_H1 = sqrt(kv(u - u_H1) + kv(v - v_H1) + (4.0 / pi) //
				* (p / ro + 2.0 * p_H1 / ro_H1));
			U_H2 = sqrt(kv(u - u_H2) + kv(v - v_H2) + (4.0 / pi) //
				* (p / ro + 2.0 * p_H2 / ro_H2));
			U_H3 = sqrt(kv(u - u_H3) + kv(v - v_H3) + (4.0 / pi) //
				* (p / ro + 2.0 * p_H3 / ro_H3));
			U_H4 = sqrt(kv(u - u_H4) + kv(v - v_H4) + (4.0 / pi) //
				* (p / ro + 2.0 * p_H4 / ro_H4));

			sigma_H1 = kv(1.0 - a_2 * log(U_M_H1)); // 0.1243
			sigma_H2 = kv(1.0 - a_2 * log(U_M_H2));
			sigma_H3 = kv(1.0 - a_2 * log(U_M_H3)); // 0.1121     a_2
			sigma_H4 = kv(1.0 - a_2 * log(U_M_H4));

			nu_H1 = ro * ro_H1 * U_M_H1 * sigma_H1;
			nu_H2 = ro * ro_H2 * U_M_H2 * sigma_H2;
			nu_H3 = ro * ro_H3 * U_M_H3 * sigma_H3;
			nu_H4 = ro * ro_H4 * U_M_H4 * sigma_H4;

			q2_1 = (n_p_LISM_ / Kn_) * (nu_H1 * (u_H1 - u) + nu_H2 * (u_H2 - u) //
				+ nu_H3 * (u_H3 - u) + nu_H4 * (u_H4 - u));
			q2_2 = (n_p_LISM_ / Kn_) * (nu_H1 * (v_H1 - v) + nu_H2 * (v_H2 - v) //
				+ nu_H3 * (v_H3 - v) + nu_H4 * (v_H4 - v));

			if (false)//(K->type == C_centr || K->type == C_1 || K->type == C_2 || K->type == C_3)//(par1.Q / par1.ro < 90)//
			{
				q3 = (n_p_LISM_ / Kn_) * (nu_H1 * ((kv(u_H1) + kv(v_H1) - kv(u) - kv(v)) / 2.0 + //
					(U_H1 / U_M_H1) * (2.0 * p_H1 / ro_H1 - p / ro)) + //
					nu_H2 * ((kv(u_H2) + kv(v_H2) - kv(u) - kv(v)) / 2.0 + //
						(U_H2 / U_M_H2) * (2.0 * p_H2 / ro_H2 - p / ro)) + //
					nu_H3 * ((kv(u_H3) + kv(v_H3) - kv(u) - kv(v)) / 2.0 + //
						(U_H3 / U_M_H3) * (2.0 * p_H3 / ro_H3 - p / ro)) + //
					nu_H4 * ((kv(u_H4) + kv(v_H4) - kv(u) - kv(v)) / 2.0 + //
						(U_H4 / U_M_H4) * (2.0 * p_H4 / ro_H4 - p / ro))) / (chi_real / chi_);
			}
			else
			{
				q3 = (n_p_LISM_ / Kn_) * (nu_H1 * ((kv(u_H1) + kv(v_H1) - kv(u) - kv(v)) / 2.0 + //
					(U_H1 / U_M_H1) * (2.0 * p_H1 / ro_H1 - p / ro)) + //
					nu_H2 * ((kv(u_H2) + kv(v_H2) - kv(u) - kv(v)) / 2.0 + //
						(U_H2 / U_M_H2) * (2.0 * p_H2 / ro_H2 - p / ro)) + //
					nu_H3 * ((kv(u_H3) + kv(v_H3) - kv(u) - kv(v)) / 2.0 + //
						(U_H3 / U_M_H3) * (2.0 * p_H3 / ro_H3 - p / ro)) + //
					nu_H4 * ((kv(u_H4) + kv(v_H4) - kv(u) - kv(v)) / 2.0 + //
						(U_H4 / U_M_H4) * (2.0 * p_H4 / ro_H4 - p / ro)));
			}


			/*q2_1 = 0.0;
			q2_2 = 0.0;
			q3 = 0.0;*/

			double ro2_p, p2_p, V1_p, V2_p, QQ2;
			double ro3, p3, u3, v3, Q33;
			double kappa = Q / ro;

			
			u = par1.u;
			v = par1.v;
			ro = par1.ro;
			p = par1.p;
			Q = par1.Q;

			if (K->type != C_centr)
			{
				//ro3 = ro * Volume / Volume2 - T[now1] * (K->Potok[0] / Volume2);  // В декартовых
				ro3 = ro * Volume / Volume2 - T[now1] * (K->Potok[0] / Volume2 + ro * v / y);  // В цилиндрических

				//Q33 = Q * Volume / Volume2 - (T[now1] / Volume2) * K->Potok[4];
				Q33 = Q * Volume / Volume2 - (T[now1] / Volume2) * K->Potok[4] - T[now1] * Q * v / y;
				if (ro3 <= 0.0)
				{
					printf("Problemsssss  ro < 0!    %lf, %lf\n", x, y);
					ro3 = 0.00001;
				}
				//u3 = (ro * u * Volume / Volume2 - T[now1] * (K->Potok[1] / Volume2  - q2_1)) / ro3;
				//v3 = (ro * v * Volume / Volume2 - T[now1] * ( (K->Potok[2] - 2.0 * p * sin(0.5 * alpha_rot * pi_ / 180.0) * K->Get_Volume()) / Volume2  - q2_2)) / ro3;
				//p3 = (((p / (ggg - 1) + ro * (u * u + v * v) * 0.5) * Volume / Volume2 - T[now1] * (K->Potok[3] / Volume2 - q3)) - //
				//	0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
				u3 = (ro * u * Volume / Volume2 - T[now1] * (K->Potok[1] / Volume2 + ro * v * u / y - q2_1)) / ro3;
				v3 = (ro * v * Volume / Volume2 - T[now1] * (K->Potok[2] / Volume2 + ro * v * v / y - q2_2)) / ro3;
				p3 = (((p / (ggg - 1) + ro * (u * u + v * v) * 0.5) * Volume / Volume2 - T[now1] * (K->Potok[3] / Volume2 + //
					+v * (ggg * p / (ggg - 1) + ro * (u * u + v * v) * 0.5) / y - q3)) - //
					0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
				if (p3 <= 0)
				{
					p3 = 0.000001;
				}

				K->par[now2].ro = ro3;
				K->par[now2].Q = Q33;
				K->par[now2].p = p3;
				K->par[now2].u = u3;
				K->par[now2].v = v3;
			}
			else
			{
				K->par[now2].ro = K->par[now1].ro;
				K->par[now2].Q = K->par[now1].Q;
				K->par[now2].p = K->par[now1].p;
				K->par[now2].u = K->par[now1].u;
				K->par[now2].v = K->par[now1].v;
			}

			if (false)//(K->type == C_centr || K->type == C_1 || K->type == C_2 || K->type == C_3)//(par1.Q / par1.ro < 90)//
			{
				u = par1.u * (chi_real / chi_);
				v = par1.v * (chi_real / chi_);
				ro = par1.ro / kv(chi_real / chi_);
				p = par1.p;
				Q = par1.Q / kv(chi_real / chi_);
			}
			else
			{
				u = par1.u;
				v = par1.v;
				ro = par1.ro;
				p = par1.p;
				Q = par1.Q;
			}


			int k1 = 0, k2 = 0, k3 = 0, k4 = 0;

			if (K->type == C_centr || K->type == C_1 || K->type == C_2 || K->type == C_3)//(par1.Q / par1.ro < 90)//
			{
				k3 = 0;
				k4 = 0;
				if (K->type == C_1 || K->type == C_centr)//(((kv(u) + kv(v)) / (ggg * p / ro) > 1.0) && ((radius <= 250.0)))
				{
					k1 = 1;
					k2 = 0;
				}
				else
				{
					k1 = 0;
					k2 = 1;
				}
			}
			else
			{
				k1 = 0;
				k2 = 0;
				if (K->type == C_4)//(p / ro > 1.8)   // ( ((y < 8.0)&&(p / ro > (y * (-0.0238) + 0.36)))||( (y >= 8.0)&&(p / ro > 0.17) ) )
				{
					k4 = 0;
					k3 = 1;
				}
				else
				{
					k3 = 0;
					k4 = 1;
				}
			}

			double S1, S2;
			double q1_H1, q1_H2, q1_H3, q1_H4;
			double q21_H1, q21_H2, q21_H3, q21_H4;
			double q22_H1, q22_H2, q22_H3, q22_H4;
			double q3_H1, q3_H2, q3_H3, q3_H4;

			S1 = nu_H1 + nu_H2 + nu_H3 + nu_H4;
			S2 = nu_H1 * ((kv(u) + kv(v)) / 2.0 + (U_H1 / U_M_H1) * (p / ro))//
				+ nu_H2 * ((kv(u) + kv(v)) / 2.0 + (U_H2 / U_M_H2) * (p / ro))//
				+ nu_H3 * ((kv(u) + kv(v)) / 2.0 + (U_H3 / U_M_H3) * (p / ro))//
				+ nu_H4 * ((kv(u) + kv(v)) / 2.0 + (U_H4 / U_M_H4) * (p / ro));

			q1_H1 = (n_H_LISM_ / Kn_) * (k1 * S1 - nu_H1);
			q1_H2 = (n_H_LISM_ / Kn_) * (k2 * S1 - nu_H2);
			q1_H3 = (n_H_LISM_ / Kn_) * (k3 * S1 - nu_H3);
			q1_H4 = (n_H_LISM_ / Kn_) * (k4 * S1 - nu_H4);

			q21_H1 = (n_H_LISM_ / Kn_) * (k1 * S1 * u - nu_H1 * u_H1);
			q21_H2 = (n_H_LISM_ / Kn_) * (k2 * S1 * u - nu_H2 * u_H2);
			q21_H3 = (n_H_LISM_ / Kn_) * (k3 * S1 * u - nu_H3 * u_H3);
			q21_H4 = (n_H_LISM_ / Kn_) * (k4 * S1 * u - nu_H4 * u_H4);

			q22_H1 = (n_H_LISM_ / Kn_) * (k1 * S1 * v - nu_H1 * v_H1);
			q22_H2 = (n_H_LISM_ / Kn_) * (k2 * S1 * v - nu_H2 * v_H2);
			q22_H3 = (n_H_LISM_ / Kn_) * (k3 * S1 * v - nu_H3 * v_H3);
			q22_H4 = (n_H_LISM_ / Kn_) * (k4 * S1 * v - nu_H4 * v_H4);

			q3_H1 = (n_H_LISM_ / Kn_) * (k1 * S2 - nu_H1 * ((kv(u_H1) + kv(v_H1)) / 2.0 + (U_H1 / U_M_H1) * 2.0 * (p_H1 / ro_H1)));
			q3_H2 = (n_H_LISM_ / Kn_) * (k2 * S2 - nu_H2 * ((kv(u_H2) + kv(v_H2)) / 2.0 + (U_H2 / U_M_H2) * 2.0 * (p_H2 / ro_H2)));
			q3_H3 = (n_H_LISM_ / Kn_) * (k3 * S2 - nu_H3 * ((kv(u_H3) + kv(v_H3)) / 2.0 + (U_H3 / U_M_H3) * 2.0 * (p_H3 / ro_H3)));
			q3_H4 = (n_H_LISM_ / Kn_) * (k4 * S2 - nu_H4 * ((kv(u_H4) + kv(v_H4)) / 2.0 + (U_H4 / U_M_H4) * 2.0 * (p_H4 / ro_H4)));


			if (K->type != C_centr)
			{
				ro3 = ro_H1 * Volume / Volume2 - T[now1] * (K->Potok_H1[0] / Volume2 + ro_H1 * v_H1 / y - q1_H1);
				//ro3 = ro_H1 * Volume / Volume2 - T[now1] * (K->Potok_H1[0] / Volume2 - q1_H1);
				if (ro3 <= 0)
				{
					printf("Problemsssss  ro H1 < 0!   %lf,  %lf\n", x, y);
					ro3 = 0.0000001;
				}
				//u3 = (ro_H1 * u_H1 * Volume / Volume2 - T[now1] * (K->Potok_H1[1] / Volume2  - q21_H1)) / ro3;
				//v3 = (ro_H1 * v_H1 * Volume / Volume2 - T[now1] * ((K->Potok_H1[2] - 2.0 * p_H1 * sin(0.5 * alpha_rot * pi_ / 180.0) * K->Get_Volume()) / Volume2  - q22_H1)) / ro3;
				//p3 = (((p_H1 / (ggg - 1) + ro_H1 * (u_H1 * u_H1 + v_H1 * v_H1) * 0.5) * Volume / Volume2 - T[now1] * (K->Potok_H1[3] / Volume2  - q3_H1)) - //
				//	0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
				u3 = (ro_H1 * u_H1 * Volume / Volume2 - T[now1] * (K->Potok_H1[1] / Volume2 + ro_H1 * v_H1 * u_H1 / y - q21_H1)) / ro3;
				v3 = (ro_H1 * v_H1 * Volume / Volume2 - T[now1] * (K->Potok_H1[2] / Volume2 + ro_H1 * v_H1 * v_H1 / y - q22_H1)) / ro3;
				p3 = (((p_H1 / (ggg - 1) + ro_H1 * (u_H1 * u_H1 + v_H1 * v_H1) * 0.5) * Volume / Volume2 - T[now1] * (K->Potok_H1[3] / Volume2 + //
					+v_H1 * (ggg * p_H1 / (ggg - 1) + ro_H1 * (u_H1 * u_H1 + v_H1 * v_H1) * 0.5) / y - q3_H1)) - //
					0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
				if (p3 <= 0)
				{
					p3 = 0.000001;
				}

				K->par[now2].ro_H1 = ro3;
				K->par[now2].p_H1 = p3;
				K->par[now2].u_H1 = u3;
				K->par[now2].v_H1 = v3;

				/*if (x > 600.0)
				{
					K->par[now2].ro_H1 = 1.0;
					K->par[now2].p_H1 = 1.0;
					K->par[now2].u_H1 = Velosity_inf;
					K->par[now2].v_H1 = 0.0;
				}*/
			}
			else
			{
				K->par[now2].ro_H1 = K->par[now1].ro_H1;
				K->par[now2].p_H1 = K->par[now1].p_H1;
				K->par[now2].u_H1 = K->par[now1].u_H1;
				K->par[now2].v_H1 = K->par[now1].v_H1;
			}


			ro3 = ro_H2 * Volume / Volume2 - T[now1] * (K->Potok_H2[0] / Volume2 + ro_H2 * v_H2 / y - q1_H2);
			//ro3 = ro_H2 * Volume / Volume2 - T[now1] * (K->Potok_H2[0] / Volume2 - q1_H2);
			if (ro3 <= 0)
			{
				printf("Problemsssss  ro H2 < 0!  %lf,  %lf\n", x, y);
				ro3 = 0.000001;
			}
			//u3 = (ro_H2 * u_H2 * Volume / Volume2 - T[now1] * (K->Potok_H2[1] / Volume2 - q21_H2)) / ro3;
			//v3 = (ro_H2 * v_H2 * Volume / Volume2 - T[now1] * ((K->Potok_H2[2] - 2.0 * p_H2 * sin(0.5 * alpha_rot * pi_ / 180.0) * K->Get_Volume()) / Volume2 - q22_H2)) / ro3;
			//p3 = (((p_H2 / (ggg - 1) + ro_H2 * (u_H2 * u_H2 + v_H2 * v_H2) * 0.5) * Volume / Volume2 - T[now1] * (K->Potok_H2[3] / Volume2  - q3_H2)) - //
			//	0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
			u3 = (ro_H2 * u_H2 * Volume / Volume2 - T[now1] * (K->Potok_H2[1] / Volume2 + ro_H2 * v_H2 * u_H2 / y - q21_H2)) / ro3;
			v3 = (ro_H2 * v_H2 * Volume / Volume2 - T[now1] * (K->Potok_H2[2] / Volume2 + ro_H2 * v_H2 * v_H2 / y - q22_H2)) / ro3;
			p3 = (((p_H2 / (ggg - 1) + ro_H2 * (u_H2 * u_H2 + v_H2 * v_H2) * 0.5) * Volume / Volume2 - T[now1] * (K->Potok_H2[3] / Volume2 + //
				+v_H2 * (ggg * p_H2 / (ggg - 1) + ro_H2 * (u_H2 * u_H2 + v_H2 * v_H2) * 0.5) / y - q3_H2)) - //
				0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
			if (p3 <= 0)
			{
				p3 = 0.000001;
			}

			K->par[now2].ro_H2 = ro3;
			K->par[now2].p_H2 = p3;
			K->par[now2].u_H2 = u3;
			K->par[now2].v_H2 = v3;


			ro3 = ro_H3 * Volume / Volume2 - T[now1] * (K->Potok_H3[0] / Volume2 + ro_H3 * v_H3 / y - q1_H3);
			//ro3 = ro_H3 * Volume / Volume2 - T[now1] * (K->Potok_H3[0] / Volume2  - q1_H3);
			if (ro3 <= 0.0)
			{
				printf("Problemsssss  ro H3 < 0!  %lf,  %lf\n", x, y);
				ro3 = 0.000001;
			}
			//u3 = (ro_H3 * u_H3 * Volume / Volume2 - T[now1] * (K->Potok_H3[1] / Volume2  - q21_H3)) / ro3;
			//v3 = (ro_H3 * v_H3 * Volume / Volume2 - T[now1] * ((K->Potok_H3[2] - 2.0 * p_H3 * sin(0.5 * alpha_rot * pi_ / 180.0) * K->Get_Volume()) / Volume2  - q22_H3)) / ro3;
			//p3 = (((p_H3 / (ggg - 1) + ro_H3 * (u_H3 * u_H3 + v_H3 * v_H3) * 0.5) * Volume / Volume2 - T[now1] * (K->Potok_H3[3] / Volume2  - q3_H3)) - //
			//	0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
			u3 = (ro_H3 * u_H3 * Volume / Volume2 - T[now1] * (K->Potok_H3[1] / Volume2 + ro_H3 * v_H3 * u_H3 / y - q21_H3)) / ro3;
			v3 = (ro_H3 * v_H3 * Volume / Volume2 - T[now1] * (K->Potok_H3[2] / Volume2 + ro_H3 * v_H3 * v_H3 / y - q22_H3)) / ro3;
			p3 = (((p_H3 / (ggg - 1) + ro_H3 * (u_H3 * u_H3 + v_H3 * v_H3) * 0.5) * Volume / Volume2 - T[now1] * (K->Potok_H3[3] / Volume2 + //
				+v_H3 * (ggg * p_H3 / (ggg - 1) + ro_H3 * (u_H3 * u_H3 + v_H3 * v_H3) * 0.5) / y - q3_H3)) - //
				0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
			if (p3 <= 0)
			{
				p3 = 0.000001;
			}

			K->par[now2].ro_H3 = ro3;
			K->par[now2].p_H3 = p3;
			K->par[now2].u_H3 = u3;
			K->par[now2].v_H3 = v3;


			ro3 = ro_H4 * Volume / Volume2 - T[now1] * (K->Potok_H4[0] / Volume2 + ro_H4 * v_H4 / y - q1_H4);
			//ro3 = ro_H4 * Volume / Volume2 - T[now1] * (K->Potok_H4[0] / Volume2  - q1_H4);
			if (ro3 <= 0)
			{
				printf("Problemsssss  ro H4 < 0!  %lf,  %lf\n", x, y);
				cout << "Setka.cpp    " << ro_H1 << " " << ro_H2 << " " << ro_H3 << " " << ro_H4 << endl;
				cout << "Setka.cpp    " << p_H1 << " " << p_H2 << " " << p_H3 << " " << p_H4 << endl;
				cout << "Setka.cpp    " << ro << " " << p << " " << u << " " << v << endl;
				ro3 = 0.000001;
			}
			//u3 = (ro_H4 * u_H4 * Volume / Volume2 - T[now1] * (K->Potok_H4[1] / Volume2 - q21_H4)) / ro3;
			//v3 = (ro_H4 * v_H4 * Volume / Volume2 - T[now1] * ((K->Potok_H4[2] - 2.0 * p_H4 * sin(0.5 * alpha_rot * pi_ / 180.0) * K->Get_Volume()) / Volume2  - q22_H4)) / ro3;
			//p3 = (((p_H4 / (ggg - 1) + ro_H4 * (u_H4 * u_H4 + v_H4 * v_H4) * 0.5) * Volume / Volume2 - T[now1] * (K->Potok_H4[3] / Volume2 - q3_H4)) - //
			//	0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
			u3 = (ro_H4 * u_H4 * Volume / Volume2 - T[now1] * (K->Potok_H4[1] / Volume2 + ro_H4 * v_H4 * u_H4 / y - q21_H4)) / ro3;
			v3 = (ro_H4 * v_H4 * Volume / Volume2 - T[now1] * (K->Potok_H4[2] / Volume2 + ro_H4 * v_H4 * v_H4 / y - q22_H4)) / ro3;
			p3 = (((p_H4 / (ggg - 1) + ro_H4 * (u_H4 * u_H4 + v_H4 * v_H4) * 0.5) * Volume / Volume2 - T[now1] * (K->Potok_H4[3] / Volume2 + //
				+v_H4 * (ggg * p_H4 / (ggg - 1) + ro_H4 * (u_H4 * u_H4 + v_H4 * v_H4) * 0.5) / y - q3_H4)) - //
				0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
			if (p3 <= 0)
			{
				p3 = 0.0000001;
			}

			K->par[now2].ro_H4 = ro3;
			K->par[now2].p_H4 = p3;
			K->par[now2].u_H4 = u3;
			K->par[now2].v_H4 = v3;

			/*if (x > 600.0)
			{
				K->par[now2].ro_H4 = 1.0;
				K->par[now2].p_H4 = 0.5;
				K->par[now2].u_H4 = Velosity_inf;
				K->par[now2].v_H4 = 0.0;
			}*/


		}


		for (auto& i : this->All_Points)
		{
			i->x = i->x2;
			i->y = i->y2;
			i->Vx = 0.0;
			i->Vy = 0.0;
			i->count = 0;
		}
	}
}

void Setka::Go_5_komponent_2(int step)
// Идея в том, чтобы использовать расчёты в угле (3д задача, сведённая к 2д)
// именно этим программа отличается от первой версии
{
	//cout << "Setka.cpp    " << "START " << step << endl;
	int now1 = 1;
	int now2 = 0;
	double T[2];
	T[0] = T[1] = 0.00000001;
	mutex mut;

	for (int st = 0; st < step; st++)
	{
		//cout << "Setka.cpp    " << st << endl;
		if (st % 10000 == 0 && st > 0)
		{
			cout << "Setka.cpp    " << st << " " << T[now2] << endl;
		}
		now1 = (now1 + 1) % 2; // Какие параметры сейчас берём
		now2 = (now2 + 1) % 2; // Какие параметры сейчас меняем
		double time = 10000000;
		T[now2] = 100000000;

		this->Move_surface(now1, T[now1]);
		this->Move_Setka_Calculate_2(T[now1]);

		//omp_set_num_threads(4);
#pragma omp parallel for
		for (int num_cell = 0; num_cell < this->All_Cells.size(); num_cell++)
		{
			auto K = this->All_Cells[num_cell];
			K->Potok[0] = K->Potok[1] = K->Potok[2] = K->Potok[3] = K->Potok[4] = 0.0;
			K->Potok_H1[0] = K->Potok_H1[1] = K->Potok_H1[2] = K->Potok_H1[3] = 0.0;
			K->Potok_H2[0] = K->Potok_H2[1] = K->Potok_H2[2] = K->Potok_H2[3] = 0.0;
			K->Potok_H3[0] = K->Potok_H3[1] = K->Potok_H3[2] = K->Potok_H3[3] = 0.0;
			K->Potok_H4[0] = K->Potok_H4[1] = K->Potok_H4[2] = K->Potok_H4[3] = 0.0;

			double P[4];
			double PP[8];
			double PQ;
			Parametr par1 = K->par[now1];
			double n1, n2;
			double dist;
			double x, y;
			K->Get_Center(x, y);
			double vs1 = K->Get_Volume();
			double vs2 = K->Get_Volume_posle();

			//y = K->Get_Volume_rotate(360.0) / (2.0 * pi_ * vs1);
			//double y2 = K->Get_Volume_posle_rotate(360.0) / (2.0 * pi_ * vs2);

			double radius = sqrt(kv(x) + kv(y));
			if (radius < R11_)
			{
				continue;
			}
			double Volume = K->Get_Volume_rotate(360.0) / (2.0 * pi_);
			//cout << "Setka.cpp    " << "Vol = " << Volume << endl;
			double Volume2 = K->Get_Volume_posle_rotate(360.0) / (2.0 * pi_);
			//cout << "Setka.cpp    " << "Vol2 = " << Volume2 << endl;
			double Scw = K->Get_Volume();
			double W = 0.0;
			bool np = true;

			for (auto& i : K->Grans)
			{
				double x2, y2, x4, y4;
				double S = i->Get_square() * (i->A->y + i->B->y) * 0.5;
				//double S = i->Get_square();
				//cout << "Setka.cpp    " << "S = " << S << endl;
				i->Get_Center(x2, y2);
				i->Get_Center_posle(x4, y4);
				dist = sqrt(kv(x - x2) + kv(y - y2));
				Parametr par2, par11;
				par11 = par1;
				i->Get_par(par2, now1);
				i->Get_normal(n1, n2);
				double Vc, Vl, Vp;
				W = ((x4 - x2) * n1 + (y4 - y2) * n2) / T[now1];
				int met = 1;  // 1
				double dis = sqrt(kv(x2) + kv(y2));

				if (i->type == Usualy)
				{
					double x3, y3;
					i->Sosed->Get_Center(x3, y3);
					double dd = sqrt(kv(x3) + kv(y3));

					if (K->type == C_1)
					{
						par11.ro = par11.ro * kv(radius) / kv(dis);
						par11.ro_H1 = par11.ro_H1 * pow(radius / dis, H_pow);
						par11.Q = par11.Q * kv(radius) / kv(dis);
						par11.p = par11.p * pow(radius / dis, 2 * ggg);
						polar_perenos(x, y, x2, y2, par11.u, par11.v);
						polar_perenos(x, y, x2, y2, par11.u_H1, par11.v_H1);

						Parametr pp1;
						i->Get_par_TVD_radial(pp1, now1);
						par11.ro = pp1.ro;
						par11.p = pp1.p;
						//par11.u_H1 = pp1.u_H1;
						//par11.v_H1 = pp1.v_H1;
						par11.p_H1 = pp1.p_H1;
						par11.ro_H1 = pp1.ro_H1;
					}
					else
					{
						i->Get_par_TVD(par11, now1);
					}

					if (i->Sosed->type == C_1)
					{
						par2.ro = par2.ro * kv(dd) / kv(dis);
						par2.ro_H1 = par2.ro_H1 * pow(dd / dis, H_pow);
						par2.Q = par2.Q * kv(dd) / kv(dis);
						par2.p = par2.p * pow(dd / dis, 2 * ggg);
						polar_perenos(x3, y3, x2, y2, par2.u, par2.v);
						polar_perenos(x3, y3, x2, y2, par2.u_H1, par2.v_H1);

						Parametr pp2;
						i->Gran_copy->Get_par_TVD_radial(pp2, now1);
						par2.ro = pp2.ro;
						par2.p = pp2.p;
						//par2.u_H1 = pp2.u_H1;
						//par2.v_H1 = pp2.v_H1;
						par2.p_H1 = pp2.p_H1;
						par2.ro_H1 = pp2.ro_H1;
					}
					else
					{
						i->Gran_copy->Get_par_TVD(par2, now1);
					}
				}
				else if (i->type == Inner_sphere)
				{
					par11.ro = par11.ro * kv(radius) / kv(dis);
					par11.ro_H1 = par11.ro_H1 * pow(radius / dis, H_pow);
					par11.Q = par11.Q * kv(radius) / kv(dis);
					par11.p = par11.p * pow(radius / dis, 2 * ggg);
				}
				else if (K->type == C_1 && i->type == Axis)
				{
					par11.ro = par11.ro * kv(radius) / kv(dis);
					par11.ro_H1 = par11.ro_H1 * pow(radius / dis, H_pow);
					par11.Q = par11.Q * kv(radius) / kv(dis);
					par11.p = par11.p * pow(radius / dis, 2 * ggg);
					polar_perenos(x, y, x2, y2, par11.u, par11.v);
					polar_perenos(x, y, x2, y2, par11.u_H1, par11.v_H1);
					//par11.v = 0.0;
					par2 = par11;
					par2.v = -par2.v;
					par2.v_H1 = -par2.v_H1;
					par2.v_H2 = -par2.v_H2;
					par2.v_H3 = -par2.v_H3;
					par2.v_H4 = -par2.v_H4;
				}
				else if (K->type != C_1 && i->type == Axis)
				{
					i->Get_par_TVD(par11, now1);
					par2 = par11;
					par2.v = -par2.v;
					par2.v_H1 = -par2.v_H1;
					par2.v_H2 = -par2.v_H2;
					par2.v_H3 = -par2.v_H3;
					par2.v_H4 = -par2.v_H4;
				}


				/*if ((i->A->type == P_Contact && i->B->type != P_Contact)||(i->A->type != P_Contact && i->B->type == P_Contact))
				{
					met = 0;
				}*/

				np = false;
				bool god = false;

				if ((i->A->type == P_Contact && i->B->type == P_Contact))
				{
					god = true;
				}
				else if ((i->A->type == P_Inner_shock && i->B->type == P_Inner_shock))
				{
					god = true;
				}
				else if ((i->A->type == P_Outer_shock && i->B->type == P_Outer_shock))
				{
					god = true;
				}


				if ((i->A->type == P_Contact && i->B->type != P_Contact))
				{
					np = true;
				}
				else if ((i->A->type != P_Contact && i->B->type == P_Contact))
				{
					np = true;
				}
				else if ((i->A->type == P_Inner_shock && i->B->type != P_Inner_shock))
				{
					np = true;
				}
				else if ((i->A->type != P_Inner_shock && i->B->type == P_Inner_shock))
				{
					np = true;
				}
				else if ((i->A->type == P_Outer_shock && i->B->type != P_Outer_shock))
				{
					np = true;
				}
				else if ((i->A->type != P_Outer_shock && i->B->type == P_Outer_shock))
				{
					np = true;
				}

				//np = false;

				/*if (y < 30 && x > 120 && x < 160)
				{
					np = true;
				}*/

				if (K->type != C_centr)
				{
					time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro_H1, par11.Q, par11.p_H1, par11.u_H1, par11.v_H1, par2.ro_H1, par2.Q, //
						par2.p_H1, par2.u_H1, par2.v_H1, W, P, PQ, n1, n2, dist, met, Vl, Vc, Vp));
					for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
					{
						K->Potok_H1[k] = K->Potok_H1[k] + P[k] * S;
					}
				}

				time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro_H2, par11.Q, par11.p_H2, par11.u_H2, par11.v_H2, par2.ro_H2, par2.Q, //
					par2.p_H2, par2.u_H2, par2.v_H2, W, P, PQ, n1, n2, dist, met, Vl, Vc, Vp));
				for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
				{
					K->Potok_H2[k] = K->Potok_H2[k] + P[k] * S;
				}

				time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro_H3, par11.Q, par11.p_H3, par11.u_H3, par11.v_H3, par2.ro_H3, par2.Q, //
					par2.p_H3, par2.u_H3, par2.v_H3, W, P, PQ, n1, n2, dist, met, Vl, Vc, Vp));
				for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
				{
					K->Potok_H3[k] = K->Potok_H3[k] + P[k] * S;
				}

				time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro_H4, par11.Q, par11.p_H4, par11.u_H4, par11.v_H4, par2.ro_H4, par2.Q, //
					par2.p_H4, par2.u_H4, par2.v_H4, W, P, PQ, n1, n2, dist, met, Vl, Vc, Vp));
				for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
				{
					K->Potok_H4[k] = K->Potok_H4[k] + P[k] * S;
				}

				if (god == true)
				{
					double Vl, Vp, VV;
					vector<double> qqq1(5);
					vector<double> qqq2(5);
					vector<double> qqq(5);
					vector<double> n(3);
					n[0] = n1;
					n[1] = n2;
					n[2] = 0.0;

					qqq1[0] = par11.ro;
					qqq1[1] = par11.p;
					qqq1[2] = par11.u;
					qqq1[3] = par11.v;
					qqq1[4] = 0.0;

					qqq2[0] = par2.ro;
					qqq2[1] = par2.p;
					qqq2[2] = par2.u;
					qqq2[3] = par2.v;
					qqq2[4] = 0.0;

					auto SS = Solvers();
					time = min(time, SS.Godunov_Solver_Alexashov(qqq1, qqq2, n, qqq, Vl, Vp, VV, W, dist));

					K->Potok[0] = K->Potok[0] + qqq[0] * S;
					K->Potok[1] = K->Potok[1] + qqq[1] * S;
					K->Potok[2] = K->Potok[2] + qqq[2] * S;
					K->Potok[3] = K->Potok[3] + qqq[4] * S;

					time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro, par11.Q, par11.p, par11.u, par11.v, par2.ro, par2.Q, //
						par2.p, par2.u, par2.v, W, P, PQ, n1, n2, dist, met, Vl, Vc, Vp, np));

					K->Potok[4] = K->Potok[4] + PQ * S;
				}
				else
				{
					time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro, par11.Q, par11.p, par11.u, par11.v, par2.ro, par2.Q, //
						par2.p, par2.u, par2.v, W, P, PQ, n1, n2, dist, met, Vl, Vc, Vp, np));
					for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
					{
						K->Potok[k] = K->Potok[k] + P[k] * S;
					}
					K->Potok[4] = K->Potok[4] + PQ * S;
				}

			}


			// Блок учёта боковых граней (которые искусственные из-за трёхмерности)

			if (false)
			{
				//cout << "Setka.cpp    " << K->Potok_H1[0] << " " << K->Potok_H1[1] << " " << K->Potok_H1[2] - 2.0 * par1.p_H1 * sin(0.5 * alpha_rot * pi_ / 180.0) * Scw << " " << K->Potok_H1[3] << endl;
				Solvers SS = Solvers();
				//-0.0436194 0.999048
				double uu = par1.v_H1, vv = 0.0;
				polar_provorot(pi_ * alpha_rot / 360.0, uu, vv);
				SS.HLLC_Aleksashov(par1.ro_H1, par1.p_H1, par1.u_H1, uu, vv,//
					par1.ro_H1, par1.p_H1, par1.u_H1, uu, vv,//
					P, 0.0, -0.0436194, 0.999048, 1.0);
				//K->Potok_H1[0] = K->Potok_H1[0] + P[0] * Scw;
				//K->Potok_H1[1] = K->Potok_H1[1] + P[1] * Scw;
				K->Potok_H1[2] = K->Potok_H1[2] + P[2] * Scw;
				//K->Potok_H1[3] = K->Potok_H1[3] + P[4] * Scw;

				uu = par1.v_H1, vv = 0.0;
				polar_provorot(-pi_ * alpha_rot / 360.0, uu, vv);
				SS.HLLC_Aleksashov(par1.ro_H1, par1.p_H1, par1.u_H1, uu, vv,//
					par1.ro_H1, par1.p_H1, par1.u_H1, uu, vv,//
					P, 0.0, -0.0436194, -0.999048, 1.0);
				//K->Potok_H1[0] = K->Potok_H1[0] + P[0] * Scw;
				//K->Potok_H1[1] = K->Potok_H1[1] + P[1] * Scw;
				K->Potok_H1[2] = K->Potok_H1[2] + P[2] * Scw;
				//K->Potok_H1[3] = K->Potok_H1[3] + P[4] * Scw;

				//cout << "Setka.cpp    " << K->Potok_H1[0] << " " << K->Potok_H1[1] << " " << K->Potok_H1[2] << " " << K->Potok_H1[3] << endl;
				//cout << "Setka.cpp    " << K->Potok_H2[0] << " " << K->Potok_H2[1] << " " << K->Potok_H2[2] - 2.0 * par1.p_H2 * sin(0.5 * alpha_rot * pi_ / 180.0) * Scw << " " << K->Potok_H2[3] << endl;
				uu = par1.v_H2, vv = 0.0;
				polar_provorot(pi_ * alpha_rot / 360.0, uu, vv);
				SS.HLLC_Aleksashov(par1.ro_H2, par1.p_H2, par1.u_H2, uu, vv,//
					par1.ro_H2, par1.p_H2, par1.u_H2, uu, vv,//
					P, 0.0, -0.0436194, 0.999048, 1.0);
				//K->Potok_H2[0] = K->Potok_H2[0] + P[0] * Scw;
				//K->Potok_H2[1] = K->Potok_H2[1] + P[1] * Scw;
				K->Potok_H2[2] = K->Potok_H2[2] + P[2] * Scw;
				//K->Potok_H2[3] = K->Potok_H2[3] + P[4] * Scw;

				uu = par1.v_H2, vv = 0.0;
				polar_provorot(-pi_ * alpha_rot / 360.0, uu, vv);
				SS.HLLC_Aleksashov(par1.ro_H2, par1.p_H2, par1.u_H2, uu, vv,//
					par1.ro_H2, par1.p_H2, par1.u_H2, uu, vv,//
					P, 0.0, -0.0436194, -0.999048, 1.0);
				//K->Potok_H2[0] = K->Potok_H2[0] + P[0] * Scw;
				//K->Potok_H2[1] = K->Potok_H2[1] + P[1] * Scw;
				K->Potok_H2[2] = K->Potok_H2[2] + P[2] * Scw;
				//K->Potok_H2[3] = K->Potok_H2[3] + P[4] * Scw;
				//cout << "Setka.cpp    " << K->Potok_H2[0] << " " << K->Potok_H2[1] << " " << K->Potok_H2[2] << " " << K->Potok_H2[3] << endl;
				//cout << "Setka.cpp    " << K->Potok_H3[0] << " " << K->Potok_H3[1] << " " << K->Potok_H3[2] - 2.0 * par1.p_H3 * sin(0.5 * alpha_rot * pi_ / 180.0) * Scw << " " << K->Potok_H3[3] << endl;
				uu = par1.v_H3, vv = 0.0;
				polar_provorot(pi_ * alpha_rot / 360.0, uu, vv);
				SS.HLLC_Aleksashov(par1.ro_H3, par1.p_H3, par1.u_H3, uu, vv,//
					par1.ro_H3, par1.p_H3, par1.u_H3, uu, vv,//
					P, 0.0, -0.0436194, 0.999048, 1.0);
				//K->Potok_H3[0] = K->Potok_H3[0] + P[0] * Scw;
				//K->Potok_H3[1] = K->Potok_H3[1] + P[1] * Scw;
				K->Potok_H3[2] = K->Potok_H3[2] + P[2] * Scw;
				//K->Potok_H3[3] = K->Potok_H3[3] + P[4] * Scw;

				uu = par1.v_H3, vv = 0.0;
				polar_provorot(-pi_ * alpha_rot / 360.0, uu, vv);
				SS.HLLC_Aleksashov(par1.ro_H3, par1.p_H3, par1.u_H3, uu, vv,//
					par1.ro_H3, par1.p_H3, par1.u_H3, uu, vv,//
					P, 0.0, -0.0436194, -0.999048, 1.0);
				//K->Potok_H3[0] = K->Potok_H3[0] + P[0] * Scw;
				//K->Potok_H3[1] = K->Potok_H3[1] + P[1] * Scw;
				K->Potok_H3[2] = K->Potok_H3[2] + P[2] * Scw;
				//K->Potok_H3[3] = K->Potok_H3[3] + P[4] * Scw;
				//cout << "Setka.cpp    " << K->Potok_H3[0] << " " << K->Potok_H3[1] << " " << K->Potok_H3[2] << " " << K->Potok_H3[3] << endl;
				//cout << "Setka.cpp    " << K->Potok_H4[0] << " " << K->Potok_H4[1] << " " << K->Potok_H4[2] - 2.0 * par1.p_H4 * sin(0.5 * alpha_rot * pi_ / 180.0) * Scw << " " << K->Potok_H4[3] << endl;
				uu = par1.v_H4, vv = 0.0;
				polar_provorot(pi_ * alpha_rot / 360.0, uu, vv);
				SS.HLLC_Aleksashov(par1.ro_H4, par1.p_H4, par1.u_H4, uu, vv,//
					par1.ro_H4, par1.p_H4, par1.u_H4, uu, vv,//
					P, 0.0, -0.0436194, 0.999048, 1.0);
				//K->Potok_H4[0] = K->Potok_H4[0] + P[0] * Scw;
				//K->Potok_H4[1] = K->Potok_H4[1] + P[1] * Scw;
				K->Potok_H4[2] = K->Potok_H4[2] + P[2] * Scw;
				//K->Potok_H4[3] = K->Potok_H4[3] + P[4] * Scw;

				uu = par1.v_H4, vv = 0.0;
				polar_provorot(-pi_ * alpha_rot / 360.0, uu, vv);
				SS.HLLC_Aleksashov(par1.ro_H4, par1.p_H4, par1.u_H4, uu, vv,//
					par1.ro_H4, par1.p_H4, par1.u_H4, uu, vv,//
					P, 0.0, -0.0436194, -0.999048, 1.0);
				//K->Potok_H4[0] = K->Potok_H4[0] + P[0] * Scw;
				//K->Potok_H4[1] = K->Potok_H4[1] + P[1] * Scw;
				K->Potok_H4[2] = K->Potok_H4[2] + P[2] * Scw;
				//K->Potok_H4[3] = K->Potok_H4[3] + P[4] * Scw;

				//cout << "Setka.cpp    " << K->Potok_H4[0] << " " << K->Potok_H4[1] << " " << K->Potok_H4[2] << " " << K->Potok_H4[3] << endl;
				//exit(-1);

				//cout << "Setka.cpp    " << K->Potok[0] << " " << K->Potok[1] << " " << K->Potok[2] - 2.0 * par1.p * sin(0.5 * alpha_rot * pi_ / 180.0) * Scw << " " << K->Potok[3] << " " << K->Potok[4] << endl;
				uu = par1.v, vv = 0.0;
				polar_provorot(pi_ * alpha_rot / 360.0, uu, vv);
				SS.HLLC_Aleksashov(par1.ro, par1.p, par1.u, uu, vv, //
					par1.ro, par1.p, par1.u, uu, vv,//
					P, 0.0, -0.0436194, 0.999048, 1.0);
				//K->Potok[0] = K->Potok[0] + PP[0] * Scw;
				//K->Potok[1] = K->Potok[1] + PP[1] * Scw;
				K->Potok[2] = K->Potok[2] + P[2] * Scw;
				//K->Potok[3] = K->Potok[3] + PP[4] * Scw;
				//K->Potok[4] = K->Potok[4] + PQ * Scw;

				uu = par1.v, vv = 0.0;
				polar_provorot(-pi_ * alpha_rot / 360.0, uu, vv);
				SS.HLLC_Aleksashov(par1.ro, par1.p, par1.u, uu, vv, //
					par1.ro, par1.p, par1.u, uu, vv,//
					P, 0.0, -0.0436194, -0.999048, 1.0);
				//K->Potok[0] = K->Potok[0] + PP[0] * Scw;
				//K->Potok[1] = K->Potok[1] + PP[1] * Scw;
				K->Potok[2] = K->Potok[2] + P[2] * Scw;
				//K->Potok[3] = K->Potok[3] + PP[4] * Scw;
				//K->Potok[4] = K->Potok[4] + PQ * Scw;

				//cout << "Setka.cpp    " << K->Potok[0] << " " << K->Potok[1] << " " << K->Potok[2] << " " << K->Potok[3] << " " << K->Potok[4] << endl;
				//cout << "Setka.cpp    " << "A  " << endl;
				//cout << "Setka.cpp    " << -2.0 * par1.p * sin(0.5 * alpha_rot * pi_ / 180.0) * Scw << endl;
				//exit(-1);
			}


			// Блок законов сохранения в ячейке

			if (time < T[now2])
			{
				mut.lock();
				T[now2] = time;
				mut.unlock();
			}



			double U_M_H1, U_M_H2, U_M_H3, U_M_H4;
			double U_H1, U_H2, U_H3, U_H4;
			double sigma_H1, sigma_H2, sigma_H3, sigma_H4;
			double nu_H1, nu_H2, nu_H3, nu_H4;
			double q2_1, q2_2, q3;
			double u, v, ro, p, Q;
			if (false)//(K->type == C_centr || K->type == C_1 || K->type == C_2 || K->type == C_3)//(par1.Q / par1.ro < 90)//
			{
				u = par1.u * (chi_real / chi_);
				v = par1.v * (chi_real / chi_);
				ro = par1.ro / kv(chi_real / chi_);
				p = par1.p;
				Q = par1.Q / kv(chi_real / chi_);
			}
			else
			{
				u = par1.u;
				v = par1.v;
				ro = par1.ro;
				p = par1.p;
				Q = par1.Q;
			}


			double u_H1 = par1.u_H1, v_H1 = par1.v_H1, ro_H1 = par1.ro_H1, p_H1 = par1.p_H1;
			double u_H2 = par1.u_H2, v_H2 = par1.v_H2, ro_H2 = par1.ro_H2, p_H2 = par1.p_H2;
			double u_H3 = par1.u_H3, v_H3 = par1.v_H3, ro_H3 = par1.ro_H3, p_H3 = par1.p_H3;
			double u_H4 = par1.u_H4, v_H4 = par1.v_H4, ro_H4 = par1.ro_H4, p_H4 = par1.p_H4;

			U_M_H1 = sqrt(kv(u - u_H1) + kv(v - v_H1) + (64.0 / (9.0 * pi)) //
				* (p / ro + 2.0 * p_H1 / ro_H1));
			U_M_H2 = sqrt(kv(u - u_H2) + kv(v - v_H2) + (64.0 / (9.0 * pi)) //
				* (p / ro + 2.0 * p_H2 / ro_H2));
			U_M_H3 = sqrt(kv(u - u_H3) + kv(v - v_H3) + (64.0 / (9.0 * pi)) //
				* (p / ro + 2.0 * p_H3 / ro_H3));
			U_M_H4 = sqrt(kv(u - u_H4) + kv(v - v_H4) + (64.0 / (9.0 * pi)) //
				* (p / ro + 2.0 * p_H4 / ro_H4));

			U_H1 = sqrt(kv(u - u_H1) + kv(v - v_H1) + (4.0 / pi) //
				* (p / ro + 2.0 * p_H1 / ro_H1));
			U_H2 = sqrt(kv(u - u_H2) + kv(v - v_H2) + (4.0 / pi) //
				* (p / ro + 2.0 * p_H2 / ro_H2));
			U_H3 = sqrt(kv(u - u_H3) + kv(v - v_H3) + (4.0 / pi) //
				* (p / ro + 2.0 * p_H3 / ro_H3));
			U_H4 = sqrt(kv(u - u_H4) + kv(v - v_H4) + (4.0 / pi) //
				* (p / ro + 2.0 * p_H4 / ro_H4));

			sigma_H1 = kv(1.0 - a_2 * log(U_M_H1)); // 0.1243
			sigma_H2 = kv(1.0 - a_2 * log(U_M_H2));
			sigma_H3 = kv(1.0 - a_2 * log(U_M_H3)); // 0.1121     a_2
			sigma_H4 = kv(1.0 - a_2 * log(U_M_H4));

			nu_H1 = ro * ro_H1 * U_M_H1 * sigma_H1;
			nu_H2 = ro * ro_H2 * U_M_H2 * sigma_H2;
			nu_H3 = ro * ro_H3 * U_M_H3 * sigma_H3;
			nu_H4 = ro * ro_H4 * U_M_H4 * sigma_H4;

			q2_1 = (n_p_LISM_ / Kn_) * (nu_H1 * (u_H1 - u) + nu_H2 * (u_H2 - u) //
				+ nu_H3 * (u_H3 - u) + nu_H4 * (u_H4 - u));
			q2_2 = (n_p_LISM_ / Kn_) * (nu_H1 * (v_H1 - v) + nu_H2 * (v_H2 - v) //
				+ nu_H3 * (v_H3 - v) + nu_H4 * (v_H4 - v));

			if (false)//(K->type == C_centr || K->type == C_1 || K->type == C_2 || K->type == C_3)//(par1.Q / par1.ro < 90)//
			{
				q3 = (n_p_LISM_ / Kn_) * (nu_H1 * ((kv(u_H1) + kv(v_H1) - kv(u) - kv(v)) / 2.0 + //
					(U_H1 / U_M_H1) * (2.0 * p_H1 / ro_H1 - p / ro)) + //
					nu_H2 * ((kv(u_H2) + kv(v_H2) - kv(u) - kv(v)) / 2.0 + //
						(U_H2 / U_M_H2) * (2.0 * p_H2 / ro_H2 - p / ro)) + //
					nu_H3 * ((kv(u_H3) + kv(v_H3) - kv(u) - kv(v)) / 2.0 + //
						(U_H3 / U_M_H3) * (2.0 * p_H3 / ro_H3 - p / ro)) + //
					nu_H4 * ((kv(u_H4) + kv(v_H4) - kv(u) - kv(v)) / 2.0 + //
						(U_H4 / U_M_H4) * (2.0 * p_H4 / ro_H4 - p / ro))) / (chi_real / chi_);
			}
			else
			{
				q3 = (n_p_LISM_ / Kn_) * (nu_H1 * ((kv(u_H1) + kv(v_H1) - kv(u) - kv(v)) / 2.0 + //
					(U_H1 / U_M_H1) * (2.0 * p_H1 / ro_H1 - p / ro)) + //
					nu_H2 * ((kv(u_H2) + kv(v_H2) - kv(u) - kv(v)) / 2.0 + //
						(U_H2 / U_M_H2) * (2.0 * p_H2 / ro_H2 - p / ro)) + //
					nu_H3 * ((kv(u_H3) + kv(v_H3) - kv(u) - kv(v)) / 2.0 + //
						(U_H3 / U_M_H3) * (2.0 * p_H3 / ro_H3 - p / ro)) + //
					nu_H4 * ((kv(u_H4) + kv(v_H4) - kv(u) - kv(v)) / 2.0 + //
						(U_H4 / U_M_H4) * (2.0 * p_H4 / ro_H4 - p / ro)));
			}


			/*q2_1 = 0.0;
			q2_2 = 0.0;
			q3 = 0.0;*/

			double ro2_p, p2_p, V1_p, V2_p, QQ2;
			double ro3, p3, u3, v3, Q33;
			double kappa = Q / ro;


			u = par1.u;
			v = par1.v;
			ro = par1.ro;
			p = par1.p;
			Q = par1.Q;

			if (K->type != C_centr)
			{
				ro3 = ro * Volume / Volume2 - T[now1] * (K->Potok[0] / Volume2);  // В декартовых
				//ro3 = ro * Volume / Volume2 - T[now1] * (K->Potok[0] / Volume2 + ro * v / y);  // В цилиндрических

				Q33 = Q * Volume / Volume2 - (T[now1] / Volume2) * K->Potok[4];
				//Q33 = Q * Volume / Volume2 - (T[now1] / Volume2) * K->Potok[4] - T[now1] * Q * v / y;
				if (ro3 <= 0.0)
				{
					printf("Problemsssss  ro < 0!    %lf, %lf\n", x, y);
					ro3 = 0.00001;
				}

				u3 = (ro * u * Volume / Volume2 - T[now1] * (K->Potok[1] / Volume2  - q2_1)) / ro3;
				v3 = (ro * v * Volume / Volume2 - T[now1] * ((K->Potok[2] - p * 0.5 * (vs1 + vs2)) / Volume2  - q2_2)) / ro3;
				//v3 = (ro * v * Volume / Volume2 - T[now1] * (K->Potok[2] / Volume2 - q2_2)) / ro3;
				p3 = (((p / (ggg - 1) + ro * (u * u + v * v) * 0.5) * Volume / Volume2 - T[now1] * (K->Potok[3] / Volume2 - q3)) - //
					0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
				
				//u3 = (ro * u * Volume / Volume2 - T[now1] * (K->Potok[1] / Volume2 + ro * v * u / y - q2_1)) / ro3;
				//v3 = (ro * v * Volume / Volume2 - T[now1] * (K->Potok[2] / Volume2 + ro * v * v / y - q2_2)) / ro3;
				//p3 = (((p / (ggg - 1) + ro * (u * u + v * v) * 0.5) * Volume / Volume2 - T[now1] * (K->Potok[3] / Volume2 + //
				//	+v * (ggg * p / (ggg - 1) + ro * (u * u + v * v) * 0.5) / y - q3)) - //
				//	0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);

				if (p3 <= 0)
				{
					p3 = 0.000001;
				}

				K->par[now2].ro = ro3;
				K->par[now2].Q = Q33;
				K->par[now2].p = p3;
				K->par[now2].u = u3;
				K->par[now2].v = v3;
			}
			else
			{
				K->par[now2].ro = K->par[now1].ro;
				K->par[now2].Q = K->par[now1].Q;
				K->par[now2].p = K->par[now1].p;
				K->par[now2].u = K->par[now1].u;
				K->par[now2].v = K->par[now1].v;
			}

			if (false)//(K->type == C_centr || K->type == C_1 || K->type == C_2 || K->type == C_3)//(par1.Q / par1.ro < 90)//
			{
				u = par1.u * (chi_real / chi_);
				v = par1.v * (chi_real / chi_);
				ro = par1.ro / kv(chi_real / chi_);
				p = par1.p;
				Q = par1.Q / kv(chi_real / chi_);
			}
			else
			{
				u = par1.u;
				v = par1.v;
				ro = par1.ro;
				p = par1.p;
				Q = par1.Q;
			}


			int k1 = 0, k2 = 0, k3 = 0, k4 = 0;

			if (K->type == C_centr || K->type == C_1 || K->type == C_2 || K->type == C_3)//(par1.Q / par1.ro < 90)//
			{
				k3 = 0;
				k4 = 0;
				if (K->type == C_1 || K->type == C_centr)//(((kv(u) + kv(v)) / (ggg * p / ro) > 1.0) && ((radius <= 250.0)))
				{
					k1 = 1;
					k2 = 0;
				}
				else
				{
					k1 = 0;
					k2 = 1;
				}
			}
			else
			{
				k1 = 0;
				k2 = 0;
				if (K->type == C_4)//(p / ro > 1.8)   // ( ((y < 8.0)&&(p / ro > (y * (-0.0238) + 0.36)))||( (y >= 8.0)&&(p / ro > 0.17) ) )
				{
					k4 = 0;
					k3 = 1;
				}
				else
				{
					k3 = 0;
					k4 = 1;
				}
			}

			double S1, S2;
			double q1_H1, q1_H2, q1_H3, q1_H4;
			double q21_H1, q21_H2, q21_H3, q21_H4;
			double q22_H1, q22_H2, q22_H3, q22_H4;
			double q3_H1, q3_H2, q3_H3, q3_H4;

			S1 = nu_H1 + nu_H2 + nu_H3 + nu_H4;
			S2 = nu_H1 * ((kv(u) + kv(v)) / 2.0 + (U_H1 / U_M_H1) * (p / ro))//
				+ nu_H2 * ((kv(u) + kv(v)) / 2.0 + (U_H2 / U_M_H2) * (p / ro))//
				+ nu_H3 * ((kv(u) + kv(v)) / 2.0 + (U_H3 / U_M_H3) * (p / ro))//
				+ nu_H4 * ((kv(u) + kv(v)) / 2.0 + (U_H4 / U_M_H4) * (p / ro));

			q1_H1 = (n_H_LISM_ / Kn_) * (k1 * S1 - nu_H1);
			q1_H2 = (n_H_LISM_ / Kn_) * (k2 * S1 - nu_H2);
			q1_H3 = (n_H_LISM_ / Kn_) * (k3 * S1 - nu_H3);
			q1_H4 = (n_H_LISM_ / Kn_) * (k4 * S1 - nu_H4);

			q21_H1 = (n_H_LISM_ / Kn_) * (k1 * S1 * u - nu_H1 * u_H1);
			q21_H2 = (n_H_LISM_ / Kn_) * (k2 * S1 * u - nu_H2 * u_H2);
			q21_H3 = (n_H_LISM_ / Kn_) * (k3 * S1 * u - nu_H3 * u_H3);
			q21_H4 = (n_H_LISM_ / Kn_) * (k4 * S1 * u - nu_H4 * u_H4);

			q22_H1 = (n_H_LISM_ / Kn_) * (k1 * S1 * v - nu_H1 * v_H1);
			q22_H2 = (n_H_LISM_ / Kn_) * (k2 * S1 * v - nu_H2 * v_H2);
			q22_H3 = (n_H_LISM_ / Kn_) * (k3 * S1 * v - nu_H3 * v_H3);
			q22_H4 = (n_H_LISM_ / Kn_) * (k4 * S1 * v - nu_H4 * v_H4);

			q3_H1 = (n_H_LISM_ / Kn_) * (k1 * S2 - nu_H1 * ((kv(u_H1) + kv(v_H1)) / 2.0 + (U_H1 / U_M_H1) * 2.0 * (p_H1 / ro_H1)));
			q3_H2 = (n_H_LISM_ / Kn_) * (k2 * S2 - nu_H2 * ((kv(u_H2) + kv(v_H2)) / 2.0 + (U_H2 / U_M_H2) * 2.0 * (p_H2 / ro_H2)));
			q3_H3 = (n_H_LISM_ / Kn_) * (k3 * S2 - nu_H3 * ((kv(u_H3) + kv(v_H3)) / 2.0 + (U_H3 / U_M_H3) * 2.0 * (p_H3 / ro_H3)));
			q3_H4 = (n_H_LISM_ / Kn_) * (k4 * S2 - nu_H4 * ((kv(u_H4) + kv(v_H4)) / 2.0 + (U_H4 / U_M_H4) * 2.0 * (p_H4 / ro_H4)));


			if (K->type != C_centr)
			{
				//ro3 = ro_H1 * Volume / Volume2 - T[now1] * (K->Potok_H1[0] / Volume2 + ro_H1 * v_H1 / y - q1_H1);
				ro3 = ro_H1 * Volume / Volume2 - T[now1] * (K->Potok_H1[0] / Volume2 - q1_H1);
				if (ro3 <= 0)
				{
					printf("Problemsssss  ro H1 < 0!   %lf,  %lf\n", x, y);
					ro3 = 0.0000001;
				}
				u3 = (ro_H1 * u_H1 * Volume / Volume2 - T[now1] * (K->Potok_H1[1] / Volume2  - q21_H1)) / ro3;
				v3 = (ro_H1 * v_H1 * Volume / Volume2 - T[now1] * ((K->Potok_H1[2] - p_H1 * 0.5 * (vs1 + vs2)) / Volume2  - q22_H1)) / ro3;
				//v3 = (ro_H1 * v_H1 * Volume / Volume2 - T[now1] * (K->Potok_H1[2] / Volume2 - q22_H1)) / ro3;
				p3 = (((p_H1 / (ggg - 1) + ro_H1 * (u_H1 * u_H1 + v_H1 * v_H1) * 0.5) * Volume / Volume2 - T[now1] * (K->Potok_H1[3] / Volume2  - q3_H1)) - //
					0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);

				//u3 = (ro_H1 * u_H1 * Volume / Volume2 - T[now1] * (K->Potok_H1[1] / Volume2 + ro_H1 * v_H1 * u_H1 / y - q21_H1)) / ro3;
				//v3 = (ro_H1 * v_H1 * Volume / Volume2 - T[now1] * (K->Potok_H1[2] / Volume2 + ro_H1 * v_H1 * v_H1 / y - q22_H1)) / ro3;
				//p3 = (((p_H1 / (ggg - 1) + ro_H1 * (u_H1 * u_H1 + v_H1 * v_H1) * 0.5) * Volume / Volume2 - T[now1] * (K->Potok_H1[3] / Volume2 + //
				//	+v_H1 * (ggg * p_H1 / (ggg - 1) + ro_H1 * (u_H1 * u_H1 + v_H1 * v_H1) * 0.5) / y - q3_H1)) - //
				//	0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
				if (p3 <= 0)
				{
					p3 = 0.000001;
				}

				K->par[now2].ro_H1 = ro3;
				K->par[now2].p_H1 = p3;
				K->par[now2].u_H1 = u3;
				K->par[now2].v_H1 = v3;

				/*if (x > 600.0)
				{
					K->par[now2].ro_H1 = 1.0;
					K->par[now2].p_H1 = 1.0;
					K->par[now2].u_H1 = Velosity_inf;
					K->par[now2].v_H1 = 0.0;
				}*/
			}
			else
			{
				K->par[now2].ro_H1 = K->par[now1].ro_H1;
				K->par[now2].p_H1 = K->par[now1].p_H1;
				K->par[now2].u_H1 = K->par[now1].u_H1;
				K->par[now2].v_H1 = K->par[now1].v_H1;
			}


			//ro3 = ro_H2 * Volume / Volume2 - T[now1] * (K->Potok_H2[0] / Volume2 + ro_H2 * v_H2 / y - q1_H2);
			ro3 = ro_H2 * Volume / Volume2 - T[now1] * (K->Potok_H2[0] / Volume2 - q1_H2);
			if (ro3 <= 0)
			{
				printf("Problemsssss  ro H2 < 0!  %lf,  %lf\n", x, y);
				ro3 = 0.000001;
			}
			u3 = (ro_H2 * u_H2 * Volume / Volume2 - T[now1] * (K->Potok_H2[1] / Volume2 - q21_H2)) / ro3;
			v3 = (ro_H2 * v_H2 * Volume / Volume2 - T[now1] * ((K->Potok_H2[2] - p_H2 * 0.5 * (vs1 + vs2)) / Volume2 - q22_H2)) / ro3;
			//v3 = (ro_H2 * v_H2 * Volume / Volume2 - T[now1] * (K->Potok_H2[2] / Volume2 - q22_H2)) / ro3;
			p3 = (((p_H2 / (ggg - 1) + ro_H2 * (u_H2 * u_H2 + v_H2 * v_H2) * 0.5) * Volume / Volume2 - T[now1] * (K->Potok_H2[3] / Volume2  - q3_H2)) - //
				0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
			
			//u3 = (ro_H2 * u_H2 * Volume / Volume2 - T[now1] * (K->Potok_H2[1] / Volume2 + ro_H2 * v_H2 * u_H2 / y - q21_H2)) / ro3;
			//v3 = (ro_H2 * v_H2 * Volume / Volume2 - T[now1] * (K->Potok_H2[2] / Volume2 + ro_H2 * v_H2 * v_H2 / y - q22_H2)) / ro3;
			//p3 = (((p_H2 / (ggg - 1) + ro_H2 * (u_H2 * u_H2 + v_H2 * v_H2) * 0.5) * Volume / Volume2 - T[now1] * (K->Potok_H2[3] / Volume2 + //
			//	+v_H2 * (ggg * p_H2 / (ggg - 1) + ro_H2 * (u_H2 * u_H2 + v_H2 * v_H2) * 0.5) / y - q3_H2)) - //
			//	0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
			if (p3 <= 0)
			{
				p3 = 0.000001;
			}

			K->par[now2].ro_H2 = ro3;
			K->par[now2].p_H2 = p3;
			K->par[now2].u_H2 = u3;
			K->par[now2].v_H2 = v3;


			//ro3 = ro_H3 * Volume / Volume2 - T[now1] * (K->Potok_H3[0] / Volume2 + ro_H3 * v_H3 / y - q1_H3);
			ro3 = ro_H3 * Volume / Volume2 - T[now1] * (K->Potok_H3[0] / Volume2  - q1_H3);
			if (ro3 <= 0.0)
			{
				printf("Problemsssss  ro H3 < 0!  %lf,  %lf\n", x, y);
				ro3 = 0.000001;
			}
			u3 = (ro_H3 * u_H3 * Volume / Volume2 - T[now1] * (K->Potok_H3[1] / Volume2  - q21_H3)) / ro3;
			v3 = (ro_H3 * v_H3 * Volume / Volume2 - T[now1] * ((K->Potok_H3[2] - p_H3 * 0.5 * (vs1 + vs2)) / Volume2  - q22_H3)) / ro3;
			//v3 = (ro_H3 * v_H3 * Volume / Volume2 - T[now1] * (K->Potok_H3[2] / Volume2 - q22_H3)) / ro3;
			p3 = (((p_H3 / (ggg - 1) + ro_H3 * (u_H3 * u_H3 + v_H3 * v_H3) * 0.5) * Volume / Volume2 - T[now1] * (K->Potok_H3[3] / Volume2  - q3_H3)) - //
				0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
			
			//u3 = (ro_H3 * u_H3 * Volume / Volume2 - T[now1] * (K->Potok_H3[1] / Volume2 + ro_H3 * v_H3 * u_H3 / y - q21_H3)) / ro3;
			//v3 = (ro_H3 * v_H3 * Volume / Volume2 - T[now1] * (K->Potok_H3[2] / Volume2 + ro_H3 * v_H3 * v_H3 / y - q22_H3)) / ro3;
			//p3 = (((p_H3 / (ggg - 1) + ro_H3 * (u_H3 * u_H3 + v_H3 * v_H3) * 0.5) * Volume / Volume2 - T[now1] * (K->Potok_H3[3] / Volume2 + //
			//	+v_H3 * (ggg * p_H3 / (ggg - 1) + ro_H3 * (u_H3 * u_H3 + v_H3 * v_H3) * 0.5) / y - q3_H3)) - //
			//	0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
			if (p3 <= 0)
			{
				p3 = 0.000001;
			}

			K->par[now2].ro_H3 = ro3;
			K->par[now2].p_H3 = p3;
			K->par[now2].u_H3 = u3;
			K->par[now2].v_H3 = v3;


			//ro3 = ro_H4 * Volume / Volume2 - T[now1] * (K->Potok_H4[0] / Volume2 + ro_H4 * v_H4 / y - q1_H4);
			ro3 = ro_H4 * Volume / Volume2 - T[now1] * (K->Potok_H4[0] / Volume2  - q1_H4);
			if (ro3 <= 0)
			{
				printf("Problemsssss  ro H4 < 0!  %lf,  %lf\n", x, y);
				cout << "Setka.cpp    " << ro_H1 << " " << ro_H2 << " " << ro_H3 << " " << ro_H4 << endl;
				cout << "Setka.cpp    " << p_H1 << " " << p_H2 << " " << p_H3 << " " << p_H4 << endl;
				cout << "Setka.cpp    " << ro << " " << p << " " << u << " " << v << endl;
				ro3 = 0.000001;
			}
			u3 = (ro_H4 * u_H4 * Volume / Volume2 - T[now1] * (K->Potok_H4[1] / Volume2 - q21_H4)) / ro3;
			v3 = (ro_H4 * v_H4 * Volume / Volume2 - T[now1] * ((K->Potok_H4[2] - p_H4 * 0.5 * (vs1 + vs2)) / Volume2  - q22_H4)) / ro3;
			//v3 = (ro_H4 * v_H4 * Volume / Volume2 - T[now1] * (K->Potok_H4[2] / Volume2 - q22_H4)) / ro3;
			p3 = (((p_H4 / (ggg - 1) + ro_H4 * (u_H4 * u_H4 + v_H4 * v_H4) * 0.5) * Volume / Volume2 - T[now1] * (K->Potok_H4[3] / Volume2 - q3_H4)) - //
				0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);

			//u3 = (ro_H4 * u_H4 * Volume / Volume2 - T[now1] * (K->Potok_H4[1] / Volume2 + ro_H4 * v_H4 * u_H4 / y - q21_H4)) / ro3;
			//v3 = (ro_H4 * v_H4 * Volume / Volume2 - T[now1] * (K->Potok_H4[2] / Volume2 + ro_H4 * v_H4 * v_H4 / y - q22_H4)) / ro3;
			//p3 = (((p_H4 / (ggg - 1) + ro_H4 * (u_H4 * u_H4 + v_H4 * v_H4) * 0.5) * Volume / Volume2 - T[now1] * (K->Potok_H4[3] / Volume2 + //
			//	+v_H4 * (ggg * p_H4 / (ggg - 1) + ro_H4 * (u_H4 * u_H4 + v_H4 * v_H4) * 0.5) / y - q3_H4)) - //
			//	0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
			if (p3 <= 0)
			{
				p3 = 0.0000001;
			}

			K->par[now2].ro_H4 = ro3;
			K->par[now2].p_H4 = p3;
			K->par[now2].u_H4 = u3;
			K->par[now2].v_H4 = v3;

			/*if (x > 600.0)
			{
				K->par[now2].ro_H4 = 1.0;
				K->par[now2].p_H4 = 0.5;
				K->par[now2].u_H4 = Velosity_inf;
				K->par[now2].v_H4 = 0.0;
			}*/


		}


		for (auto& i : this->All_Points)
		{
			i->x = i->x2;
			i->y = i->y2;
			i->Vx = 0.0;
			i->Vy = 0.0;
			i->count = 0;
		}
	}
}

void Setka::Go_5_komponent_MK(int step, bool dvig)
{
	//cout << "Setka.cpp    " << "START " << step << endl;
	int now1 = 1;
	int now2 = 0;
	double T[2];
	T[0] = T[1] = 0.00000001;
	mutex mut;

	for (int st = 0; st < step; st++)
	{
		/*if (st % 1000 == 0 && st > 0)
		{
			cout << "Setka.cpp    " << st << " " << T[now2] << endl;
		}*/
		now1 = (now1 + 1) % 2; // Какие параметры сейчас берём
		now2 = (now2 + 1) % 2; // Какие параметры сейчас меняем
		double time = 10000000;
		T[now2] = 100000000;


		if (dvig == true)
		{
			this->Move_surface(now1, T[now1]);
			this->Move_Setka_Calculate_2(T[now1]);
		}

		//omp_set_num_threads(4);
#pragma omp parallel for
		for (int num_cell = 0; num_cell < this->All_Cells.size(); num_cell++)
		{
			auto K = this->All_Cells[num_cell];
			K->Potok[0] = K->Potok[1] = K->Potok[2] = K->Potok[3] = K->Potok[4] = 0.0;

			double P[4];
			double PQ;
			Parametr par1 = K->par[now1];
			double n1, n2;
			double dist;
			double x, y;
			K->Get_Center(x, y);
			double radius = sqrt(kv(x) + kv(y));
			if (radius < R11_)
			{
				continue;
			}
			double Volume = K->Get_Volume(); // K->Get_Volume_rotate(alpha_rot);
			//cout << "Setka.cpp    " << "Vol = " << Volume << endl;
			double Volume2 = K->Get_Volume_posle(); // K->Get_Volume_posle_rotate(alpha_rot);
			//cout << "Setka.cpp    " << "Vol2 = " << Volume2 << endl;
			double W = 0.0;
			bool np = true;

			double q2_1, q2_2, q3;
			q2_1 = 0.0;
			q2_2 = 0.0;
			q3 = 0.0;

			for (auto& i : K->Grans)
			{
				double x2, y2, x4, y4;
				//double S = i->Get_square_rotate(alpha_rot); 
				double S = i->Get_square();
				//cout << "Setka.cpp    " << "S = " << S << endl;
				i->Get_Center(x2, y2);
				i->Get_Center_posle(x4, y4);
				dist = sqrt(kv(x - x2) + kv(y - y2));
				Parametr par2, par11;
				par11 = par1;
				i->Get_par(par2, now1);
				i->Get_normal(n1, n2);
				double Vc, Vl, Vp;
				W = ((x4 - x2) * n1 + (y4 - y2) * n2) / T[now1];
				int met = i->metod_HLLC;
				double dis = sqrt(kv(x2) + kv(y2));


				if (i->type == Usualy)
				{
					double x3, y3;
					i->Sosed->Get_Center(x3, y3);
					double dd = sqrt(kv(x3) + kv(y3));

					if (K->type == C_1)
					{
						par11.ro = par11.ro * kv(radius) / kv(dis);
						par11.Q = par11.Q * kv(radius) / kv(dis);
						par11.p = par11.p * pow(radius / dis, 2 * ggg);
						polar_perenos(x, y, x2, y2, par11.u, par11.v);

						Parametr pp1;
						i->Get_par_TVD_radial(pp1, now1);
						par11.ro = pp1.ro;
						par11.p = pp1.p;
					}
					else
					{
						i->Get_par_TVD(par11, now1);
					}

					if (i->Sosed->type == C_1)
					{
						par2.ro = par2.ro * kv(dd) / kv(dis);
						par2.Q = par2.Q * kv(dd) / kv(dis);
						par2.p = par2.p * pow(dd / dis, 2 * ggg);
						polar_perenos(x3, y3, x2, y2, par2.u, par2.v);
						polar_perenos(x3, y3, x2, y2, par2.u_H1, par2.v_H1);

						Parametr pp2;
						i->Gran_copy->Get_par_TVD_radial(pp2, now1);
						par2.ro = pp2.ro;
						par2.p = pp2.p;
					}
					else
					{
						i->Gran_copy->Get_par_TVD(par2, now1);
					}
				}
				else if (i->type == Inner_sphere)
				{
					par11.ro = par11.ro * kv(radius) / kv(dis);
					par11.Q = par11.Q * kv(radius) / kv(dis);
					par11.p = par11.p * pow(radius / dis, 2 * ggg);
				}
				else if (K->type == C_1 && i->type == Axis)
				{
					par11.ro = par11.ro * kv(radius) / kv(dis);
					par11.Q = par11.Q * kv(radius) / kv(dis);
					par11.p = par11.p * pow(radius / dis, 2 * ggg);
					polar_perenos(x, y, x2, y2, par11.u, par11.v);
					//par11.v = 0.0;
					par2 = par11;
					par2.v = -par2.v;
				}
				else if (K->type != C_1 && i->type == Axis)
				{
					i->Get_par_TVD(par11, now1);
					par2 = par11;
					par2.v = -par2.v;
				}


				if ((i->A->type == P_Contact && i->B->type != P_Contact)||(i->A->type != P_Contact && i->B->type == P_Contact))
				{
					met = 0;
				}

				if ((i->A->type == P_Inner_shock && i->B->type != P_Inner_shock) || (i->A->type != P_Inner_shock && i->B->type == P_Inner_shock))
				{
					met = 0;
				}

				np = true;  // Сглаживание
				if ((i->A->type == P_Contact && i->B->type == P_Contact))
				{
					np = false;
					met = 1;
				}
				//np = true;  // Сглаживание

				if (x < -2.0 && y < 0.5)
				{
					met = 0;
				}
				
				//met = 0;
				
				if (false)//(K->type != C_centr && np == false)
				{
					vector<double> qqq1(5);
					vector<double> qqq2(5);
					vector<double> qqq(5);
					vector<double> n(3);

					auto SS = Solvers();
					n[0] = n1;
					n[1] = n2;
					n[2] = 0.0;

					qqq1[0] = par11.ro;
					qqq1[1] = par11.p;
					qqq1[2] = par11.u;
					qqq1[3] = par11.v;
					qqq1[4] = 0.0;

					qqq2[0] = par2.ro;
					qqq2[1] = par2.p;
					qqq2[2] = par2.u;
					qqq2[3] = par2.v;
					qqq2[4] = 0.0;


					time = min(time, SS.Godunov_Solver_Alexashov(qqq1, qqq2, n, qqq, Vl, Vp, Vc, W, dist));
					PQ = 0.0;
					P[0] = qqq[0];
					P[1] = qqq[1];
					P[2] = qqq[2];
					P[3] = qqq[3];
					for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
					{
						K->Potok[k] = K->Potok[k] + P[k] * S;
					}
					K->Potok[4] = K->Potok[4] + PQ * S;

				}
				else if (K->type != C_centr)
				{
					
					time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro, par11.Q, par11.p, par11.u, par11.v, par2.ro, par2.Q, //
						par2.p, par2.u, par2.v, W, P, PQ, n1, n2, dist, met, Vl, Vc, Vp, np));
					for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
					{
						K->Potok[k] = K->Potok[k] + P[k] * S;
					}
					K->Potok[4] = K->Potok[4] + PQ * S;
				}


				// Для источников при Kn = 0
				if (i->type == Usualy && i->Sosed->type == C_5 && (i->A->type == P_Outer_shock && i->B->type == P_Outer_shock))//(i->Sosed->type == C_4 || i->Sosed->type == C_5)
				{
					double Vn_0 = -(n1 * par2.u + n2 * par2.v);
					double c_2 = sqrt(par2.p / par2.ro);
					double k_0 = 0.5 * (1.0 + erf(Vn_0 / c_2));

					//Vn_0 = Vn_0 + c_2 * exp(-kv(Vn_0 / c_2)) / (sqrtpi_ * (1.0 + erf(Vn_0 / c_2)));

					//q2_1 = n_p_LISM_ * k_0 * par2.ro * (par2.u - par11.u) * Vn_0 * S / Volume2;
					//q2_2 = n_p_LISM_ * k_0 * par2.ro * (par2.v - par11.v) * Vn_0 * S / Volume2;
					//q3 = n_p_LISM_ * 0.5 * k_0 * par2.ro * (kvv(par2.u, par2.v, 0.0) - kvv(par11.u, par11.v, 0.0)) * Vn_0 * S / Volume2;
				
				}

				//if (i->type == Usualy && i->Sosed->type == C_4 && (i->A->type == P_Outer_shock && i->B->type == P_Outer_shock))//(i->Sosed->type == C_4 || i->Sosed->type == C_5)
				//{
				//	double Vn_0 = -(n1 * par2.u + n2 * par2.v);
				//	double c_2 = sqrt(par2.p / par2.ro);
				//	double k_0 = 0.5 * (1.0 + erf(Vn_0 / c_2));

				//	q2_1 = n_p_LISM_ * k_0 * par2.ro * (-c_2 * n1 - par11.u) * c_2 * S / Volume2;
				//	q2_2 = n_p_LISM_ * k_0 * par2.ro * (-c_2 * n2 - par11.v) * c_2 * S / Volume2;
				//	q3 = n_p_LISM_ * 0.5 * k_0 * par2.ro * (kvv(c_2 * n1, c_2 * n2, 0.0) - kvv(par11.u, par11.v, 0.0)) * c_2 * S / Volume2;
				//}


				if (i->type == Usualy && i->Sosed->type == C_4 && (i->A->type == P_Contact && i->B->type == P_Contact))//(i->Sosed->type == C_4 || i->Sosed->type == C_5)
				{
					double Vn_0 = 0.0; // -(n1 * par2.u + n2 * par2.v);
					double c_2 = sqrt(par2.p / par2.ro);
					double k_0 = 0.5 * (1 + erf(Vn_0 / c_2));

					Vn_0 = Vn_0 + c_2 * exp(-kv(Vn_0 / c_2)) / (sqrtpi_ * (1.0 + erf(Vn_0 / c_2)));

					//q2_1 = n_p_LISM_ * k_0 * par2.ro * (-Vn_0 * n1 - par11.u * (chi_real / chi_)) * Vn_0 * S / Volume2;
					//q2_2 = n_p_LISM_ * k_0 * par2.ro * (-Vn_0 * n2 - par11.v * (chi_real / chi_)) * Vn_0 * S / Volume2;
					//q3 = n_p_LISM_ * 0.5 * k_0 * par2.ro * (kvv(Vn_0, 0.0, 0.0) - kvv(par11.u * (chi_real / chi_), par11.v * (chi_real / chi_), 0.0)) * Vn_0 * S / Volume2 / (chi_real / chi_);
				
				}

				if (i->type == Usualy && i->Sosed->type == C_3 && (i->A->type == P_Contact && i->B->type == P_Contact))//(i->Sosed->type == C_4 || i->Sosed->type == C_5)
				{
					double Vn_0 = 0.0; //-(n1 * par2.u + n2 * par2.v);
					double Vn_1 = 0.0;// (n1* par11.u + n2 * par11.v);
					double c_2 = sqrt(par2.p / par2.ro) * (chi_real / chi_);
					double c_1 = sqrt(par11.p / par11.ro);
					double k_0 = 0.5 * (1 + erf(Vn_0 / c_2));
					double k_1 = 0.5 * (1 + erf(Vn_1 / c_1));

					Vn_0 = Vn_0 + c_2 * exp(-kv(Vn_0 / c_2)) / (sqrtpi_ * (1.0 + erf(Vn_0 / c_2)));
					Vn_1 = Vn_1 + c_1 * exp(-kv(Vn_1 / c_1)) / (sqrtpi_ * (1.0 + erf(Vn_1 / c_1)));

					/*q2_1 = 0.0 * n_p_LISM_ * k_1 * k_0 * par1.ro * (-Vn_0 * n1 - par11.u) * Vn_1 * S / Volume2;
					q2_2 = 0.0 * n_p_LISM_ * k_1 * k_0 * par1.ro * (-Vn_0 * n2 - par11.v) * Vn_1 * S / Volume2;
					q3 = 0.0 * n_p_LISM_ * 0.5 * k_1 * k_0 * par1.ro * (kvv(Vn_0, 0.0, 0.0) - kvv(par11.u, par11.v, 0.0)) * Vn_1 * S / Volume2;
					*/
					//cout << x2 << "  " << q2_1 << "  " << q2_2 << "  " << q3 << endl;
				
				}
			}

			if (time < T[now2])
			{
				mut.lock();
				T[now2] = time;
				mut.unlock();
			}

			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			if (false)
			{
				double U_M_H1, U_M_H2, U_M_H3, U_M_H4;
				double U_H1, U_H2, U_H3, U_H4;
				double sigma_H1, sigma_H2, sigma_H3, sigma_H4;
				double nu_H1, nu_H2, nu_H3, nu_H4;
				double q2_1, q2_2, q3;
				double u, v, ro, p, Q;

				double u_H1 = K->par[0].H_u[0], v_H1 = K->par[0].H_v[0], ro_H1 = K->par[0].H_n[0], p_H1 = 0.5 * K->par[0].H_T[0] * K->par[0].H_n[0];
				double u_H2 = K->par[0].H_u[1], v_H2 = K->par[0].H_v[1], ro_H2 = K->par[0].H_n[1], p_H2 = 0.5 * K->par[0].H_T[1] * K->par[0].H_n[1];
				double u_H3 = K->par[0].H_u[2], v_H3 = K->par[0].H_v[2], ro_H3 = K->par[0].H_n[2], p_H3 = 0.5 * K->par[0].H_T[2] * K->par[0].H_n[2];
				double u_H4 = K->par[0].H_u[3], v_H4 = K->par[0].H_v[3], ro_H4 = K->par[0].H_n[3], p_H4 = 0.5 * K->par[0].H_T[3] * K->par[0].H_n[3];


				if (false)//(par1.Q / par1.ro < 90)// (K->type == C_centr || K->type == C_1 || K->type == C_2 || K->type == C_3)
				{
					u = par1.u * (chi_real / chi_);
					v = par1.v * (chi_real / chi_);
					ro = par1.ro / kv(chi_real / chi_);
					p = par1.p;
					Q = par1.Q / kv(chi_real / chi_);
				}
				else
				{
					u = par1.u;
					v = par1.v;
					ro = par1.ro;
					p = par1.p;
					Q = par1.Q;

					/*if (K->pui_ == true)
					{
						ro = par1.ro - K->par[0].npui;
						p = K->par[0].pp;
					}
					else
					{
						ro = par1.ro;
						p = par1.p;
					}*/
				}

				if (ro <= 0.0)
				{
					ro = 0.0000001;
					p = 0.0;
				}
				if (ro_H1 <= 0.0)
				{
					ro_H1 = 0.0000001;
					p_H1 = 0.0;
				}
				if (ro_H2 <= 0.0)
				{
					ro_H2 = 0.0000001;
					p_H2 = 0.0;
				}
				if (ro_H3 <= 0.0)
				{
					ro_H3 = 0.0000001;
					p_H3 = 0.0;
				}
				if (ro_H4 <= 0.0)
				{
					ro_H4 = 0.0000001;
					p_H4 = 0.0;
				}

				U_M_H1 = sqrt(kv(u - u_H1) + kv(v - v_H1) + (64.0 / (9.0 * pi_)) //
					* (p / ro + 2.0 * p_H1 / ro_H1));
				U_M_H2 = sqrt(kv(u - u_H2) + kv(v - v_H2) + (64.0 / (9.0 * pi_)) //
					* (p / ro + 2.0 * p_H2 / ro_H2));
				U_M_H3 = sqrt(kv(u - u_H3) + kv(v - v_H3) + (64.0 / (9.0 * pi_)) //
					* (p / ro + 2.0 * p_H3 / ro_H3));
				U_M_H4 = sqrt(kv(u - u_H4) + kv(v - v_H4) + (64.0 / (9.0 * pi_)) //
					* (p / ro + 2.0 * p_H4 / ro_H4));

				U_H1 = sqrt(kv(u - u_H1) + kv(v - v_H1) + (4.0 / pi_) //
					* (p / ro + 2.0 * p_H1 / ro_H1));
				U_H2 = sqrt(kv(u - u_H2) + kv(v - v_H2) + (4.0 / pi_) //
					* (p / ro + 2.0 * p_H2 / ro_H2));
				U_H3 = sqrt(kv(u - u_H3) + kv(v - v_H3) + (4.0 / pi_) //
					* (p / ro + 2.0 * p_H3 / ro_H3));
				U_H4 = sqrt(kv(u - u_H4) + kv(v - v_H4) + (4.0 / pi_) //
					* (p / ro + 2.0 * p_H4 / ro_H4));

				sigma_H1 = kv(1.0 - a_2 * log(U_M_H1)); // 0.1243
				sigma_H2 = kv(1.0 - a_2 * log(U_M_H2));
				sigma_H3 = kv(1.0 - a_2 * log(U_M_H3)); // 0.1121     a_2
				sigma_H4 = kv(1.0 - a_2 * log(U_M_H4));

				nu_H1 = ro * ro_H1 * U_M_H1 * sigma_H1;
				nu_H2 = ro * ro_H2 * U_M_H2 * sigma_H2;
				nu_H3 = ro * ro_H3 * U_M_H3 * sigma_H3;
				nu_H4 = ro * ro_H4 * U_M_H4 * sigma_H4;

				K->par[0].M_u = (n_p_LISM_ / Kn_) * (nu_H1 * (u_H1 - u) + nu_H2 * (u_H2 - u) //
					+ nu_H3 * (u_H3 - u) + nu_H4 * (u_H4 - u));
				K->par[0].M_v = (n_p_LISM_ / Kn_) * (nu_H1 * (v_H1 - v) + nu_H2 * (v_H2 - v) //
					+ nu_H3 * (v_H3 - v) + nu_H4 * (v_H4 - v));

				if (false)//(par1.Q / par1.ro < 90)//(K->type == C_centr || K->type == C_1 || K->type == C_2 || K->type == C_3)//
				{
					// Нужно изменить источник, так как он подставляется в уравнения, где параметры плазмы перенормированы 
					K->par[0].M_T = (n_p_LISM_ / Kn_) * (nu_H1 * ((kv(u_H1) + kv(v_H1) - kv(u) - kv(v)) / 2.0 + //
						(U_H1 / U_M_H1) * (2.0 * p_H1 / ro_H1 - p / ro)) + //
						nu_H2 * ((kv(u_H2) + kv(v_H2) - kv(u) - kv(v)) / 2.0 + //
							(U_H2 / U_M_H2) * (2.0 * p_H2 / ro_H2 - p / ro)) + //
						nu_H3 * ((kv(u_H3) + kv(v_H3) - kv(u) - kv(v)) / 2.0 + //
							(U_H3 / U_M_H3) * (2.0 * p_H3 / ro_H3 - p / ro)) + //
						nu_H4 * ((kv(u_H4) + kv(v_H4) - kv(u) - kv(v)) / 2.0 + //
							(U_H4 / U_M_H4) * (2.0 * p_H4 / ro_H4 - p / ro))) / (chi_real / chi_);
				}
				else
				{
					K->par[0].M_T = (n_p_LISM_ / Kn_) * (nu_H1 * ((kv(u_H1) + kv(v_H1) - kv(u) - kv(v)) / 2.0 + //
						(U_H1 / U_M_H1) * (2.0 * p_H1 / ro_H1 - p / ro)) + //
						nu_H2 * ((kv(u_H2) + kv(v_H2) - kv(u) - kv(v)) / 2.0 + //
							(U_H2 / U_M_H2) * (2.0 * p_H2 / ro_H2 - p / ro)) + //
						nu_H3 * ((kv(u_H3) + kv(v_H3) - kv(u) - kv(v)) / 2.0 + //
							(U_H3 / U_M_H3) * (2.0 * p_H3 / ro_H3 - p / ro)) + //
						nu_H4 * ((kv(u_H4) + kv(v_H4) - kv(u) - kv(v)) / 2.0 + //
							(U_H4 / U_M_H4) * (2.0 * p_H4 / ro_H4 - p / ro)));
				}

				q2_1 = K->par[0].M_u * K->par[0].k_u;
				q2_2 = K->par[0].M_v * K->par[0].k_v;
				q3 = K->par[0].M_T * K->par[0].k_T;

				//cout << "Setka.cpp    " << q2_1 << " " << q2_2 << " " << q3 << endl;

				u = par1.u;
				v = par1.v;
				ro = par1.ro;
				p = par1.p;
				Q = par1.Q;

				double ro3, Q33, u3, v3, p3;
			}

			double u = par1.u;
			double v = par1.v;
			double ro = par1.ro;
			double p = par1.p;
			double Q = par1.Q;

			double a1, a2, a3;
			

			if (par1.Q / par1.ro < Q_gran)// (K->type == C_centr || K->type == C_1 || K->type == C_2 || K->type == C_3)
			{
				u = par1.u * (chi_real / chi_);
				v = par1.v * (chi_real / chi_);
				ro = par1.ro / kv(chi_real / chi_);
				p = par1.p;
				Q = par1.Q / kv(chi_real / chi_);
			}


			K->Get_Sourse_MK1(a1, a2, a3, u, v, ro, p, true);
			q2_1 = a1;
			q2_2 = a2;
			q3 = a3;


			if (par1.Q / par1.ro < Q_gran)// (K->type == C_centr || K->type == C_1 || K->type == C_2 || K->type == C_3)
			{
				q3 = q3 / (chi_real / chi_);
			}
			
			u = par1.u;
			v = par1.v;
			ro = par1.ro;
			p = par1.p;
			Q = par1.Q;

			//q2_1 = q2_2 = q3 = 0.0;

			double ro3, Q33, u3, v3, p3;

			if (K->type != C_centr)
			{
				ro3 = ro * Volume / Volume2 - T[now1] * (K->Potok[0] / Volume2 + ro * v / y);
				Q33 = Q * Volume / Volume2 - (T[now1] / Volume2) * K->Potok[4] - T[now1] * Q * v / y;
				if (ro3 <= 0)
				{
					printf("Problemsssss  ro < 0!   %lf,  %lf\n", x, y);
					ro3 = 0.001;
				}
				u3 = (ro * u * Volume / Volume2 - T[now1] * (K->Potok[1] / Volume2 + ro * v * u / y - q2_1)) / ro3;
				v3 = (ro * v * Volume / Volume2 - T[now1] * (K->Potok[2] / Volume2 + ro * v * v / y - q2_2)) / ro3;
				p3 = (((p / (ggg - 1) + ro * (u * u + v * v) * 0.5) * Volume / Volume2 - T[now1] * (K->Potok[3] / Volume2 + //
					+v * (ggg * p / (ggg - 1) + ro * (u * u + v * v) * 0.5) / y - q3)) - //
					0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);
				if (p3 <= 0)
				{
					cout << "problem p  " << x << " " << y << " " << p3 << endl;
					p3 = 0.000001;
				}

				

				K->par[now2].ro = ro3;
				K->par[now2].Q = Q33;
				K->par[now2].p = p3;
				K->par[now2].u = u3;
				K->par[now2].v = v3;
			}
			else
			{
				K->par[now2].ro = K->par[now1].ro;
				K->par[now2].Q = K->par[now1].Q;
				K->par[now2].p = K->par[now1].p;
				K->par[now2].u = K->par[now1].u;
				K->par[now2].v = K->par[now1].v;
			}

		}

		if (dvig == true)
		{
			for (auto& i : this->All_Points)
			{
				i->x = i->x2;
				i->y = i->y2;
				i->Vx = 0.0;
				i->Vy = 0.0;
				i->count = 0;
			}
		}
	}
}

void Setka::Go_5_komponent__MK2(int step)
// Идея в том, чтобы использовать расчёты в угле (3д задача, сведённая к 2д)
// именно этим программа отличается от первой версии
{
	//cout << "Setka.cpp    " << "START " << step << endl;
	int now1 = 1;
	int now2 = 0;
	double T[2];
	T[0] = T[1] = 0.00000001;
	mutex mut;

	for (int st = 0; st < step; st++)
	{
		//cout << "Setka.cpp    " << st << endl;
		if (st % 10000 == 0 && st > 0)
		{
			cout << "Setka.cpp    " << st << " " << T[now2] << endl;
		}
		now1 = (now1 + 1) % 2; // Какие параметры сейчас берём
		now2 = (now2 + 1) % 2; // Какие параметры сейчас меняем
		double time = 10000000;
		T[now2] = 100000000;

		this->Move_surface(now1, T[now1]);
		this->Move_Setka_Calculate_2(T[now1]);

		//omp_set_num_threads(4);
#pragma omp parallel for
		for (int num_cell = 0; num_cell < this->All_Cells.size(); num_cell++)
		{
			auto K = this->All_Cells[num_cell];
			K->Potok[0] = K->Potok[1] = K->Potok[2] = K->Potok[3] = K->Potok[4] = 0.0;

			double P[4];
			double PP[8];
			double PQ;
			Parametr par1 = K->par[now1];
			double n1, n2;
			double dist;
			double x, y;
			K->Get_Center(x, y);
			double vs1 = K->Get_Volume();
			double vs2 = K->Get_Volume_posle();

			//y = K->Get_Volume_rotate(360.0) / (2.0 * pi_ * vs1);
			//double y2 = K->Get_Volume_posle_rotate(360.0) / (2.0 * pi_ * vs2);

			double radius = sqrt(kv(x) + kv(y));
			if (radius < R11_)
			{
				continue;
			}

			if (false)//(radius < 65.0/RR_)                                                    // ДЛЯ ОБЫЧНОЙ ГАЗОВОЙ ДИНАМИКИ
			{
				double dist = sqrt(x * x + y * y);
				double ro = 1.0 / (kv(chi_real) * dist * dist);
				double P_E = ro * chi_real * chi_real / (ggg * 10.0 * 10.0);
				K->par[now2].ro = ro * kv(chi_real / chi_);
				K->par[now2].Q = ro * kv(chi_real / chi_);
				K->par[now2].p = P_E;
				K->par[now2].u = chi_real * x / dist / (chi_real / chi_);
				K->par[now2].v = chi_real * y / dist / (chi_real / chi_);
				continue;
			}

			double Volume = K->Get_Volume_rotate(360.0) / (2.0 * pi_);
			//cout << "Setka.cpp    " << "Vol = " << Volume << endl;
			double Volume2 = K->Get_Volume_posle_rotate(360.0) / (2.0 * pi_);
			//cout << "Setka.cpp    " << "Vol2 = " << Volume2 << endl;
			double Scw = K->Get_Volume();
			double W = 0.0;
			bool np = true;

			for (auto& i : K->Grans)
			{
				double x2, y2, x4, y4;
				double S = i->Get_square() * (i->A->y + i->B->y) * 0.5;
				//double S = i->Get_square();
				//cout << "Setka.cpp    " << "S = " << S << endl;
				i->Get_Center(x2, y2);
				i->Get_Center_posle(x4, y4);
				dist = sqrt(kv(x - x2) + kv(y - y2));
				Parametr par2, par11;
				par11 = par1;
				i->Get_par(par2, now1);
				i->Get_normal(n1, n2);
				double Vc, Vl, Vp;
				W = ((x4 - x2) * n1 + (y4 - y2) * n2) / T[now1];
				int met = 0;  // 1                                                                 метод распадника
				double dis = sqrt(kv(x2) + kv(y2));

				if (i->type == Usualy)
				{
					double x3, y3;
					i->Sosed->Get_Center(x3, y3);
					double dd = sqrt(kv(x3) + kv(y3));

					if (K->type == C_1)
					{
						par11.ro = par11.ro * kv(radius) / kv(dis);
						par11.Q = par11.Q * kv(radius) / kv(dis);
						par11.p = par11.p * pow(radius / dis, 2 * ggg);
						polar_perenos(x, y, x2, y2, par11.u, par11.v);
						polar_perenos(x, y, x2, y2, par11.u_H1, par11.v_H1);

						Parametr pp1;
						i->Get_par_TVD_radial(pp1, now1);
						par11.ro = pp1.ro;
						par11.p = pp1.p;
					}
					else
					{
						i->Get_par_TVD(par11, now1);
					}

					if (i->Sosed->type == C_1)
					{
						par2.ro = par2.ro * kv(dd) / kv(dis);
						par2.Q = par2.Q * kv(dd) / kv(dis);
						par2.p = par2.p * pow(dd / dis, 2 * ggg);
						polar_perenos(x3, y3, x2, y2, par2.u, par2.v);
						polar_perenos(x3, y3, x2, y2, par2.u_H1, par2.v_H1);

						Parametr pp2;
						i->Gran_copy->Get_par_TVD_radial(pp2, now1);
						par2.ro = pp2.ro;
						par2.p = pp2.p;
					}
					else
					{
						i->Gran_copy->Get_par_TVD(par2, now1);
					}
				}
				else if (i->type == Inner_sphere)
				{
					par11.ro = par11.ro * kv(radius) / kv(dis);
					par11.Q = par11.Q * kv(radius) / kv(dis);
					par11.p = par11.p * pow(radius / dis, 2 * ggg);
				}
				else if (K->type == C_1 && i->type == Axis)
				{
					par11.ro = par11.ro * kv(radius) / kv(dis);
					par11.Q = par11.Q * kv(radius) / kv(dis);
					par11.p = par11.p * pow(radius / dis, 2 * ggg);
					polar_perenos(x, y, x2, y2, par11.u, par11.v);
					//par11.v = 0.0;
					par2 = par11;
					par2.v = -par2.v;
				}
				else if (K->type != C_1 && i->type == Axis)
				{
					i->Get_par_TVD(par11, now1);
					par2 = par11;
					par2.v = -par2.v;
				}


				/*if ((i->A->type == P_Contact && i->B->type != P_Contact)||(i->A->type != P_Contact && i->B->type == P_Contact))
				{
					met = 0;
				}*/

				np = false;
				bool god = false;

				if ((i->A->type == P_Contact && i->B->type == P_Contact))
				{
					god = true;
				}
				else if ((i->A->type == P_Inner_shock && i->B->type == P_Inner_shock))
				{
					god = true;
				}
				else if ((i->A->type == P_Outer_shock && i->B->type == P_Outer_shock))
				{
					god = true;
				}


				if ((i->A->type == P_Contact && i->B->type != P_Contact))
				{
					np = true;
				}
				else if ((i->A->type != P_Contact && i->B->type == P_Contact))
				{
					np = true;
				}
				else if ((i->A->type == P_Inner_shock && i->B->type != P_Inner_shock))
				{
					np = true;
				}
				else if ((i->A->type != P_Inner_shock && i->B->type == P_Inner_shock))
				{
					np = true;
				}
				else if ((i->A->type == P_Outer_shock && i->B->type != P_Outer_shock))
				{
					np = true;
				}
				else if ((i->A->type != P_Outer_shock && i->B->type == P_Outer_shock))
				{
					np = true;
				}

				//np = true;

				/*if (y < 30 && x > 120 && x < 160)
				{
					np = true;
				}*/

				//god = false;

				if (god == true)
				{
					double Vl, Vp, VV;
					vector<double> qqq1(5);
					vector<double> qqq2(5);
					vector<double> qqq(5);
					vector<double> n(3);
					n[0] = n1;
					n[1] = n2;
					n[2] = 0.0;

					qqq1[0] = par11.ro;
					qqq1[1] = par11.p;
					qqq1[2] = par11.u;
					qqq1[3] = par11.v;
					qqq1[4] = 0.0;

					qqq2[0] = par2.ro;
					qqq2[1] = par2.p;
					qqq2[2] = par2.u;
					qqq2[3] = par2.v;
					qqq2[4] = 0.0;

					auto SS = Solvers();
					time = min(time, SS.Godunov_Solver_Alexashov(qqq1, qqq2, n, qqq, Vl, Vp, VV, W, dist));

					K->Potok[0] = K->Potok[0] + qqq[0] * S;
					K->Potok[1] = K->Potok[1] + qqq[1] * S;
					K->Potok[2] = K->Potok[2] + qqq[2] * S;
					K->Potok[3] = K->Potok[3] + qqq[4] * S;

					time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro, par11.Q, par11.p, par11.u, par11.v, par2.ro, par2.Q, //
						par2.p, par2.u, par2.v, W, P, PQ, n1, n2, dist, met, Vl, Vc, Vp, np));

					K->Potok[4] = K->Potok[4] + PQ * S;
				}
				else
				{
					time = min(time, this->HLLC_2d_Korolkov_b_s(par11.ro, par11.Q, par11.p, par11.u, par11.v, par2.ro, par2.Q, //
						par2.p, par2.u, par2.v, W, P, PQ, n1, n2, dist, met, Vl, Vc, Vp, np));
					for (int k = 0; k < 4; k++)  // Суммируем все потоки в ячейке
					{
						K->Potok[k] = K->Potok[k] + P[k] * S;
					}
					K->Potok[4] = K->Potok[4] + PQ * S;
				}

			}



			// Блок законов сохранения в ячейке

			if (time < T[now2])
			{
				mut.lock();
				T[now2] = time;
				mut.unlock();
			}

			double u = par1.u;
			double v = par1.v;
			double ro = par1.ro;
			double p = par1.p;
			double Q = par1.Q;
			
			double a1, a2, a3;
			double b1, b2, b3;
			double c1, c2, c3;
			double q2_1, q2_2, q3;

			if (polusum == true)
			{
				K->Get_Sourse_MK1(a1, a2, a3, u, v, ro, p);
				K->Get_Sourse_MK2(b1, b2, b3, u, v, ro, p);
				K->Get_Sourse_MK3(c1, c2, c3, u, v, ro, p);
				q2_1 = (a1 + b1 + c1) / 3.0;
				q2_2 = (a2 + b2 + c2) / 3.0;
				q3 = (a3 + b3 + c3) / 3.0;
			}
			else
			{
				K->Get_Sourse_MK1(a1, a2, a3, u, v, ro, p);
				q2_1 = a1;
				q2_2 = a2;
				q3 = a3;
			}

			double ro3, Q33, u3, v3, p3;

			if (K->type != C_centr)
			{
				ro3 = ro * Volume / Volume2 - T[now1] * (K->Potok[0] / Volume2);  // В декартовых
				//ro3 = ro * Volume / Volume2 - T[now1] * (K->Potok[0] / Volume2 + ro * v / y);  // В цилиндрических

				Q33 = Q * Volume / Volume2 - (T[now1] / Volume2) * K->Potok[4];
				//Q33 = Q * Volume / Volume2 - (T[now1] / Volume2) * K->Potok[4] - T[now1] * Q * v / y;
				if (ro3 <= 0.0)
				{
					printf("Problemsssss  ro < 0!    %lf, %lf\n", x, y);
					ro3 = 0.00001;
				}

				u3 = (ro * u * Volume / Volume2 - T[now1] * (K->Potok[1] / Volume2 - q2_1)) / ro3;
				v3 = (ro * v * Volume / Volume2 - T[now1] * ((K->Potok[2] - p * 0.5 * (vs1 + vs2)) / Volume2 - q2_2)) / ro3;
				//v3 = (ro * v * Volume / Volume2 - T[now1] * (K->Potok[2] / Volume2 - q2_2)) / ro3;
				p3 = (((p / (ggg - 1) + ro * (u * u + v * v) * 0.5) * Volume / Volume2 - T[now1] * (K->Potok[3] / Volume2 - q3)) - //
					0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);

				//u3 = (ro * u * Volume / Volume2 - T[now1] * (K->Potok[1] / Volume2 + ro * v * u / y - q2_1)) / ro3;
				//v3 = (ro * v * Volume / Volume2 - T[now1] * (K->Potok[2] / Volume2 + ro * v * v / y - q2_2)) / ro3;
				//p3 = (((p / (ggg - 1) + ro * (u * u + v * v) * 0.5) * Volume / Volume2 - T[now1] * (K->Potok[3] / Volume2 + //
				//	+v * (ggg * p / (ggg - 1) + ro * (u * u + v * v) * 0.5) / y - q3)) - //
				//	0.5 * ro3 * (u3 * u3 + v3 * v3)) * (ggg - 1);

				if (p3 <= 0)
				{
					p3 = 0.000001;
				}

				K->par[now2].ro = ro3;
				K->par[now2].Q = Q33;
				K->par[now2].p = p3;
				K->par[now2].u = u3;
				K->par[now2].v = v3;
			}
			else
			{
				K->par[now2].ro = K->par[now1].ro;
				K->par[now2].Q = K->par[now1].Q;
				K->par[now2].p = K->par[now1].p;
				K->par[now2].u = K->par[now1].u;
				K->par[now2].v = K->par[now1].v;
			}

		}


		for (auto& i : this->All_Points)
		{
			i->x = i->x2;
			i->y = i->y2;
			i->Vx = 0.0;
			i->Vy = 0.0;
			i->count = 0;
		}
	}
}

double Setka::HLLC_2d_Korolkov_b_s(const double& ro_L, const double& Q_L, const double& p_L, const double& v1_L, const double& v2_L,//
	const double& ro_R, const double& Q_R, const double& p_R, const double& v1_R, const double& v2_R, const double& W, //
	double* P, double& PQ, const double& n1, const double& n2, const double& rad, int metod, double& Vl, double& Vc, double& Vp, bool nul_potok)
	// BestSeries
	// Лучший работающий 2д распадник
	//
	//  Вывод:
	// P[1]       // Скорости
	// P[2]
	// P[0]       // Масса
	// P[3]       // Энергия
{
	double t1 = n2;
	double t2 = -n1;

	double u1, v1, u2, v2;
	u1 = v1_L * n1 + v2_L * n2;
	v1 = v1_L * t1 + v2_L * t2;
	u2 = v1_R * n1 + v2_R * n2;
	v2 = v1_R * t1 + v2_R * t2;

	double sqrtroL = sqrt(ro_L);
	double sqrtroR = sqrt(ro_R);
	double cL = sqrt(ggg * p_L / ro_L);
	double cR = sqrt(ggg * p_R / ro_R);

	double uu_L = (kv(v1_L) + kv(v2_L)) / 2.0;
	double uu_R = (kv(v1_R) + kv(v2_R)) / 2.0;


	double SL = min(u1, u2) - max(cL, cR);
	double SR = max(u1, u2) + max(cL, cR);

	//cout << "Sl ==  " << SL << " " << SR << endl;


	Vl = SL;
	Vp = SR;

	double suR = (SR - u2);
	double suL = (SL - u1);

	double SM = (suR * ro_R * u2 - suL * ro_L * u1 - p_R + p_L) //
		/ (suR * ro_R - suL * ro_L);

	Vc = SM;

	//cout << "SM = " << SM << endl;

	double PTT = (suR * ro_R * p_L - suL * ro_L * p_R + ro_L * ro_R * suR * suL * (u2 - u1)) / (suR * ro_R - suL * ro_L);

	double UU = max(fabs(SL), fabs(SR));
	double time = kurant * rad / UU;

	double FL[5], FR[5], UL[5], UR[5];

	double e1 = p_L / g1 + ro_L * uu_L;
	double e2 = p_R / g1 + ro_R * uu_R;


	FL[0] = ro_L * u1;
	FL[1] = ro_L * u1 * u1 + p_L;
	FL[2] = ro_L * u1 * v1;
	FL[3] = (e1 + p_L) * u1;
	FL[4] = Q_L * u1;

	FR[0] = ro_R * u2;
	FR[1] = ro_R * u2 * u2 + p_R;
	FR[2] = ro_R * u2 * v2;
	FR[3] = (e2 + p_R) * u2;
	FR[4] = Q_R * u2;


	UL[0] = ro_L;
	UL[1] = ro_L * u1;
	UL[2] = ro_L * v1;
	UL[3] = e1;
	UL[4] = Q_L;

	UR[0] = ro_R;
	UR[1] = ro_R * u2;
	UR[2] = ro_R * v2;
	UR[3] = e2;
	UR[4] = Q_R;

	if (SL >= W)
	{
		//cout << "Setka.cpp    " << "A = " << 1 << endl;
		P[1] = n1 * (FL[1] - W * UL[1]) + t1 * (FL[2] - W * UL[2]);     // Скорости
		P[2] = n2 * (FL[1] - W * UL[1]) + t2 * (FL[2] - W * UL[2]);
		P[0] = FL[0] - W * UL[0];                       // Масса
		P[3] = FL[3] - W * UL[3];                       // Энергия
		PQ = FL[4] - W * UL[4];
		return time;
	}

	if (SR <= W)
	{
		//cout << "Setka.cpp    " << "A = " << 4 << endl;
		P[1] = n1 * (FR[1] - W * UR[1]) + t1 * (FR[2] - W * UR[2]);     // Скорости
		P[2] = n2 * (FR[1] - W * UR[1]) + t2 * (FR[2] - W * UR[2]);
		P[0] = FR[0] - W * UR[0];                       // Масса
		P[3] = FR[3] - W * UR[3];                       // Энергия
		PQ = FR[4] - W * UR[4];
		return time;
	}

	double ro_LL = ro_L * (SL - u1) / (SL - SM);
	double ro_RR = ro_R * (SR - u2) / (SR - SM);
	double Q_LL = Q_L * (SL - u1) / (SL - SM);
	double Q_RR = Q_R * (SR - u2) / (SR - SM);


	double UZ0 = (SR * UR[0] - SL * UL[0] + FL[0] - FR[0]) / (SR - SL);
	double UZ1 = (SR * UR[1] - SL * UL[1] + FL[1] - FR[1]) / (SR - SL);
	double UZ2 = (SR * UR[2] - SL * UL[2] + FL[2] - FR[2]) / (SR - SL);
	double UZ3 = (SR * UR[3] - SL * UL[3] + FL[3] - FR[3]) / (SR - SL);
	double UZ4 = (SR * UR[4] - SL * UL[4] + FL[4] - FR[4]) / (SR - SL);
	double vzL, vzR, vLL, vRR, ppLR, ee1, ee2;

	if (metod == 0)
	{
		double  PO[5];
		for (int i = 0; i < 5; i++)
		{
			PO[i] = (SR * FL[i] - SL * FR[i] + SR * SL * (UR[i] - UL[i])) / (SR - SL);
		}

		P[1] = n1 * (PO[1] - W * UZ1) + t1 * (PO[2] - W * UZ2);     // Скорости
		P[2] = n2 * (PO[1] - W * UZ1) + t2 * (PO[2] - W * UZ2);
		P[0] = PO[0] - W * UZ0;                       // Масса
		P[3] = PO[3] - W * UZ3;                       // Энергия
		PQ = PO[4] - W * UZ4;
		return time;
	}


	double suRm = suR / (SR - SM);
	double suLm = suL / (SL - SM);
	double rzR = ro_R * suRm;
	double rzL = ro_L * suLm;

	double ptzR = p_R + ro_R * suR * (SM - u2);
	double ptzL = p_L + ro_L * suL * (SM - u1);
	double ptz = (ptzR + ptzL) / 2.0;



	/*if( fabs(v1 - v2) > 0.1)
	{
		vLL = v1;
		vRR = v2;
	}
	else
	{
		vRR = UZ2 / UZ0;
		vLL = vRR;
	}*/


	if (nul_potok == true)   // Некое сглаживание
	{
		vRR = UZ2 / UZ0;
		vLL = vRR;
	}
	else
	{
		vLL = v1;
		vRR = v2;
	}



	ee2 = e2 * suRm + (ptz * SM - p_R * u2) / (SR - SM);
	ee1 = e1 * suLm + (ptz * SM - p_L * u1) / (SL - SM);


	double  ULL[5], URR[5], PO[5];
	ULL[0] = ro_LL;
	ULL[1] = ro_LL * SM;
	ULL[2] = ro_LL * vLL;
	ULL[3] = ee1;
	ULL[4] = Q_LL;

	URR[0] = ro_RR;
	URR[1] = ro_RR * SM;
	URR[2] = ro_RR * vRR;
	URR[3] = ee2;
	URR[4] = Q_RR;

	//cout << ee1 << " " << ee2 << endl;

	if (SL < W && SM >= W)
	{
		//cout << "Setka.cpp    " << "A = " << 2 << endl;
		for (int i = 0; i < 5; i++)
		{
			PO[i] = FL[i] + SL * ULL[i] - SL * UL[i] - W * ULL[i];
		}
	}
	else if (SR > W && SM < W)
	{
		//cout << "Setka.cpp    " << "A = " <<  3 << endl;
		for (int i = 0; i < 5; i++)
		{
			PO[i] = FR[i] + SR * URR[i] - SR * UR[i] - W * URR[i];
		}
	}

	P[1] = n1 * PO[1] + t1 * PO[2];     // Скорости
	P[2] = n2 * PO[1] + t2 * PO[2];
	P[0] = PO[0];                       // Масса
	P[3] = PO[3];                       // Энергия
	PQ = PO[4];

	return time;
}

double Setka::HLLC_2d_Korolkov_b_s_2(const double& ro_L, const double& Q_L, const double& p_L, const double& v1_L, const double& v2_L,//
	const double& ro_R, const double& Q_R, const double& p_R, const double& v1_R, const double& v2_R, const double& W, //
	double* P, double& PQ, const double& n1, const double& n2, const double& rad, int metod, double& Vl, double& Vc, double& Vp, bool nul_potok)
	// BestSeries
	// Лучший работающий 2д распадник
	//
	//  Вывод:
	// P[1]       // Скорости
	// P[2]
	// P[0]       // Масса
	// P[3]       // Энергия
{
	double t1 = n2;
	double t2 = -n1;

	double u1, v1, u2, v2;
	u1 = v1_L * n1 + v2_L * n2;
	v1 = v1_L * t1 + v2_L * t2;
	u2 = v1_R * n1 + v2_R * n2;
	v2 = v1_R * t1 + v2_R * t2;

	double sqrtroL = sqrt(ro_L);
	double sqrtroR = sqrt(ro_R);
	double cL = sqrt(ggg * p_L / ro_L);
	double cR = sqrt(ggg * p_R / ro_R);


	double uu_L = (kv(v1_L) + kv(v2_L)) / 2.0;
	double uu_R = (kv(v1_R) + kv(v2_R)) / 2.0;



	double SL = min(u1, u2) - max(cL, cR);
	double SR = max(u1, u2) + max(cL, cR);

	Vl = SL;
	Vp = SR;

	double suR = (SR - u2);
	double suL = (SL - u1);

	double SM = (suR * ro_R * u2 - suL * ro_L * u1 - p_R + p_L) //
		/ (suR * ro_R - suL * ro_L);

	Vc = SM;

	double PTT = (suR * ro_R * p_L - suL * ro_L * p_R + ro_L * ro_R * suR * suL * (u2 - u1)) / (suR * ro_R - suL * ro_L);

	double UU = max(fabs(SL), fabs(SR));
	double time = kurant * rad / UU;

	double FL[5], FR[5], UL[5], UR[5];

	double e1 = p_L / g1 + ro_L * uu_L;
	double e2 = p_R / g1 + ro_R * uu_R;


	FL[0] = ro_L * u1;
	FL[1] = ro_L * u1 * u1 + p_L;
	FL[2] = ro_L * u1 * v1;
	FL[3] = (e1 + p_L) * u1;
	FL[4] = Q_L * u1;

	FR[0] = ro_R * u2;
	FR[1] = ro_R * u2 * u2 + p_R;
	FR[2] = ro_R * u2 * v2;
	FR[3] = (e2 + p_R) * u2;
	FR[4] = Q_R * u2;


	UL[0] = ro_L;
	UL[1] = ro_L * u1;
	UL[2] = ro_L * v1;
	UL[3] = e1;
	UL[4] = Q_L;

	UR[0] = ro_R;
	UR[1] = ro_R * u2;
	UR[2] = ro_R * v2;
	UR[3] = e2;
	UR[4] = Q_R;

	if (SL >= W)
	{
		P[1] = n1 * (FL[1] - W * UL[1]) + t1 * (FL[2] - W * UL[2]);     // Скорости
		P[2] = n2 * (FL[1] - W * UL[1]) + t2 * (FL[2] - W * UL[2]);
		P[0] = FL[0] - W * UL[0];                       // Масса
		P[3] = FL[3] - W * UL[3];                       // Энергия
		PQ = FL[4] - W * UL[4];
		return time;
	}

	if (SR <= W)
	{
		P[1] = n1 * (FR[1] - W * UR[1]) + t1 * (FR[2] - W * UR[2]);     // Скорости
		P[2] = n2 * (FR[1] - W * UR[1]) + t2 * (FR[2] - W * UR[2]);
		P[0] = FR[0] - W * UR[0];                       // Масса
		P[3] = FR[3] - W * UR[3];                       // Энергия
		PQ = FR[4] - W * UR[4];
		return time;
	}


	double ro_LL = ro_L * (SL - u1) / (SL - SM);
	double ro_RR = ro_R * (SR - u2) / (SR - SM);
	double Q_LL = Q_L * (SL - u1) / (SL - SM);
	double Q_RR = Q_R * (SR - u2) / (SR - SM);


	double UZ0 = (SR * UR[0] - SL * UL[0] + FL[0] - FR[0]) / (SR - SL);
	double UZ1 = (SR * UR[1] - SL * UL[1] + FL[1] - FR[1]) / (SR - SL);
	double UZ2 = (SR * UR[2] - SL * UL[2] + FL[2] - FR[2]) / (SR - SL);
	double UZ3 = (SR * UR[3] - SL * UL[3] + FL[3] - FR[3]) / (SR - SL);
	double UZ4 = (SR * UR[4] - SL * UL[4] + FL[4] - FR[4]) / (SR - SL);
	double vzL, vzR, vLL, vRR, ppLR, ee1, ee2;

	if (metod == 0)
	{
		double  PO[5];
		for (int i = 0; i < 5; i++)
		{
			PO[i] = (SR * FL[i] - SL * FR[i] + SR * SL * (UR[i] - UL[i])) / (SR - SL);
		}

		P[1] = n1 * (PO[1] - W * UZ1) + t1 * (PO[2] - W * UZ2);     // Скорости
		P[2] = n2 * (PO[1] - W * UZ1) + t2 * (PO[2] - W * UZ2);
		P[0] = PO[0] - W * UZ0;                       // Масса
		P[3] = PO[3] - W * UZ3;                       // Энергия
		PQ = PO[4] - W * UZ4;
		return time;
	}


	double suRm = suR / (SR - SM);
	double suLm = suL / (SL - SM);
	double rzR = ro_R * suRm;
	double rzL = ro_L * suLm;

	double ptzR = p_R + ro_R * suR * (SM - u2);
	double ptzL = p_L + ro_L * suL * (SM - u1);
	double ptz = (ptzR + ptzL) / 2.0;


	/*if( fabs(v1 - v2) > 0.1)
	{
		vLL = v1;
		vRR = v2;
	}
	else
	{
		vRR = UZ2 / UZ0;
		vLL = vRR;
	}*/


	if (nul_potok == true)   // Некое сглаживание
	{
		vRR = UZ2 / UZ0;
		vLL = vRR;
	}
	else
	{
		vLL = v1;
		vRR = v2;
	}



	ee2 = e2 * suRm + (ptz * SM - p_R * u2) / (SR - SM);
	ee1 = e1 * suLm + (ptz * SM - p_L * u1) / (SL - SM);


	double  ULL[5], URR[5], PO[5];
	ULL[0] = ro_LL;
	ULL[1] = ro_LL * SM;
	ULL[2] = ro_LL * vLL;
	ULL[3] = ee1;
	ULL[4] = Q_LL;

	URR[0] = ro_RR;
	URR[1] = ro_RR * SM;
	URR[2] = ro_RR * vRR;
	URR[3] = ee2;
	URR[4] = Q_RR;

	if (SL < W && SM >= W)
	{
		for (int i = 0; i < 5; i++)
		{
			PO[i] = FL[i] + SL * ULL[i] - SL * UL[i] - W * ULL[i];
		}
	}
	else if (SR > W && SM < W)
	{
		for (int i = 0; i < 5; i++)
		{
			PO[i] = FR[i] + SR * URR[i] - SR * UR[i] - W * URR[i];
		}
	}

	P[1] = n1 * PO[1] + t1 * PO[2];     // Скорости
	P[2] = n2 * PO[1] + t2 * PO[2];
	P[0] = PO[0];                       // Масса
	P[3] = PO[3];                       // Энергия
	PQ = PO[4];

	return time;
}


void Setka::Init_Velosity(Sensor* sens, const double& A2, vector <double>& mu, vector <double>& Wt, vector <double>& Wp, vector <double>& Wr, const double& the)
{
	double Y = fabs(Velosity_inf);
	double Yr = -Y * cos(the); 
	double Yt = Y * sin(the); 
	double e1, e2, e3;

	e1 = sqrtpi_ * kv(Yr);
	e2 = 2.0 * fabs(Yr);
	e3 = 0.5 * sqrtpi_;

	double p1, p2, p4, p5;

	p1 = e1 / (e1 + e2 + e3);
	p2 = e2 / (e1 + e2 + e3);

	double ksi1, ksi2, ksi3, ksi4, ksi5, ksi6, ksi7, ksi8;
	double z = 0.0;
	double h = 0.0;
	double Wa = 0.0;

	for (int i = 1; i < I_; i++)
	{
		//cout << "Setka.cpp    " << 1 << endl;
		do
		{
			ksi1 = sens->MakeRandom();

			if (p1 > ksi1)
			{
				ksi2 = sens->MakeRandom();
				ksi3 = sens->MakeRandom();
				z = sqrt(-log(ksi2)) * cos(pi_ * ksi3);
			}
			else if (p1 + p2 > ksi1)
			{
				ksi2 = sens->MakeRandom();
				if (ksi2 <= 0.5)
				{
					z = -sqrt(-log(2.0 * ksi2));
				}
				else
				{
					z = sqrt(-log(2.0 * (1.0 - ksi2)));
				}
			}
			else
			{
				ksi2 = sens->MakeRandom();
				ksi3 = sens->MakeRandom();
				ksi4 = sens->MakeRandom();
				ksi5 = sens->MakeRandom();
				z = sign(ksi5 - 0.5) * sqrt(-log(ksi2) - log(ksi3) * kv(cos(pi_ * ksi4)));
			}

			Wr[i - 1] = z + Yr;
			h = kv(Yr + z) / kv(fabs(Yr) + fabs(z));
			ksi6 = sens->MakeRandom();
		} while (h <= ksi6 || z > -Yr);

		ksi7 = sens->MakeRandom();
		ksi8 = sens->MakeRandom();
		double V = 2.0 * pi_ * ksi7;
		if (i >= 2)
		{
			Wa = sqrt((gam(i - 2) + ksi8 * ((gam(i - 1) - gam(i - 2)))) * kv(Wr[i - 1]));
		}
		else
		{
			Wa = sqrt((ksi8 * ((gam(i - 1)))) * kv(Wr[i - 1]));
		}
		Wt[i - 1] = Wa * cos(V);
		Wp[i - 1] = Wa * sin(V);
		if (i >= 2)
		{
			mu[i - 1] = ((F_mk(gam(i - 1), Yr) - F_mk(gam(i - 2), Yr)) / A2) * fabs(Wr[i - 1]) * exp(-kv(Wt[i - 1] - Yt) - kv(Wp[i - 1]));
		}
		else
		{
			mu[i - 1] = ((F_mk(gam(i - 1), Yr)) / A2) * fabs(Wr[i - 1]) * exp(-kv(Wt[i - 1] - Yt) - kv(Wp[i - 1]));
		}
		//cout << "Setka.cpp    " << gam(i - 2) << endl;
	}



	p4 = (sqrtpi_ * fabs(Yr)) / (1.0 + sqrtpi_ * fabs(Yr));
	p5 = 1.0 - p4;
	double gamma_ = 0.0;
	if (I_ >= 2)
	{
		gamma_ = gam(I_ - 2);
	}

	do
	{
		//cout << "Setka.cpp    " << 2 << endl;
		do
		{
			//cout << "Setka.cpp    " << 3 << endl;
			ksi1 = sens->MakeRandom();
			ksi2 = sens->MakeRandom();
			ksi3 = sens->MakeRandom();
			if (p4 > ksi1)
			{
				z = sqrt(-log(ksi2)) * cos(pi_ * ksi3);
			}
			else
			{
				if (ksi2 <= 0.5)
				{
					z = -sqrt(-log(2.0 * ksi2));
				}
				else
				{
					z = sqrt(-log(2.0 * (1.0 - ksi2)));
				}
			}

			Wr[I_ - 1] = z + Yr;
			h = fabs(Yr + z) / (fabs(Yr) + fabs(z));
			ksi4 = sens->MakeRandom();
		} while (h <= ksi4 || z >= -Yr);

		ksi5 = sens->MakeRandom();
		ksi6 = sens->MakeRandom();
		//ksi7 = sens->MakeRandom();

		Wt[I_ - 1] = Yt + sqrt(-log(ksi5)) * cos(2.0 * pi_ * ksi6);
		Wp[I_ - 1] = sqrt(-log(ksi5)) * sin(2.0 * pi_ * ksi6);
	} while (kv(Wt[I_ - 1]) + kv(Wp[I_ - 1]) <= gamma_ * kv(Wr[I_ - 1]));

	mu[I_ - 1] = 1.0;
	for (int k = 0; k < I_ - 1; k++)
	{
		mu[I_ - 1] = mu[I_ - 1] - mu[k];
	}

	if (mu[I_ - 1] <= 0.0)
	{
		cout << "Setka.cpp    " << "6337 ERROR friuhuyvcyvawecsrygvirserccgr  " << mu[I_ - 1] << endl;
		exit(-1);
	}

}

void Setka::Velosity_initial2(Sensor* s, double& Vx, double& Vy, double& Vz)
{
	double ksi1 = s->MakeRandom();
	double ksi2 = s->MakeRandom();
	double a = sqrt(-log(1.0 - ksi2));
	Vy = a * cos(2.0 * pi_ * ksi1);
	Vz = a * sin(2.0 * pi_ * ksi1);
	//cout << "Setka.cpp    " << Vy << endl;
	double ksi3, ksi4, ksi5, ksi6;
	double z = 0;
	double p1 = 0.5 * fabs(Velosity_inf) * sqrtpi_ / (0.5 + 0.5 * fabs(Velosity_inf) * sqrtpi_);

	do
	{
		ksi3 = s->MakeRandom();
		ksi4 = s->MakeRandom();
		ksi5 = s->MakeRandom();
		ksi6 = s->MakeRandom();

		if (p1 > ksi3)
		{
			z = cos(pi_ * ksi5) * sqrt(-log(ksi4));
		}
		else
		{
			z = sqrt(-log(1.0 - ksi4));

		}
	} while (fabs(z + Velosity_inf) / (fabs(Velosity_inf) + fabs(z)) <= ksi6 || z < -Velosity_inf);

	Vx = z + Velosity_inf;
	if (Vx <= 0)
	{
		cout << "Setka.cpp    " << "dfEEERR 32424442" << endl;
	}
	return;
}

template<typename Random_type, typename Distribution_type, typename accuracy_type>
void Setka::Velosity_initial2(Random_type& gen, Distribution_type& dis, accuracy_type& Vx, accuracy_type& Vy, accuracy_type& Vz)
{
	double ksi1 = (dis(gen));
	double ksi2 = (dis(gen));
	double a = sqrt(-log(1.0 - ksi2));
	Vy = a * cos(2.0 * pi_ * ksi1);
	Vz = a * sin(2.0 * pi_ * ksi1);
	//cout << "Setka.cpp    " << Vy << endl;
	double ksi3, ksi4, ksi5, ksi6;
	double z = 0;
	double p1 = 0.5 * fabs(Velosity_inf) * sqrtpi_ / (0.5 + 0.5 * fabs(Velosity_inf) * sqrtpi_);

	do
	{
		ksi3 = (dis(gen));
		ksi4 = (dis(gen));
		ksi5 = (dis(gen));
		ksi6 = (dis(gen));

		if (p1 > ksi3)
		{
			z = cos(pi_ * ksi5) * sqrt(-log(ksi4));
		}
		else
		{
			z = sqrt(-log(1.0 - ksi4));

		}
	} while (fabs(z + Velosity_inf) / (fabs(Velosity_inf) + fabs(z)) <= ksi6 || z < -Velosity_inf);

	Vx = z + Velosity_inf;
	if (Vx <= 0)
	{
		cout << "Setka.cpp    " << "dfEEERR 32424442" << endl;
	}
	return;
}

void Setka::Velosity_initial(Sensor* s, double& Vx, double& Vy, double& Vz)
{
	double ksi1 = s->MakeRandom();
	double ksi2 = s->MakeRandom();
	double a = sqrt(-log(ksi2));
	Vy = a * cos(2.0 * pi_ * ksi1);
	Vz = a * sin(2.0 * pi_ * ksi1);
	//cout << "Setka.cpp    " << Vy << endl;
	double ksi3, ksi4, ksi5, ksi6;
	double z = 0;
	double p1 = fabs(Velosity_inf) * sqrtpi_ / (1.0 + fabs(Velosity_inf) * sqrtpi_);

	do
	{
		ksi3 = s->MakeRandom();
		ksi4 = s->MakeRandom();
		ksi5 = s->MakeRandom();
		ksi6 = s->MakeRandom();

		if (p1 > ksi3)
		{
			z = cos(pi_ * ksi5) * sqrt(-log(ksi4));
		}
		else
		{
			if (ksi4 <= 0.5)
			{
				z = -sqrt(-log(2.0 * ksi4));
			}
			else
			{
				z = sqrt(-log(2.0 * (1.0 - ksi4)));
			}
		}
	} while (fabs(z + Velosity_inf) / (fabs(Velosity_inf) + fabs(z)) < ksi6 || z > -Velosity_inf);

	Vx = z + Velosity_inf;
	return;
}

template<typename Random_type, typename Distribution_type, typename accuracy_type>
void Setka::Velosity_initial(Random_type& gen, Distribution_type& dis, accuracy_type& Vx, accuracy_type& Vy, accuracy_type& Vz)
{
	double ksi1 = (dis(gen));
	double ksi2 = (dis(gen));
	double a = sqrt(-log(ksi2));
	Vy = a * cos(2.0 * pi_ * ksi1);
	Vz = a * sin(2.0 * pi_ * ksi1);
	//cout << "Setka.cpp    " << Vy << endl;
	double ksi3, ksi4, ksi5, ksi6;
	double z = 0;
	double p1 = fabs(Velosity_inf) * sqrtpi_ / (1.0 + fabs(Velosity_inf) * sqrtpi_);

	do
	{
		ksi3 = (dis(gen));
		ksi4 = (dis(gen));
		ksi5 = (dis(gen));
		ksi6 = (dis(gen));

		if (p1 > ksi3)
		{
			z = cos(pi_ * ksi5) * sqrt(-log(ksi4));
		}
		else
		{
			if (ksi4 <= 0.5)
			{
				z = -sqrt(-log(2.0 * ksi4));
			}
			else
			{
				z = sqrt(-log(2.0 * (1.0 - ksi4)));
			}
		}
	} while (fabs(z + Velosity_inf) / (fabs(Velosity_inf) + fabs(z)) < ksi6 || z > -Velosity_inf);

	Vx = z + Velosity_inf;
	return;
}


double Setka::F_mk(const double& gamma, const double& Yr)   // Нужно задать
{
	return gamma * ((0.5 + kv(Yr)) * (1.0 + erf(-Yr)) - Yr * exp(-kv(Yr)) / sqrtpi_);
}

void Setka::Init_Pozision(Sensor* sens, const double& A1, double& phi, double& the)
{
	double ksi0, ksi1, ksi2, ksi3, ksi4, ksi5, ksi6;
	double X, h;
	double Y = fabs(Velosity_inf);

	ksi0 = sens->MakeRandom();
	ksi1 = sens->MakeRandom();
	phi = 2.0 * pi_ * ksi0;
	double p1 = erf(Y) / (A1 * kv(Y));

	if (p1 >= ksi1)
	{
		do
		{
			ksi2 = sens->MakeRandom();
			ksi3 = sens->MakeRandom();
			X = (1.0 / Y) * sqrt(-log(ksi2)) * cos(ksi3 * pi_ / 2.0);
		} while (X >= 1.0);
	}
	else
	{
		do
		{
			ksi2 = sens->MakeRandom();
			ksi3 = sens->MakeRandom();
			X = sqrt(ksi2);
			h = (1.0 + erf(X * Y)) / (1.0 + erf(Y));
		} while (h <= ksi3);
	}

	the = acos(X);
}

Cell* Setka::Find_cell(int& b, const double& x, const double& y)
{
	for (auto& i : this->All_Cells)
	{
		if (i->belong(x, y))
		{
			b = 1;
			return i;
		}
	}

	b = 0;
	return nullptr;
}
 

Cell* Setka::Belong_point(int b, const double& x, const double& y)
{
	if (b == 1)  // Внешняя сфера
	{
		for (auto& i : this->Cell_sphere)
		{
			if (i->belong(x, y))
			{
				return i;
			}
		}
	}
	else if (b == 2)  // Боковая граница
	{
		for (auto& i : this->Cell_side)
		{
			if (i->belong(x, y))
			{
				return i;
			}
		}
	}
	else if (b == 3)  // Задняя граница
	{
		for (auto& i : this->Cell_back)
		{
			if (i->belong(x, y))
			{
				return i;
			}
		}
	}
	else if (b == 4)  // Передний диск
	{
		for (auto& i : this->Cell_disk)
		{
			if (i->belong(x, y))
			{
				return i;
			}
		}
	}

	cout << "Setka.cpp    " << "ERRORRRORORfireubfvwcefrvjywgkvygdcwkug324h324334" << endl;
	cout << "Setka.cpp    " << "FFF  " << x << " " << y << " " << b << endl;
	//exit(-1);
	for (auto& i : this->Cell_side)
	{
		double x1, y1;
		i->Get_Center(x1, y1);
		cout << "Setka.cpp    " << "DDD  +  " << i->belong(x, y) << " " << x1 << " " << y1 << endl;
	}
	exit(-1);
	return nullptr;
}

void Setka::M_K_prepare(void)
{
	double xx1, yy1;
	double xx2, yy2;
	bool bb1 = false;
	bool bb2 = false;
	bool bb3 = false;


	// Заполним граничные ячейки
	//for (auto& i : this->All_Cells)    
	//{
	//	for (auto& j : i->Grans)
	//	{
	//		if (j->type == Input)
	//		{
	//			this->Cell_sphere.push_back(i);
	//			//cout << "Setka.cpp    " << j->A->x << " " << j->A->y << endl;
	//			break;
	//		}
	//	}
	//}


	for (auto& i : this->All_Cells)    // Заполняем массивы начальными гранями
	{
		for (int ijk = 0; ijk < n_S; ijk++)
		{
			i->S_m[ijk] = 0.0;
			i->S_p[ijk] = 0.0;
		}
		bb1 = false;
		bb2 = false;
		i->renew();
		for (auto& j : i->contour)
		{

			if (sqrt(kvv(j->x, j->y, 0.0)) <= (Rmax_ + 10.0 / RR_) && j->x > -10.1 / RR_)
			{
				bb1 = true;
			}

			if (sqrt(kvv(j->x, j->y, 0.0)) >= (Rmax_ - 10.0 / RR_) && j->x > -10.1 / RR_)
			{
				bb2 = true;
			}
		}

		if (bb1 == true && bb2 == true)
		{
			this->Cell_sphere.push_back(i);
			//double x1, y1;
			//i->Get_Center(x1, y1);
			//cout << "Setka.cpp    " << x1 << " " << y1 << endl;
		}

		if (bb2 == true)
		{
			this->Cell_other.push_back(i);
		}

		i->par[0].F_n = 0.0;
		i->par[0].F_u = 0.0;
		i->par[0].F_v = 0.0;
		i->par[0].F_T = 0.0;
		i->par[0].I_u = 0.0;
		i->par[0].I_v = 0.0;
		i->par[0].I_T = 0.0;
		i->par[0].II_u = 0.0;
		i->par[0].II_v = 0.0;
		i->par[0].II_T = 0.0;
		i->par[0].H_n[0] = 0.0;
		i->par[0].H_n[1] = 0.0;
		i->par[0].H_n[2] = 0.0;
		i->par[0].H_n[3] = 0.0;
		i->par[0].H_u[0] = 0.0;
		i->par[0].H_u[1] = 0.0;
		i->par[0].H_u[2] = 0.0;
		i->par[0].H_u[3] = 0.0;
		i->par[0].H_v[0] = 0.0;
		i->par[0].H_v[1] = 0.0;
		i->par[0].H_v[2] = 0.0;
		i->par[0].H_v[3] = 0.0;
		i->par[0].H_T[0] = 0.0;
		i->par[0].H_T[1] = 0.0;
		i->par[0].H_T[2] = 0.0;
		i->par[0].H_T[3] = 0.0;

		for (int aa = 0; aa < pop_atom; aa++)
		{
			i->par[0].H_n[aa] = 0.0;
			i->par[0].H_u[aa] = 0.0;
			i->par[0].H_v[aa] = 0.0;
			i->par[0].H_T[aa] = 0.0;
			i->par[0].H_uu[aa] = 0.0;
			i->par[0].H_uv[aa] = 0.0;
			i->par[0].H_vv[aa] = 0.0;
			i->par[0].H_uuu[aa] = 0.0;
		}

		for (int li = 0; li < 7; li++)
		{
			i->par[0].w_m[li] = 0.0;
		}

		double xxx, yyy;
		i->Get_Center(xxx, yyy);
		i->x_center = xxx;
		i->y_center = yyy;
		i->alf_center = polar_angle(xxx, yyy);
		i->y_ax = yyy;
		i->axis_ = false;
		for (auto& j : i->Grans)
		{
			if (j->type == Axis)
			{
				i->axis_ = true;
			}
		}

	}
	//exit(-1);

	cout << "Setka.cpp    " << "this->Cell_sphere = " << this->Cell_sphere.size() << endl;

	short int ikl = 0;
	this->Cell_side.clear();
	for (auto& i : this->All_Cells)    // Заполняем массивы начальными гранями
	{
		if (i->Grans.size() != 4)
		{
			cout << "ERROR   " << i->Grans.size() << endl;
			exit(-2);
		}

		double x1, y1;
		i->Get_Center(x1, y1);
		/*if (sqrt(kvv((x1 + 0.05), (y1 - 5.58), 0.0)) < 0.02)
		{
			cout << "Cell =  " << x1 << "  " << y1 << endl;
		}*/


		for (auto& j : i->Grans)
		{
			/*if (sqrt(kvv((x1 + 0.05), (y1 - 5.58), 0.0)) < 0.02)
			{
				cout << "type = " << j->type << endl;
				cout << x1 << "  " << y1 << endl;
			}*/

			if (j->type == Gran_type::Upper_wall)
			{
				this->Cell_side.push_back(i);
				//cout << "CC = " << x1 << "  " << y1 << endl;
				ikl++;
				break;
			}
		}
	}
	cout << "ikl = " << ikl << endl;

	for (auto& i : this->All_Cells)    // Заполняем массивы начальными ячейками
	{
		bb1 = false;
		bb2 = false;
		bb3 = true;

		for (auto& j : i->contour)
		{
			if (j->x <= 10.0 / RR_ && j->y > Rmax_ - 10.0 / RR_)
			{
				bb1 = true;
			}
			else if (j->x >= -10.0 / RR_ && j->y > Rmax_ - 10.0 / RR_)
			{
				bb2 = true;
			}

			if (j->x <= 3.0 / RR_ && j->x >= -3.0 / RR_ && j->y > Rmax_ - 10.0 / RR_)
			{
				bb3 = true;
			}
		}

		if ((bb1 == true && bb2 == true)||(bb3 == true))
		{
			this->Cell_disk.push_back(i);
		}
	}

	cout << "Setka.cpp    " << "this->Cell_disk = " << this->Cell_disk.size() << endl;

	for (auto& i : this->All_Cells)    // Заполняем массивы начальными ячейками
	{
		for (auto& j : i->Grans)
		{
			if (j->type == Extern)
			{
				this->Cell_back.push_back(i);
				break;
			}
		}
	}

	for (auto& i : this->All_Gran)    // Обновим уравнения граней для правильного нахождения пересечения траекторий с ними
	{
		i->renew();
	}

	for (auto& i : this->All_Gran_copy)
	{
		i->renew();
	}

	for (auto& i : this->All_Cells)
	{
		//if (i->type == C_centr || i->type == C_1 || i->type == C_2 || i->type == C_3)
		//{
		//	i->par[0].u = i->par[0].u * (chi_real / chi_);       // Перенормировка
		//	i->par[0].v = i->par[0].v * (chi_real / chi_);
		//	i->par[0].ro = i->par[0].ro / kv(chi_real / chi_);
		//	i->par[0].Q = i->par[0].Q / kv(chi_real / chi_);
		//}

		for (auto& j : i->Grans)
		{
			if (j->Sosed != nullptr)
			{
				i->L = min(i->L, j->Sosed->L);
			}
		}
	}


	// Блок загрузки датчиков случайных чисел
	ifstream fin2;
	fin2.open("rnd_my.txt");
	if (fin2.is_open() == false)
	{
		cout << "ERROR open  rnd_my.txt " << endl;
		exit(-100);
	}
	double d, a1, b1, c;
	for (int i = 0; i < 1021; i++)
	{
		fin2 >> a1 >> b1 >> c;
		auto s = new Sensor(a1, b1, c);
		this->Sensors.push_back(s);
	}
	fin2.close();

	fin2.open("Sensors276405.txt");
	unsigned int v1, v2, v3, v4;
	for (int i = 0; i < 276405; i++)
	{
		fin2 >> v1 >> v2 >> v3 >> v4;
		if (i % 1000 == 0)
		{
			auto s = new sensor2(v1, v2, v3, v4);
			this->Sensors2.push_back(s);
		}
	}
	fin2.close();


	double Y = fabs(Velosity_inf);
	this->sqv_1 = (Rmax_ - 3.0 / RR_) * (0.5 * (Rmax_ - 3.0 / RR_) * pi_ * Y * (erf(Y) * (1.0 + 1.0 / (2.0 * kv(Y))) + 1.0 + exp(-kv(Y)) / (Y * sqrtpi_)));
	this->sqv_2 = sqrtpi_ * (R5_ - 2.0 / RR_) * fabs(Left_ + 2.5 / RR_);
	this->sqv_3 = pi_ * kv(R5_ - 2.0 / RR_) * exp(-kv(Velosity_inf)) * (1.0 + exp(kv(Velosity_inf)) * sqrtpi_ * Velosity_inf * (1.0 + erf(Velosity_inf))) / (2.0 * sqrtpi_);
	this->sqv_4 = (sqrtpi_ / 2.0) * (kv((Rmax_ - 3.0 / RR_)) - kv(R5_ - 2.0 / RR_)) * exp(-kv(Velosity_inf)) * (sqrtpi_ * Velosity_inf * erfc(Velosity_inf) * exp(kv(Velosity_inf)) - 1.0);
	cout << "Setka.cpp    " << "this->sqv_1 = " << this->sqv_1 << endl;
	cout << "Setka.cpp    " << "this->sqv_4 = " << this->sqv_4 << endl;
	this->sum_s = this->sqv_1 + this->sqv_2 + this->sqv_3 + this->sqv_4;
	short int nmn = N_traektor; //10 * 3;  // 1 - это 3 минуты (7.15 минут)
	this->Number1 = 411 * 1000 * nmn;// * 154;// * 1071;// * 1440 * 2;//0 * 45;// *4000 * 2; // 411 * 20 * 353;// * 5 * 10 * 16;// * 100;// * 80;// * 20;// *36; // 280 * 250 * 7 * 20; // 280 * 62;// *23; // 280 * 120 * 2 * 3;// * 100; // 411 * 25 * 144; // 0 * 15;// *12 or 38;// *150;// *250; // * 10; // * 50; //250  6000;  250 * 50   411
	this->Number2 = 411 * 30 * nmn;// * 20;// *36; // 280 * 250 * 7; //280 * 62;// *23; // 280 * 10 * 2 * 3; //411 * 30; // * 30; // * 30; // 30; // 30;
	this->Number3 = 411 * 30 * nmn;// * 20; // 280 * 250; //280 * 62;// *23; // 280 * 3 * 2 * 3; //411 * 5; // * 5; // * 10; // 10;
	this->Number4 = 411 * 30 * nmn;// * 20;// *36; // 280 * 250 * 7; //280 * 62;// *23; // 280 * 80 * 2 * 3; //411 * 200; // * 200; // * 30; // 200; //  300  411 * 1650; // 135 * 40; // 30;
	this->AllNumber = ((this->Number1) + (this->Number2) + (this->Number3) + (this->Number4));
	cout << "Setka.cpp    " << "this->AllNumber " << this->AllNumber << endl;

	Ri.resize(I_ + 1);
	//Ri[0] = R5_ - 2.0;
	Ri[0] = 1.0 / RR_;                 // Здесь задаются радиусы
	Ri[1] = 2.61 / RR_;
	Ri[2] = 6.84 / RR_;
	Ri[3] = 17.92 / RR_;
	Ri[4] = 46.9 / RR_;
	Ri[5] = 122.7 / RR_;
	Ri[6] = 321.23 / RR_;                 // Здесь задаются радиусы
	Ri[7] = 840.65 / RR_;
	Ri[8] = (Rmax_ - 3.0 / RR_);

	for (int ik = 0; ik < J_; ++ik)
	{
		SINKR[ik] = max(fabs(sin((ik + 1) * pi_ / J_)), fabs(sin(ik * pi_ / J_)));
		//cout << "SINKR = " << SINKR[ik] << endl;
	}


	//Ri[4] = 540.0; // 216.0;
	//Ri[5] = Rmax_;
	//Ri[0] = 1.0;                 // Здесь задаются радиусы
	//Ri[1] = 2.0;
	//Ri[2] = 4.0;


	for (auto& i : this->All_Cells)    // Заполняем зоны ячейки
	{
		double x, y;
		i->Get_Center(x, y);
		double r = sqrt(kvv(x, y, 0.0));
		i->zona = geo_zones(r);
		i->zona_alpha = alpha_zones(x, y);
	}

	double mu1;
	mu1 = ((this->sqv_1) / this->sum_s) * (1.0 * this->AllNumber / this->Number1);  // Учёт, что для большого количества траекторий исходный вес меньше


	double kas = 10.0; // 6.0; // 20.0    0.5 для Kn = 3
	cout << "Setka.cpp    " << "kas = " << kas << endl;
	double ss1 = 10.0; // 1.0 * 0.5;
	double ss2 = 0.5; // 0.02; // * 0.5;    0.03
	double ss3 = 0.4; // 0.4; // 0.4;

	double ss4 = 1.0; // 0.005   0.01

	//double kas = 0.25 * mu1 * 100.0 * 10.0 * 0.25; // 0.00025; //0.001;
	//double ss1 = 1.0;
	//double ss2 = 1.0;
	//double ss3 = 1.0;

	// Основной сорт
	/*Mu[3][0] = 0.0001 * kas;
	Mu[3][1] = 0.001 * kas;
	Mu[3][2] = 0.01 * kas;
	Mu[3][3] = 0.1 * kas;
	Mu[3][4] = 1.0 * kas;
	Mu[3][5] = 10.0 * kas;
	Mu[3][6] = 100.0 * kas;
	Mu[3][7] = 1000.0 * kas;
	Mu[3][8] = 10000.0 * kas;*/

	Mu[3][0] = kv(Ri[0] / Rmax_) * kas * ss4;
	Mu[3][1] = kv(Ri[1] / Rmax_) * kas * ss4;
	Mu[3][2] = kv(Ri[2] / Rmax_) * kas * ss4;
	Mu[3][3] = kv(Ri[3] / Rmax_) * kas * ss4;
	Mu[3][4] = kv(Ri[4] / Rmax_) * kas * ss4; //0.2 * ss4;
	Mu[3][5] = kv(Ri[5] / Rmax_) * kas * ss4; //0.05 * ss4;
	Mu[3][6] = kv(Ri[6] / Rmax_) * kas * ss4; //0.13 * ss4;
	Mu[3][7] = kv(Ri[7] / Rmax_) * kas * ss4; //0.1 * ss4;
	Mu[3][8] = kv(Ri[8] / Rmax_) * kas * ss4; //0.1 * ss4;

	//cout << "All Weyghts = " << Mu[3][0] << " " << Mu[3][1] << " " << Mu[3][2] << " " << Mu[3][3] << //
	//	" " << Mu[3][4] << " " << Mu[3][5] << " " << Mu[3][6] << " " << Mu[3][7] << " " << Mu[3][8] << endl;

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < I_ + 1; j++)
		{
			this->Mu_stat[i][j] = 0.0;
			this->I_stat[i][j] = 0;
		}
	}

	/*Mu[3][0] = 0.0001 * kas;
	Mu[3][1] = 0.0001 * kas;
	Mu[3][2] = 0.0001 * kas;
	Mu[3][3] = 0.0001 * kas;
	Mu[3][4] = 0.0001 * kas;
	Mu[3][5] = 0.0001 * kas;
	Mu[3][6] = 0.0001 * kas;
	Mu[3][7] = 0.0001 * kas;
	Mu[3][8] = 0.0001 * kas;*/


	Mu[0][0] = Mu[3][0] * ss1;
	Mu[0][1] = Mu[3][1] * ss1;
	Mu[0][2] = Mu[3][2] * ss1;
	Mu[0][3] = Mu[3][3] * ss1;
	Mu[0][4] = Mu[3][4] * ss1;
	Mu[0][5] = Mu[3][5] * ss1;
	Mu[0][6] = Mu[3][6] * ss1;
	Mu[0][7] = Mu[3][7] * ss1;
	Mu[0][8] = Mu[3][8] * ss1;

	Mu[1][0] = Mu[3][0] * ss2;
	Mu[1][1] = Mu[3][1] * ss2;
	Mu[1][2] = Mu[3][2] * ss2;
	Mu[1][3] = Mu[3][3] * ss2;
	Mu[1][4] = Mu[3][4] * ss2;
	Mu[1][5] = Mu[3][5] * ss2;
	Mu[1][6] = Mu[3][6] * ss2 * 2.0;
	Mu[1][7] = Mu[3][7] * ss2 * 2.0;
	Mu[1][8] = Mu[3][8] * ss2 * 2.0;

	Mu[2][0] = Mu[3][0] * ss3;
	Mu[2][1] = Mu[3][1] * ss3;
	Mu[2][2] = Mu[3][2] * ss3;
	Mu[2][3] = Mu[3][3] * ss3;
	Mu[2][4] = Mu[3][4] * ss3;
	Mu[2][5] = Mu[3][5] * ss3 * 2.0;
	Mu[2][6] = Mu[3][6] * ss3 * 2.0;
	Mu[2][7] = Mu[3][7] * ss3 * 2.0;
	Mu[2][8] = Mu[3][8] * ss3 * 2.0;

	//return;

	fin2.open(parameter_6);
	if (fin2.is_open() == false)
	{
		cout << "ERROR open  stat_do.txt " << endl;
		exit(-100);
	}
	int x1;
	double cc;

	double koeff[4];
	koeff[0] = 0.5; // 0.1
	koeff[1] = 0.2; // 0.05
	koeff[2] = 0.5; // 0.5
	koeff[3] = 0.2; // 0.5

	double KKK = 5.0;


	if (true)
	{
		for (size_t k = 0; k < 4; k++)
		{
			for (size_t i = 0; i < I_; i++)
			{
				for (size_t j = 0; j < J_; j++)
				{
					fin2 >> x1 >> x1 >> x1 >> cc;
					Mu_[k][i][j] = max(cc * koeff[k] * KKK, 0.000001);  // 0.000001
				}
			}
		}
	}

	fin2.close();

	this->V_r_stat = new double[100000];
	this->V_t_stat = new double[100000];
	this->V_p_stat = new double[100000];
	this->mu_stat = new double[100000];
	this->phi_stat = new double[100000];
	this->num_stat = new int[100000];
	this->number_stat = 0;
	this->number_stat_2 = 0;
	//exit(-1);
}

int Setka::geo_zones(const double& r, const double& k)
// В какой геометрической зоне сейчас находится атом (зоны считаются с нуля)
{
	if (r < 0.0)
	{
		cout << "Setka.cpp    " << "Error  7151 geo_zones" << endl;
	}

	for (int i = 0; i < I_; i++)
	{
		if (r < k * this->Ri[i])
		{
			return i;
		}
	}
	return I_ - 1;
}

int Setka::alpha_zones(const double& x, const double& y)
// В какой геометрической зоне сейчас находится атом (зоны считаются с нуля)
{
	double al = polar_angle(x, y);

	for (int i = 0; i < J_; i++)
	{
		if (al < (i + 1) * pi_ / (J_ + 1))
		{
			return i;
		}
	}
	return J_ - 1;
}

double Setka::distination(const double& x0, const double& y0, const double& z0, 
	const double& Vx, const double& Vy, const double& Vz, const int& sort, int& to_iii, int& to_i, int& to_j)
{
	double time_do_peregel = (-Vx * x0 - Vy * y0 - Vz * z0) / kvv(Vx, Vy, Vz);
	double peregel = sqrt(kvv(x0 + Vx * time_do_peregel, y0 + Vy * time_do_peregel, z0 + Vz * time_do_peregel));
	int ii_z = this->geo_zones(peregel);                     // Номер зоны перегелия атома
	int  ii_alp  = alpha_zones(x0 + Vx * time_do_peregel, sqrt(kvv(y0 + Vy * time_do_peregel, z0 + Vz * time_do_peregel, 0.0)));
	
	to_iii = ii_z;
	to_i = ii_z;
	to_j = ii_alp;
	return sqrt(kvv(x0 + Vx * time_do_peregel, y0 + Vy * time_do_peregel, z0 + Vz * time_do_peregel));

	if (false)
	{
		time_do_peregel = time_to_peregel_axis(y0, z0, Vy, Vz);
		peregel = sqrt(kvv(x0 + Vx * time_do_peregel, y0 + Vy * time_do_peregel, z0 + Vz * time_do_peregel));
		int ii_z2 = this->geo_zones(peregel);                     // Номер зоны перегелия атома
		int ii_alp2 = alpha_zones(x0 + Vx * time_do_peregel, sqrt(kvv(y0 + Vy * time_do_peregel, z0 + Vz * time_do_peregel, 0.0)));
		to_iii = ii_z;

		if (Mu_[sort][ii_z][ii_alp] < Mu_[sort][ii_z2][ii_alp2])
		{
			to_i = ii_z;
			to_j = ii_alp;
			return 0.0; // !!!!!!!!!!!!!!!!!!!
		}
		else
		{
			to_i = ii_z2;
			to_j = ii_alp2;
			return 0.0;// !!!!!!!!!!!!!!!!!!!
		}
	}
}

void Setka::GD_prepare(void)
{
	for (auto& i : this->All_Cells)
	{
		if (i->pui_ == true)
		{
			i->par[0].ro = i->par[0].ro + i->par[0].npui;
			double dd = i->par[0].p;
			i->par[0].p = i->par[0].pp;

			i->par[0].pp = dd;

			//i->par[1] = i->par[0];
		}
	}
}

double Setka::get_w_init(const int& k)
{
	//return this->wmin * pow(this->wmax / this->wmin, 1.0 * k / (this->Nw - 1));  // Для степенной зависимости
	return this->wmin +  (this->wmax - this->wmin) * k / (this->Nw - 1);    // Для линейной зависимости
} 

void Setka::culc_PUI(void)
{
	ifstream fout;
	fout.open("fpui.txt");

	double wL, wR;
	int NN;

	fout >> wL >> wR >> NN;
	this->wmin = wL * 36.1275;
	this->wmax = wR * 36.1275;
	this->Nw = NN;

	cout << this->wmin << "  min - max  " << this->wmax << endl;

	double a;
	for (int i = 0; i < NN; i++)
	{
		fout >> a;
	}

	double nd = 5.36502e-20;

	for (auto& i : this->All_Cells)
	{
		double x, y, d, Tsw;
		int ii, iii;

		fout >> ii >> x >> y >> iii >> i->par[0].npui >> i->par[0].Tpui >> d >> d >> d >> Tsw;
		i->par[0].npui = i->par[0].npui / 0.06;
		i->par[0].Tpui = i->par[0].Tpui / 6527.0;
		i->par[0].ppui = i->par[0].npui * i->par[0].Tpui;
		Tsw = Tsw / 6527.0;

		i->par[0].pp = i->par[0].p;

		i->fpui_max = 0.0;
		for (int k = 0; k < NN; k++)
		{
			fout >> a;
			i->fpui.push_back(a / nd);
			i->fpui_max = max(i->fpui_max, a / nd);
		}

		i->Wmin = this->wmin;
		//for (int k = 0; k < NN; k++)
		//{
		//	if (i->fpui[k] > i->fpui_max / 1000.0)
		//	{
		//		//i->Wmin = this->wmin * pow(this->wmax / this->wmin, 1.0 * k / (this->Nw - 1));
		//		i->Wmin = this->get_w_init(k);
		//		break;
		//	}
		//}

		i->Wmax = this->wmax;
		//for (int k = NN - 1; k >= 0; k--)
		//{
		//	if (i->fpui[k] > i->fpui_max / 10000.0)
		//	{
		//		//i->Wmax = this->wmin * pow(this->wmax / this->wmin, 1.0 * k / (this->Nw - 1));
		//		i->Wmax = this->get_w_init(k);
		//		break;
		//	}
		//}

		double SSS = 0.0;
		for (int ik = 0; ik < this->Nw - 1; ik++)
		{
			double L = this->get_w_init(ik);
			double R = this->get_w_init(ik + 1);
			double f1 = i->fpui[ik] * L * L * L * L;
			double f2 = i->fpui[ik + 1] * R * R * R * R;
			SSS = SSS + 4.0 * pi_ * 0.5 * (f1 + f2) * (R - L);
		}

		double SS = 0.0;
		for (int ik = 0; ik < this->Nw - 1; ik++)
		{
			double L = this->get_w_init(ik);
			double R = this->get_w_init(ik + 1);
			double f1 = i->fpui[ik] * L * L * L * L;
			double f2 = i->fpui[ik + 1] * R * R * R * R;
			SS = SS + 4.0 * pi_ * 0.5 * (f1 + f2) * (R - L);
			if (SS > SSS * 0.999)
			{
				i->Wmax = R;
				break;
			}
		}
		//cout << i->par[0].npui << "  " << SS << "  " << i->number << endl;
		//cout << i->contour[0]->x << "  " << i->contour[0]->y << "  " << i->Wmin << "   " << i->Wmax << endl;

		if (i->par[0].npui > 0.001 * i->par[0].ro)
		{
			if (i->par[0].npui > i->par[0].ro)
			{
				cout << "pui > ro   " << i->contour[0]->x << "  " << i->contour[0]->y << "  " << i->par[0].npui  << "  " << i->par[0].ro << endl;
				double d = (1.01 * i->par[0].npui) / i->par[0].ro;
				i->par[0].npui = i->par[0].npui / d;
				i->par[0].ppui = i->par[0].ppui / d;
				for (auto& j : i->fpui)
				{
					j = j / d;
				}
			}
			i->par[0].ro = i->par[0].ro - i->par[0].npui;
			i->par[0].p = 2.0 * i->par[0].ro * Tsw;
			i->pui_ = true;

			if (i->par[0].p <= 0.0)
			{
				double rr = 1.0 / RR_;
				double xx, yy;
				i->Get_Center(xx, yy);
				double dist = sqrt(kv(xx) + kv(yy));
				i->par[0].p = 1109.31 * pow(rr / dist, 2.0 * ggg);
			}

		}

		//i->par[1] = i->par[0];


		/*if (i->pui_ == true)
		{
			cout << this->wmin << "  " << this->wmax << endl;
			cout << i->fpui[0] << "  " << i->fpui[1] << "  " << i->fpui[2] << endl;
			cout << i->get_fpui(10.0, this->wmin, this->wmax) << endl;
			exit(-1);
		}*/


	}


	// Надо разделить функцию распредления на области
	for (auto& i : this->All_Cells)
	{
		double mi = i->fpui[0];
		double ma = i->fpui[0];
		i->Wmin_sort.push_back(i->Wmin);
		bool b1 = false;
		for (int k = 0; k < i->fpui.size(); k++)
		{
			double L = this->get_w_init(k);
			if (L < i->Wmin && b1 == false)
			{
				b1 = true;
				continue;
			}
			if (mi > i->fpui[k])
			{
				mi = i->fpui[k];
			}
			if (ma < i->fpui[k])
			{
				ma = i->fpui[k];
			}

			if (mi * 100 / ma < 3.0 || L >= i->Wmax)
			{
				i->fpui_max_sort.push_back(ma);
				i->Wmax_sort.push_back(L);
				mi = i->fpui[k];
				ma = i->fpui[k];
				if (L >= i->Wmax)
				{
					break;
				}
				i->Wmin_sort.push_back(L);
			}
		}

		i->i_pui = i->Wmin_sort.size();
	}

	ofstream fin;
	fin.open("PUIIII_1130.txt");

	for (int i = 0; i < this->Nw; i++)
	{
		//double L = this->wmin * pow(this->wmax / this->wmin, 1.0 * i / (this->Nw - 1));
		double L = this->get_w_init(i);
		fin << L << " " << this->All_Cells[1130]->fpui[i] << endl;
	}

	for (int ik = 0; ik < this->All_Cells[1130]->Wmin_sort.size(); ik++)
	{
		cout << this->All_Cells[1130]->Wmin_sort[ik] << " " << this->All_Cells[1130]->Wmax_sort[ik] << " " << this->All_Cells[1130]->fpui_max_sort[ik] << endl;
	}
	cout << endl;
	fin.close();

	fin.open("PUIIII_23.txt");

	for (int i = 0; i < this->Nw; i++)
	{
		//double L = this->wmin * pow(this->wmax / this->wmin, 1.0 * i / (this->Nw - 1));
		double L = this->get_w_init(i);
		fin << L << " " << this->All_Cells[23]->fpui[i] << endl;
	}
	for (int ik = 0; ik < this->All_Cells[23]->Wmin_sort.size(); ik++)
	{
		cout << this->All_Cells[23]->Wmin_sort[ik] << " " << this->All_Cells[23]->Wmax_sort[ik] << " " << this->All_Cells[23]->fpui_max_sort[ik] << endl;
	}
	cout << endl;
	fin.close();

	//cout << this->All_Cells[1130]->get_fpui(38.0, this->wmin, this->wmax) << " " << endl;;
	//cout << this->All_Cells[1130]->Wmin << " " << this->All_Cells[1130]->Wmax << endl;
	//exit(-1);

	// Считаем всякие полезные интеграллы для последующих перезарядок
	if (false)
	{



//#pragma omp parallel for
			for (auto& i : this->All_Cells)
			{
				double L, R;
				double f1, f2;
				double ff1, ff2;
				double fff1, fff2;
				double d;
				int n1 = 100;
				double dthe = pi_ / (n1 - 1);

				if (i->pui_ == true)
				{
					i->nu_pui.resize(i->i_pui);
					i->nu2_pui.resize(i->i_pui);
					i->nu3_pui.resize(i->i_pui);
					if (i->number % 55 == 0)
					{
						cout << i->number << endl;
					}

					vector <double> pui_I1;
					vector <double> pui_I2;
					vector <double> pui_I3;

					for (int k = 0; k < k_Igor; k++)
					{
						double LL = L_Igor / (k_Igor - 1) * k;
						double Int1 = 0.0;
						double Int2 = 0.0;
						double Int3 = 0.0;
						int kj = 0;
						for (int ik = 0; ik < this->Nw - 1; ik++)
						{
							L = this->get_w_init(ik);
							R = this->get_w_init(ik + 1);
							if (R > i->Wmax)
							{
								break;
							}
							for (int ij = 0; ij < n1; ij++)
							{
								double the = dthe * ij;
								d = sqrt(L * L + LL * LL - 2.0 * L * LL * cos(the));
								f1 = i->fpui[ik] * d * sigma(d) * L * L * sin(the);
								ff1 = i->fpui[ik] * d * sigma(d) * L * L * sin(the) * (LL - L * cos(the));
								fff1 = i->fpui[ik] * d * d * d * sigma(d) * L * L * sin(the);
								d = sqrt(R * R + LL * LL - 2.0 * R * LL * cos(the));
								f2 = i->fpui[ik + 1] * d * sigma(d) * R * R * sin(the);
								ff2 = i->fpui[ik + 1] * d * sigma(d) * R * R * sin(the) * (LL - R * cos(the));
								fff2 = i->fpui[ik + 1] * d * d * d * sigma(d) * R * R * sin(the);
								Int1 = Int1 + 0.5 * (R - L) * (f2 + f1) * dthe;
								Int2 = Int2 + 0.5 * (R - L) * (ff2 + ff1) * dthe;
								Int3 = Int3 + 0.5 * (R - L) * (fff2 + fff1) * dthe;
							}
							if (R >= i->Wmax_sort[kj])
							{
								i->nu_pui[kj].push_back(Int1 * 2.0 * pi_);
								i->nu2_pui[kj].push_back(Int2 * 2.0 * pi_);
								i->nu3_pui[kj].push_back(Int3 * 2.0 * pi_);
								Int1 = 0.0;
								Int2 = 0.0;
								Int3 = 0.0;
								kj++;
							}
						}
					}
					//cout << IM << endl;
				}
			}

			/*for (int ik = 0; ik < this->All_Cells[23]->i_pui; ik++)
			{
				cout << this->All_Cells[23]->nu_pui[ik][19] << " " << this->All_Cells[23]->nu2_pui[ik][19] << " " << this->All_Cells[23]->nu3_pui[ik][19] << endl;
			}*/

		
	}
}

void Setka::proverka_Istok(int ni)
{
	auto i = this->All_Cells[ni];
	double L, R;
	double f1, f2;
	double ff1, ff2;
	double fff1, fff2;
	double d;
	int n1 = 1000;
	double dthe = pi_ / (n1 - 1);

	vector <double> pui_I1;
	vector <double> pui_I2;
	vector <double> pui_I3;

	double Vh1 = 1.0;
	double Vh2 = -2.0;
	double Vh3 = 3.0;
	double Vp1 = 1.0;
	double Vp2 = 3.0;
	double Vp3 = 1.0;

	double Lx = Vh1 - Vp1;
	double Ly = Vh2 - Vp2;
	double Lz = Vh3 - Vp3;
	double LL = sqrt(kv(Lx) + kv(Ly) + kv(Lz));


	double Int1 = 0.0;
	double Int2 = 0.0;
	double Int3 = 0.0;
	int kj = 0;
	for (int ik = 0; ik < this->Nw - 1; ik++)
	{
		L = this->get_w_init(ik);
		R = this->get_w_init(ik + 1);

		for (int ij = 0; ij < n1; ij++)
		{
			double the = dthe * ij;
			d = sqrt(L * L + LL * LL - 2.0 * L * LL * cos(the));
			f1 = i->fpui[ik] * d * sigma(d) * L * L * sin(the);
			ff1 = i->fpui[ik] * d * sigma(d) * L * L * sin(the) * (LL - L * cos(the));
			fff1 = i->fpui[ik] * d * d * d * sigma(d) * L * L * sin(the);
			d = sqrt(R * R + LL * LL - 2.0 * R * LL * cos(the));
			f2 = i->fpui[ik + 1] * d * sigma(d) * R * R * sin(the);
			ff2 = i->fpui[ik + 1] * d * sigma(d) * R * R * sin(the) * (LL - R * cos(the));
			fff2 = i->fpui[ik + 1] * d * d * d * sigma(d) * R * R * sin(the);
			Int1 = Int1 + 0.5 * (R - L) * (f2 + f1) * dthe;
			Int2 = Int2 + 0.5 * (R - L) * (ff2 + ff1) * dthe;
			Int3 = Int3 + 0.5 * (R - L) * (fff2 + fff1) * dthe;
		}
	}

	double skalar2 = (Vh1 * Lx + Vh2 * Ly + Vh3 * Lz) / LL;
	Int3 = (-0.5 * Int3 + Int2 * skalar2);
	
	cout << Int1 * 2.0 * pi_ << "  " << Int2 * 2.0 * pi_ << "  " << Int3 * 2.0 * pi_ << endl << endl;
	Int1 = 0.0;
	Int2 = 0.0;
	Int3 = 0.0;

	double Int2x = 0.0;
	double Int2y = 0.0;
	double Int2z = 0.0;

	double dx = 0.05;
	double dy = 0.05;
	double dz = 0.05;

	for (double x = -50.0; x < 50.0; x = x + dx)
	{
		for (double y = -50.0; y < 50.0; y = y + dy)
		{
			for (double z = -50.0; z < 50.0; z = z + dz)
			{
				double d1 = sqrt(kv(Vh1 - x) + kv(Vh2 - y) + kv(Vh3 - z));
				double d2 = sqrt(kv(Vp1 - x) + kv(Vp2 - y) + kv(Vp3 - z));
				Int1 = Int1 + d1 * sigma(d1) * dx * dy * dz * i->get_fpui(d2, this->wmin, this->wmax);
				Int2x = Int2x + d1 * sigma(d1) * dx * dy * dz * i->get_fpui(d2, this->wmin, this->wmax) * (Vh1 - x);
				Int2y = Int2y + d1 * sigma(d1) * dx * dy * dz * i->get_fpui(d2, this->wmin, this->wmax) * (Vh2 - y);
				Int2z = Int2z + d1 * sigma(d1) * dx * dy * dz * i->get_fpui(d2, this->wmin, this->wmax) * (Vh3 - z);
				Int3 = Int3 + d1 * sigma(d1) * dx * dy * dz * i->get_fpui(d2, this->wmin, this->wmax) * 0.5 * (kvv(Vh1, Vh2, Vh3) - kvv(x, y, z));
			}
		}
	}



	Int2 = sqrt(kv(Int2x) + kv(Int2y) + kv(Int2z));

	cout << Int1 << "  " << Int2 << "  " << Int3 << endl << endl << endl;
	
}

void Setka::MK_start(void)
{

	mutex mut_1;
	double Y = fabs(Velosity_inf);
	int st = 1;

#pragma omp parallel for
	for (int index = 0; index < 135; index++)
	{
		double A1 = 1.0 + (1.0 + 1.0 / (2.0 * kv(Y))) * erf(Y) + exp(-kv(Y)) / (Y * sqrtpi_);
		double A2 = 0.0; // Нужно задать после нахождения угла theta
		double mu1, mu2, mu3, mu4;
		double Vx, Vy, Vz;
		mu1 = ((this->sqv_1) / this->sum_s) * (1.0 * this->AllNumber / this->Number1);
		mu2 = ((this->sqv_2) / this->sum_s) * (1.0 * this->AllNumber / this->Number2);
		mu3 = ((this->sqv_3) / this->sum_s) * (1.0 * this->AllNumber / this->Number3);
		mu4 = ((this->sqv_4) / this->sum_s) * (1.0 * this->AllNumber / this->Number4);
		if (index == 0)
		{
			mut_1.lock();
			cout << "Setka.cpp    " << "Wright weyght   " << mu1 << " " << mu2 << " " << mu3 << " " << mu4 << endl;
			mut_1.unlock();
		}
		Sensor* sens1 = Sensors[2 * index];
		Sensor* sens2 = Sensors[2 * index + 1];
		MKmethod MK = MKmethod();
		mut_1.lock();
		cout << "Setka.cpp    " << st << " potok  is  135;  index = " << index << endl;
		st++;
		mut_1.unlock();
		double x, y, z;
		double ksi1, ksi2, ksi3, ksi4, ksi5, ksi6, ksi7;
		double phi, Vphi, Vr;
		double a, b, c;
		double r;

		for (int ii = 0; ii < Number1 / 135; ii++)  //  Вылет с полусферы
		{
			/*if (ii % 100000 == 0 && ii > 1 )
			{
				cout << "Setka.cpp    " << "ii =  " << ii << endl;
			}*/
			double phi, the;
			// Нахождение начальной траектории
			this->Init_Pozision(sens1, A1, phi, the);
			// перевод позиции в декартовы координаты
			x = (R5_ - 2.0) * cos(the);
			y = (R5_ - 2.0) * sin(the) * cos(phi);
			z = (R5_ - 2.0) * sin(the) * sin(phi);
			//cout << "Setka.cpp    " << x << " " << sqrt(kv(y)+ kv(z)) << endl;

			Cell* Point = Belong_point(1, x, sqrt(kv(z) + kv(y)));  // Находит ячейку, которой принадлежит точка
			vector <double> mu(I_);
			vector <double> Wt(I_);
			vector <double> Wp(I_);
			vector <double> Wr(I_);

			double X = cos(the);
			A2 = exp(-kv(Y) * kv(X)) / sqrtpi_ + Y * X * (1.0 + erf(Y * X));
			Init_Velosity(sens1, A2, mu, Wt, Wp, Wr, the);

			//Velosity_initial(sens1, Vx, Vy, Vz);
			//Fly_exchenge_Imit(sens2, x, y, z, Vx, Vy, Vz, Point, mu1, -log(1.0 - sens1->MakeRandom()), 0.0, 3, mu1, I_ - 1, 5);

			//cout << "Setka.cpp    " << "Start" << endl;
			if (true)
			{
				for (int i = 0; i < I_; i++)
				{
					dekard_skorost2((R5_ - 2.0), the, phi, Wr[i], Wp[i], Wt[i], Vy, Vz, Vx);
					/*for (int k = 0; k < 60000; k++)
					{
						cout << "Setka.cpp    " << x + 0.3 * k * Vx << " " << sqrt(kv(y + 0.3 * k * Vy) + kv(z + 0.3 * k * Vz)) << endl;
					}*/
					//cout << "Setka.cpp    " << "START" << endl;
					//Fly_exchenge(sens2, x, y, z, Vx, Vy, Vz, Point, mu[i] * mu1, mu[i] * mu1, false);
					//Fly_exchenge_Split(sens2, x, y, z, Vx, Vy, Vz, Point, mu[i] * mu1, mu[i] * mu1, false, i);
					//cout << "Setka.cpp    " << "A" << endl;

					//Fly_exchenge_Imit_Korol(sens2, x, y, z, Vx, Vy, Vz, Point, mu[i] * mu1,  3, false, mu1, I_ - 1);
					//cout << "Setka.cpp    " << ii << endl;
					Fly_exchenge_Imit(MK, sens2, x, y, z, Vx, Vy, Vz, Point, mu[i] * mu1, -log(1.0 - sens1->MakeRandom()), 0.0, 3, mu[i] * mu1, I_ - 1, ii);

				}
			}
		}
		for (int ii = 0; ii < Number2 / 135; ii++)  //  Вылет сверху 
		{
			//cout << "Setka.cpp    " << "A" << endl;
			ksi1 = sens1->MakeRandom();
			ksi2 = sens1->MakeRandom();
			ksi3 = sens1->MakeRandom();
			ksi4 = sens1->MakeRandom();
			ksi5 = sens1->MakeRandom();

			x = (Left_ + 0.1) + ksi1 * (- 0.2 - Left_);
			phi = ksi2 * 2.0 * pi_;
			Vphi = cos(2.0 * pi_ * ksi3) * sqrt(-log(1.0 - ksi4));
			Vx = Velosity_inf + sin(2.0 * pi_ * ksi3) * sqrt(-log(1.0 - ksi4));
			Vr = -sqrt(-log(ksi5));
			y = (R5_ - 0.2) * cos(phi);
			z = (R5_ - 0.2) * sin(phi);

			Cell* Point = Belong_point(2, x, sqrt(kv(z) + kv(y)));  // Находит ячейку, которой принадлежит точка
			//Fly_exchenge_Imit_Korol(sens2, x, y, z, Vx, cos(phi) * Vr - sin(phi) * Vphi, sin(phi) * Vr + cos(phi) * Vphi,//
			//	Point, mu2, 3, false, mu2, I_ - 1);
			Fly_exchenge_Imit(MK, sens2, x, y, z, Vx, cos(phi) * Vr - sin(phi) * Vphi, sin(phi) * Vr + cos(phi) * Vphi,//
				Point, mu2, -log(1.0 - sens1->MakeRandom()), 0.0, 3, mu2, I_ - 1, ii);
		}
		for (int ii = 0; ii < Number3 / 135; ii++)  //  Вылет сзади
		{
			//cout << "Setka.cpp    " << "B" << endl;
			Velosity_initial2(sens1, a, b, c);
			ksi1 = sens1->MakeRandom();
			ksi2 = sens1->MakeRandom();

			r = sqrt(ksi1 * (R5_ - 2.0) * (R5_ - 2.0));
			phi = ksi2 * 2.0 * pi_;
			y = r * cos(phi);
			z = r * sin(phi);

			Cell* Point = Belong_point(3, Left_ + 2.0, sqrt(kv(z) + kv(y)));  // Находит ячейку, которой принадлежит точка
			//cout << "Setka.cpp    " << "Start" << endl;
			//Fly_exchenge_Imit_Korol(sens2, Left_ + 2.0, y, z, a, b, c, Point, mu3, 3, false, mu3, I_ - 1);
			Fly_exchenge_Imit(MK, sens2, Left_ + 2.0, y, z, a, b, c, Point, mu3, -log(1.0 - sens1->MakeRandom()), 0.0, //
				3, mu3, I_ - 1, ii);
			//cout << "Setka.cpp    " << "Stop" << endl;
		}
		for (int ii = 0; ii < Number4 / 135; ii++)  //Number4 / 135
		{
			double a, b, c;
			Velosity_initial(sens1, a, b, c);
			ksi1 = sens1->MakeRandom();
			ksi2 = sens1->MakeRandom();

			r = sqrt(-(ksi1 * (kv(Rmax_) - kv(R5_)) - kv(Rmax_)));
			phi = ksi2 * 2.0 * pi_;
			y = r * cos(phi);
			z = r * sin(phi);

			Cell* Point = Belong_point(4, -0.1, sqrt(kv(z) + kv(y)));

			Fly_exchenge_Imit(MK, sens2, -0.1, y, z, a, b, c, Point, mu4, -log(1.0 - sens1->MakeRandom()), 0.0, //
				3, mu4, I_ - 1, ii);
		}
	}


	for (auto& k : this->All_Cells)
	{

		/*for (int il = 0; il < 7; il++)
		{
			if (k->par[0].num_atoms > 0)
			{
				k->par[0].w_m[il] /= k->par[0].num_atoms;
			}
		}*/


		double no = (1.0 * AllNumber * k->Get_Volume_rotate(360.0));

		k->par[0].I_u = sum_s * k->par[0].I_u / no;
		k->par[0].I_v = sum_s * k->par[0].I_v / no;
		k->par[0].I_T = sum_s * k->par[0].I_T / no;

		k->par[0].II_u = sum_s * k->par[0].II_u / no;
		k->par[0].II_v = sum_s * k->par[0].II_v / no;
		k->par[0].II_T = sum_s * k->par[0].II_T / no;

		if (k->par[0].F_n > 0)
		{
			/*k->par[0].I_u = k->par[0].I_u / k->par[0].F_n;
			k->par[0].I_v = k->par[0].I_v / k->par[0].F_n;
			k->par[0].I_T = k->par[0].I_T / k->par[0].F_n;*/

			k->par[0].F_u = k->par[0].F_u / k->par[0].F_n;
			k->par[0].F_v = k->par[0].F_v / k->par[0].F_n;
			k->par[0].F_T = (2.0 / 3.0) * (k->par[0].F_T / k->par[0].F_n - kvv(k->par[0].F_u, k->par[0].F_v, 0.0));
		}
		else
		{
			/*k->par[0].I_u = 0.0;
			k->par[0].I_v = 0.0;
			k->par[0].I_T = 0.0;*/

			k->par[0].F_n = 0.0;
			k->par[0].F_u = 0.0;
			k->par[0].F_v = 0.0;
			k->par[0].F_T = 0.0;

		}

		if (k->par[0].F_n > 0)
		{
			for (int i = 0; i < 4; i++)
			{
				if (k->par[0].H_n[i] > 0)
				{
					k->par[0].H_u[i] = k->par[0].H_u[i] / k->par[0].H_n[i];
					k->par[0].H_v[i] = k->par[0].H_v[i] / k->par[0].H_n[i];
					//k->par[0].H_T[i] = (2.0 / 3.0) * (k->par[0].H_T[i] / k->par[0].H_n[i] - kvv(k->par[0].F_u, k->par[0].F_v, 0.0));
					k->par[0].H_T[i] = (2.0 / 3.0) * (k->par[0].H_T[i] / k->par[0].H_n[i] - kvv(k->par[0].H_u[i], k->par[0].H_v[i], 0.0));
				}
				else
				{
					k->par[0].H_u[i] = 0.0;
					k->par[0].H_v[i] = 0.0;
					k->par[0].H_T[i] = 0.0;
				}
			}
		}

		k->par[0].F_n = sum_s * k->par[0].F_n / no;
		for (int i = 0; i < 4; i++)
		{
			k->par[0].H_n[i] = sum_s * k->par[0].H_n[i] / no;
		}

		/*k->par[0].I_u = (n_p_LISM_ / Kn_) * (-k->par[0].ro * k->par[0].I_u);
		k->par[0].I_v = (n_p_LISM_ / Kn_) * (-k->par[0].ro * k->par[0].I_v);
		k->par[0].I_T = (n_p_LISM_ / Kn_) * (-(1.0 / 2.0) * k->par[0].ro * k->par[0].I_T);*/

		/*k->par[0].I_u = (n_p_LISM_ / Kn_) * (k->par[0].I_u);
		k->par[0].I_v = (n_p_LISM_ / Kn_) * (k->par[0].I_v);
		k->par[0].I_T = (n_p_LISM_ / Kn_) * (k->par[0].I_T);*/
		k->par[0].I_u = (n_p_LISM_) * (k->par[0].I_u);
		k->par[0].I_v = (n_p_LISM_) * (k->par[0].I_v);
		k->par[0].I_T = (n_p_LISM_) * (k->par[0].I_T);
		k->par[0].II_u = (n_p_LISM_) * (k->par[0].II_u);
		k->par[0].II_v = (n_p_LISM_) * (k->par[0].II_v);
		k->par[0].II_T = (n_p_LISM_) * (k->par[0].II_T);

		// Считаем мультифдюидные источники

		double U_M_H1, U_M_H2, U_M_H3, U_M_H4;
		double U_H1, U_H2, U_H3, U_H4;
		double sigma_H1, sigma_H2, sigma_H3, sigma_H4;
		double nu_H1, nu_H2, nu_H3, nu_H4;
		double q2_1, q2_2, q3;
		double u, v, ro, p, Q;

		double u_H1 = k->par[0].H_u[0], v_H1 = k->par[0].H_v[0], ro_H1 = k->par[0].H_n[0], p_H1 = 0.5 * k->par[0].H_T[0] * k->par[0].H_n[0];
		double u_H2 = k->par[0].H_u[1], v_H2 = k->par[0].H_v[1], ro_H2 = k->par[0].H_n[1], p_H2 = 0.5 * k->par[0].H_T[1] * k->par[0].H_n[1];
		double u_H3 = k->par[0].H_u[2], v_H3 = k->par[0].H_v[2], ro_H3 = k->par[0].H_n[2], p_H3 = 0.5 * k->par[0].H_T[2] * k->par[0].H_n[2];
		double u_H4 = k->par[0].H_u[3], v_H4 = k->par[0].H_v[3], ro_H4 = k->par[0].H_n[3], p_H4 = 0.5 * k->par[0].H_T[3] * k->par[0].H_n[3];
		u = k->par[0].u;
		v = k->par[0].v;
		ro = k->par[0].ro;
		p = k->par[0].p;

		if (ro <= 0.0)
		{
			ro = 0.0000001;
			p = 0.0;
		}
		if (ro_H1 <= 0.0)
		{
			ro_H1 = 0.0000001;
			p_H1 = 0.0;
		}
		if (ro_H2 <= 0.0)
		{
			ro_H2 = 0.0000001;
			p_H2 = 0.0;
		}
		if (ro_H3 <= 0.0)
		{
			ro_H3 = 0.0000001;
			p_H3 = 0.0;
		}
		if (ro_H4 <= 0.0)
		{
			ro_H4 = 0.0000001;
			p_H4 = 0.0;
		}

		U_M_H1 = sqrt(kv(u - u_H1) + kv(v - v_H1) + (64.0 / (9.0 * pi_)) //
			* (p / ro + 2.0 * p_H1 / ro_H1));
		U_M_H2 = sqrt(kv(u - u_H2) + kv(v - v_H2) + (64.0 / (9.0 * pi_)) //
			* (p / ro + 2.0 * p_H2 / ro_H2));
		U_M_H3 = sqrt(kv(u - u_H3) + kv(v - v_H3) + (64.0 / (9.0 * pi_)) //
			* (p / ro + 2.0 * p_H3 / ro_H3));
		U_M_H4 = sqrt(kv(u - u_H4) + kv(v - v_H4) + (64.0 / (9.0 * pi_)) //
			* (p / ro + 2.0 * p_H4 / ro_H4));

		U_H1 = sqrt(kv(u - u_H1) + kv(v - v_H1) + (4.0 / pi_) //
			* (p / ro + 2.0 * p_H1 / ro_H1));
		U_H2 = sqrt(kv(u - u_H2) + kv(v - v_H2) + (4.0 / pi_) //
			* (p / ro + 2.0 * p_H2 / ro_H2));
		U_H3 = sqrt(kv(u - u_H3) + kv(v - v_H3) + (4.0 / pi_) //
			* (p / ro + 2.0 * p_H3 / ro_H3));
		U_H4 = sqrt(kv(u - u_H4) + kv(v - v_H4) + (4.0 / pi_) //
			* (p / ro + 2.0 * p_H4 / ro_H4));

		sigma_H1 = kv(1.0 - a_2 * log(U_M_H1)); // 0.1243
		sigma_H2 = kv(1.0 - a_2 * log(U_M_H2));
		sigma_H3 = kv(1.0 - a_2 * log(U_M_H3)); // 0.1121     a_2
		sigma_H4 = kv(1.0 - a_2 * log(U_M_H4));

		nu_H1 = ro * ro_H1 * U_M_H1 * sigma_H1;
		nu_H2 = ro * ro_H2 * U_M_H2 * sigma_H2;
		nu_H3 = ro * ro_H3 * U_M_H3 * sigma_H3;
		nu_H4 = ro * ro_H4 * U_M_H4 * sigma_H4;

		k->par[0].M_u = (n_p_LISM_ / Kn_) * (nu_H1 * (u_H1 - u) + nu_H2 * (u_H2 - u) //
			+ nu_H3 * (u_H3 - u) + nu_H4 * (u_H4 - u));
		k->par[0].M_v = (n_p_LISM_ / Kn_) * (nu_H1 * (v_H1 - v) + nu_H2 * (v_H2 - v) //
			+ nu_H3 * (v_H3 - v) + nu_H4 * (v_H4 - v));

		
		k->par[0].M_T = (n_p_LISM_ / Kn_) * (nu_H1 * ((kv(u_H1) + kv(v_H1) - kv(u) - kv(v)) / 2.0 + //
			(U_H1 / U_M_H1) * (2.0 * p_H1 / ro_H1 - p / ro)) + //
			nu_H2 * ((kv(u_H2) + kv(v_H2) - kv(u) - kv(v)) / 2.0 + //
				(U_H2 / U_M_H2) * (2.0 * p_H2 / ro_H2 - p / ro)) + //
			nu_H3 * ((kv(u_H3) + kv(v_H3) - kv(u) - kv(v)) / 2.0 + //
				(U_H3 / U_M_H3) * (2.0 * p_H3 / ro_H3 - p / ro)) + //
			nu_H4 * ((kv(u_H4) + kv(v_H4) - kv(u) - kv(v)) / 2.0 + //
				(U_H4 / U_M_H4) * (2.0 * p_H4 / ro_H4 - p / ro)));

		
		k->par[0].k_u = k->par[0].I_u / k->par[0].M_u;
		if (fabs(k->par[0].k_u) > 10.0)
		{
			k->par[0].k_u = 1.0;
		}
		k->par[0].k_v = k->par[0].I_v / k->par[0].M_v;
		if (fabs(k->par[0].k_v) > 10.0)
		{
			k->par[0].k_v = 1.0;
		}
		k->par[0].k_T = k->par[0].I_T / k->par[0].M_T;
		if (fabs(k->par[0].k_T) > 10.0)
		{
			k->par[0].k_T = 1.0;
		}

		
	}
}

void Setka::culc_K_Istok(void)
{

	for (auto& k : this->All_Cells)
	{
		double U_M_H[4];
		double U_H[4];
		double sigma_H[4];
		double nu_H[4];
		double q2_1, q2_2, q3;
		double u, v, ro, p, Q;

		double u_H[4];
		double v_H[4];
		double ro_H[4];
		double p_H[4];

		for (int i = 0; i < 4; i++)
		{
			p_H[i] = 0.5 * k->par[0].H_T[i] * k->par[0].H_n[i];
			if (k->par[0].H_n[i] <= 0.0)
			{
				k->par[0].H_n[i] = 0.0000001;
				p_H[i] = 0.0;
			}
		}

		u = k->par[0].u;
		v = k->par[0].v;
		ro = k->par[0].ro;
		p = k->par[0].p;

		if (ro <= 0.0)
		{
			ro = 0.0000001;
			p = 0.0;
		}


		for (int i = 0; i < 4; i++)
		{
			U_M_H[i] = sqrt(kv(u - k->par[0].H_u[i]) + kv(v - k->par[0].H_v[i]) + (64.0 / (9.0 * pi_)) //
				* (p / ro + 2.0 * p_H[i] / k->par[0].H_n[i]));

			U_H[i] = sqrt(kv(u - k->par[0].H_u[i]) + kv(v - k->par[0].H_v[i]) + (4.0 / (pi_)) //
				* (p / ro + 2.0 * p_H[i] / k->par[0].H_n[i]));

			sigma_H[i] = kv(1.0 - a_2 * log(U_M_H[i]));

			nu_H[i] = ro * k->par[0].H_n[i] * U_M_H[i] * sigma_H[i];
		}

		k->par[0].M_u = 0.0;
		k->par[0].M_v = 0.0;
		k->par[0].M_T = 0.0;

		for (int i = 0; i < 4; i++)
		{
			k->par[0].M_u = k->par[0].M_u + (n_p_LISM_ / Kn_) * nu_H[i] * (k->par[0].H_u[i] - u);
			k->par[0].M_v = k->par[0].M_v + (n_p_LISM_ / Kn_) * nu_H[i] * (k->par[0].H_v[i] - v);
			k->par[0].M_T = k->par[0].M_T + (n_p_LISM_ / Kn_) * nu_H[i] * ((kv(k->par[0].H_u[i]) + kv(k->par[0].H_v[i]) - kv(u) - kv(v)) / 2.0 + //
				(U_H[i] / U_M_H[i]) * (2.0 * p_H[i] / k->par[0].H_n[i] - p / ro));
		}


		if (fabs(k->par[0].M_u) > 0.000001)
		{
			k->par[0].k_u = (k->par[0].I_u) / k->par[0].M_u;
			if (k->par[0].k_u > 30.0 || k->par[0].k_u < 0.01)
			{
				k->par[0].k_u = 1.0;
			}
		}
		else
		{
			k->par[0].k_u = 1.0;
		}

		if (fabs(k->par[0].M_v) > 0.000001)
		{
			k->par[0].k_v = (k->par[0].I_v) / k->par[0].M_v;
			if (k->par[0].k_v > 30.0 || k->par[0].k_v < 0.01)
			{
				k->par[0].k_v = 1.0;
			}
		}
		else
		{
			k->par[0].k_v = 1.0;
		}

		if (fabs(k->par[0].M_T) > 0.000001)
		{
			k->par[0].k_T = (k->par[0].I_T) / k->par[0].M_T;
			if (k->par[0].k_T > 30.0 || k->par[0].k_T < 0.01)
			{
				k->par[0].k_T = 1.0;
			}
		}
		else
		{
			k->par[0].k_T = 1.0;
		}

	}
}

void Setka::MK_start_new(void)
{

	mutex mut_1;
	double Y = fabs(Velosity_inf);
	int st = 1;

	ofstream fout;
	fout.open(parameter_7);

	ofstream fout2;
	fout2.open("stat_moments.txt");



	//omp_set_num_threads(12);

#pragma omp parallel for schedule(dynamic) // num_threads(12)
		for (int index = 0; index < 411; index++)
		{
			double A1 = 1.0 + (1.0 + 1.0 / (2.0 * kv(Y))) * erf(Y) + exp(-kv(Y)) / (Y * sqrtpi_);
			double A2 = 0.0; // Нужно задать после нахождения угла theta
			double mu1, mu2, mu3, mu4;
			double Vx, Vy, Vz;
			mu1 = ((this->sqv_1) / this->sum_s) * (1.0 * this->AllNumber / this->Number1);
			mu2 = ((this->sqv_2) / this->sum_s) * (1.0 * this->AllNumber / this->Number2);
			mu3 = ((this->sqv_3) / this->sum_s) * (1.0 * this->AllNumber / this->Number3);
			mu4 = ((this->sqv_4) / this->sum_s) * (1.0 * this->AllNumber / this->Number4);
			if (index == 0)
			{
				mut_1.lock();
				cout << "Setka.cpp    " << "Mk_start_new   " << "Wright weyght   " << mu1 << " " << mu2 << " " << mu3 << " " << mu4 << endl;
				mut_1.unlock();
			}

			// Печать статистики в файл
			//if(func_stat)
			//{
			//	if (omp_get_thread_num() == 0)
			//	{
			//		this->mut_stat.lock();
			//		for (int i = this->number_stat_2; i < this->number_stat; i++)
			//		{
			//			fout << V_r_stat[i] << " " << V_t_stat[i] << " " << V_p_stat[i] << " " << mu_stat[i] << " " << phi_stat[i] << " " << //
			//				num_stat[i] << endl;
			//		}

			//		this->number_stat = 0;
			//		this->number_stat_2 = this->number_stat;
			//		this->mut_stat.unlock();
			//	}
			//}

			Sensor* sens1 = Sensors[2 * index];
			Sensor* sens2 = Sensors[2 * index + 1];
			MKmethod MK = MKmethod();

			int s1 = sens2->a1_;
			int s2 = sens2->a2_;
			int s3 = sens2->a3_;

			mut_1.lock();
			cout << "Setka.cpp    " << "Start Mk_start_new   " << st << " potok  is  411;  index = " << index << "  nomer = " << omp_get_thread_num() << "   k1 = " << this->number_stat << endl;
			st++;
			mut_1.unlock();

			double x, y, z;
			double ksi1, ksi2, ksi3, ksi4, ksi5, ksi6, ksi7;
			double phi, Vphi, Vr;
			double a, b, c;
			double r;

			for (int ii = 0; ii < Number1 / 411; ii++)  //  Вылет с полусферы   Number1 / 135
			{
				if (true)
				{
					vector <double> mu(I_ + 1);
					vector <double> Wt(I_ + 1);
					vector <double> Wp(I_ + 1);
					vector <double> Wr(I_ + 1);
					vector <double> X(I_ + 1);
					double phi, sin_;
					//bool bb = MK.Init_Parametrs(sens1, mu, Wt, Wp, Wr, X);
					//cout << "Setka.cpp    " << "Mk_start_new   " << "A " << endl;
					bool bb = MK.Init_Parametrs(sens1, mu, Wt, Wp, Wr, X);
					//bb = true;
					//mu[I_] = 1.0;
					//cout << "Setka.cpp    " << "Mk_start_new   " << "B " << endl;

					for (int i = 0; i <= I_; i++) // I_  от 0 было
					{
						sin_ = sqrt(1.0 - kv(X[i]));
						x = (Rmax_ - 3.0 / RR_) * X[i];
						phi = 2.0 * pi_ * sens1->MakeRandom();
						y = (Rmax_ - 3.0 / RR_) * sin_ * cos(phi);
						z = (Rmax_ - 3.0 / RR_) * sin_ * sin(phi);
						Cell* Point = Belong_point(1, x, sqrt(kv(z) + kv(y)));
						//dekard_skorost2(Rmax_, acos(X[i]), phi, Wr[i], Wp[i], Wt[i], Vy, Vz, Vx);
						dekard_skorost(y, z, x, Wr[i], Wp[i], Wt[i], Vy, Vz, Vx);

						//Fly_exchenge_Imit_Korol(MK, sens2, x, y, z, Vx, Vy, Vz, Point, mu1 * mu[i], 3, false, mu1 * mu[i], i);

						if (i != I_ || bb == true)
						{
							double time_do_peregel, peregel;
							int ii_z, ii_alp;
							if (Vx * x + Vy * y + Vz * z < 0.0)
							{
								time_do_peregel = (-Vx * x - Vy * y - Vz * z) / kvv(Vx, Vy, Vz);
								peregel = sqrt(kvv(x + Vx * time_do_peregel, y + Vy * time_do_peregel, z + Vz * time_do_peregel));
								ii_z = this->geo_zones(peregel);                     // Номер зоны перегелия атома
								ii_alp = alpha_zones(x + Vx * time_do_peregel, sqrt(kvv(y + Vy * time_do_peregel, z + Vz * time_do_peregel, 0.0)));
							}
							else
							{
								ii_z = I_;
								ii_alp = alpha_zones(x, sqrt(kvv(y, z, 0.0)));
							}


							/*bool** BZ = new bool* [I_];
							for (size_t i2 = 0; i2 < I_; i2++)
							{
								BZ[i2] = new bool[J_];
							}

							for (size_t i2 = 0; i2 < I_; i2++)
							{
								for (size_t j2 = 0; j2 < J_; j2++)
								{
									BZ[i2][j2] = false;
								}
							}*/


							// Расчёт весов
							//cout << "A" << endl;
							//Fly_exchenge_Imit_Korol(MK, sens2, BZ, x, y, z, Vx, Vy, Vz, Point, mu1 * mu[i], 3, false, mu1 * mu[i], ii_z, ii_alp, true);
							//cout << "B" << endl;

							/*for (int i = 0; i < I_; ++i) {
								delete[] BZ[i];
							}
							delete[] BZ;*/


							// Эта функция, если веса посчитаны
							Fly_exchenge_Imit_Korol_auto_weight(MK, s1, s2, s3, x, y, z, Vx, Vy, Vz, Point, mu1 * mu[i],
								3, false, mu1 * mu[i], ii_z, ii_alp, true);


							// Работающая с заданными аналитически весами
							//Fly_exchenge_Imit_Korol(MK, s1, s2, s3, x, y, z, Vx, Vy, Vz, Point, mu1 * mu[i], 3, false, mu1 * mu[i], ii_z, ii_alp, true);


							//Fly_exchenge_Imit_Korol_PUI(MK, sens2, x, y, z, Vx, Vy, Vz, Point, mu1 * mu[i], 3, false, mu1, ii_z, ii_alp, true);  // mu1 * mu[i]


							// Эта для иммитационного метода
							//Fly_exchenge_Imit_Korol_2(MK, sens2, x, y, z, Vx, Vy, Vz,//
							//	Point, mu1 * mu[i], -log(1.0 - sens1->MakeRandom()), 0.0, 3, mu1 * mu[i]);

							//Fly_exchenge_Imit(MK, sens2, x, y, z, Vx, Vy, Vz, Point, mu1 * mu[i], -log(1.0 - sens1->MakeRandom()), 0.0, 3, mu1, i, ii); // i


							

						}
					}
				}

			}
			for (int ii = 0; ii < Number2 / 411; ii++)  //  Вылет сверху 
			{
				//cout << "Setka.cpp    " << "Mk_start_new   " << "A" << endl;
				ksi1 = sens1->MakeRandom();
				ksi2 = sens1->MakeRandom();
				ksi3 = sens1->MakeRandom();
				ksi4 = sens1->MakeRandom();
				ksi5 = sens1->MakeRandom();

				double ll = (Left_ + 2.0 / RR_);
				double rr = -0.5 / RR_;
				x = ll + ksi1 * (rr - ll);
				phi = ksi2 * 2.0 * pi_;
				Vphi = cos(2.0 * pi_ * ksi3) * sqrt(-log(1.0 - ksi4));
				Vx = Velosity_inf + sin(2.0 * pi_ * ksi3) * sqrt(-log(1.0 - ksi4));
				Vr = -sqrt(-log(ksi5));
				y = (R5_ - 2.0 / RR_) * cos(phi);
				z = (R5_ - 2.0 / RR_) * sin(phi);

				Cell* Point = Belong_point(2, x, sqrt(kv(z) + kv(y)));  // Находит ячейку, которой принадлежит точка
				//Fly_exchenge_Imit_Korol(sens2, x, y, z, Vx, cos(phi) * Vr - sin(phi) * Vphi, sin(phi) * Vr + cos(phi) * Vphi,//
				//	Point, mu2, 3, false, mu2, I_ - 1);

				//Fly_exchenge_Imit_Korol(MK, sens2, x, y, z, Vx, cos(phi) * Vr - sin(phi) * Vphi, sin(phi) * Vr + cos(phi) * Vphi,//
				//	Point, mu2, 3, false, mu2, I_ - 1, J_ - 1, false);

				Fly_exchenge_Imit_Korol_2(MK, sens2, x, y, z, Vx, cos(phi) * Vr - sin(phi) * Vphi, sin(phi) * Vr + cos(phi) * Vphi,//
					Point, mu2, -log(1.0 - sens1->MakeRandom()), 0.0, 3, mu2);

				//Fly_exchenge_Imit(MK, sens2, x, y, z, Vx, cos(phi) * Vr - sin(phi) * Vphi, sin(phi) * Vr + cos(phi) * Vphi,//
				//	Point, mu2, -log(1.0 - sens1->MakeRandom()), 0.0, 3, mu2, I_ - 1, ii);
			}
			for (int ii = 0; ii < Number3 / 411; ii++)  //  Вылет сзади
			{
				//cout << "Setka.cpp    " << "Mk_start_new   " << "B" << endl;
				Velosity_initial2(sens1, a, b, c);
				ksi1 = sens1->MakeRandom();
				ksi2 = sens1->MakeRandom();

				r = sqrt(ksi1 * (R5_ - 2.0 / RR_) * (R5_ - 2.0 / RR_));
				phi = ksi2 * 2.0 * pi_;
				y = r * cos(phi);
				z = r * sin(phi);

				Cell* Point = Belong_point(3, Left_ + 2.0 / RR_, sqrt(kv(z) + kv(y)));  // Находит ячейку, которой принадлежит точка
				//cout << "Setka.cpp    " << "Mk_start_new   " << "Start" << endl;
				//Fly_exchenge_Imit_Korol(sens2, Left_ + 2.0, y, z, a, b, c, Point, mu3, 3, false, mu3, I_ - 1);

				//Fly_exchenge_Imit_Korol_PUI(MK, sens2, Left_ + 2.0 / RR_, y, z, a, b, c, Point, mu3, 3, false, mu3, I_ - 1, J_ - 1, false);

				//Fly_exchenge_Imit(MK, sens2, Left_ + 2.0/ RR_, y, z, a, b, c, Point, mu3, -log(1.0 - sens1->MakeRandom()), 0.0, //
				//	3, mu3, I_ - 1, ii);

				Fly_exchenge_Imit_Korol_2(MK, sens2, Left_ + 2.0 / RR_, y, z, a, b, c,//
					Point, mu3, -log(1.0 - sens1->MakeRandom()), 0.0, 3, mu3);
			}
			for (int ii = 0; ii < Number4 / 411; ii++)  //Number4 / 411
			{
				double a, b, c;
				Velosity_initial(sens1, a, b, c);
				ksi1 = sens1->MakeRandom();
				ksi2 = sens1->MakeRandom();

				r = sqrt(ksi1 * (kv(R5_ - 2.0 / RR_) - kv(Rmax_ + 1.0 / RR_)) + kv(Rmax_ + 1.0 / RR_));  // 0.1
				phi = ksi2 * 2.0 * pi_;
				y = r * cos(phi);
				z = r * sin(phi);

				Cell* Point = Belong_point(4, -0.5 / RR_, sqrt(kv(z) + kv(y)));

				//Fly_exchenge_Imit_Korol(MK, sens2, -0.1/RR_, y, z, a, b, c, Point, mu4, 3, false, mu4, I_ - 1, J_ - 1, true);

				Fly_exchenge_Imit_Korol_2(MK, sens2, -0.5 / RR_, y, z, a, b, c, Point, mu4, -log(1.0 - sens1->MakeRandom()), //
					0.0, 3, mu4);
				//Fly_exchenge_Imit(MK, sens2, -0.01 / RR_, y, z, a, b, c, Point, mu4, -log(1.0 - sens1->MakeRandom()), 0.0, //
				//	3, mu4, I_ - 1, ii);
			}

			/*mut_1.lock();
			cout << "Setka.cpp    " << "End Mk_start_new   " << st << " potok  is  411;  index = " << index << "  nomer = " << omp_get_thread_num() << "   k1 = " << this->number_stat << endl;
			mut_1.unlock();*/
		}

	
	// Делаем веса монотонными по радиусу
	for (size_t k = 0; k < 4; k++)
	{
		for (size_t i = 1; i < I_; i++)
		{
			for (size_t j = 0; j < J_; j++)
			{
				if (Mu_statistic[k][i][j] < Mu_statistic[k][i - 1][j])
				{
					Mu_statistic[k][i][j] = Mu_statistic[k][i - 1][j];
				}
			}
		}
	}

	for (size_t k = 0; k < 4; k++)
	{
		for (size_t i = 0; i < I_; i++)
		{
			for (size_t j = 0; j < J_; j++)
			{
				std::cout << std::fixed;
				std::cout << std::setprecision(13);
				fout << k << " " << i << " " << j << " " << Mu_statistic[k][i][j] * kv(Rmax_/Ri[i]) * J_/ this->AllNumber << endl;
			}
		}
	}

	fout.close();

	//for (auto& k : this->All_Cells)

#pragma omp parallel for
		for(int ikld = 0; ikld < this->All_Cells.size(); ikld++)
	{
		/*for (int il = 0; il < 7; il++)
		{
			if (k->par[0].num_atoms > 0)
			{
				k->par[0].w_m[il] /= k->par[0].num_atoms;
			}
		}*/

		auto k = this->All_Cells[ikld];

		double no = (1.0 * AllNumber * k->Get_Volume_rotate(360.0));

		k->par[0].I_u = sum_s * k->par[0].I_u / no;
		k->par[0].I_v = sum_s * k->par[0].I_v / no;
		k->par[0].I_T = sum_s * k->par[0].I_T / no;

		k->par[0].II_u = sum_s * k->par[0].II_u / no;
		k->par[0].II_v = sum_s * k->par[0].II_v / no;
		k->par[0].II_T = sum_s * k->par[0].II_T / no;



		if (k->par[0].F_n > 0)
		{
			/*k->par[0].I_u = k->par[0].I_u / k->par[0].F_n;
			k->par[0].I_v = k->par[0].I_v / k->par[0].F_n;
			k->par[0].I_T = k->par[0].I_T / k->par[0].F_n;*/

			k->par[0].F_u = k->par[0].F_u / k->par[0].F_n;
			k->par[0].F_v = k->par[0].F_v / k->par[0].F_n;
			k->par[0].F_T = (2.0 / 3.0) * (k->par[0].F_T / k->par[0].F_n - kvv(k->par[0].F_u, k->par[0].F_v, 0.0));

			for (int ij = 0; ij < n_S; ij++)
			{
				double  w2 = (ij + 1) * (max_S / n_S);
				double  w1 = (ij) * (max_S / n_S);
				k->S_p[ij] = sum_s * k->S_p[ij] / (no * 4.0 * pi_ * (1.0 / 3.0) * (pow(w2, 3) - pow(w1, 3)));
				k->S_m[ij] = sum_s * k->S_m[ij] / (no);
			}
		}
		else
		{
			/*k->par[0].I_u = 0.0;
			k->par[0].I_v = 0.0;
			k->par[0].I_T = 0.0;*/

			k->par[0].F_n = 0.0;
			k->par[0].F_u = 0.0;
			k->par[0].F_v = 0.0;
			k->par[0].F_T = 0.0;

			for (int ij = 0; ij < n_S; ij++)
			{
				k->S_m[ij] = 0.0;
				k->S_p[ij] = 0.0;
			}

		}

		// Теперь посчиатем S- для игоря ///////////////////////////////////////////////
		if (false)//(k->type == C_1 || k->type == C_2 || k->type == C_3)
		{
			double S_min[n_S];
			for (int ij = 0; ij < n_S; ij++)
			{
				S_min[ij] = 0.0;
			}

			double dthe = pi_ / 40.0;

			for (int ij = 0; ij < n_S; ij++)
			{
				double w = ((ij + 1.0) * max_S / n_S + ij * max_S / n_S) / 2.0;
				double w2 = ((ij + 1.0) * max_S / n_S);
				double w1 = ((ij) * max_S / n_S);
				for (int ik = 0; ik < n_S; ik++)
				{
					double Vh = ((ik + 1.0) * max_S / n_S + ik * max_S / n_S) / 2.0;
					double ff = k->S_m[ik];
					if (ff > 0.0)
					{
						for (double the = 0.0; the <= pi_; the = the + dthe)
						{
							double d = sqrt(kv(Vh) + kv(w) - 2.0 * w * Vh * cos(the));
							if (d > 0.000000001)
							{
								S_min[ij] = S_min[ij] + ff * d * sigma(d) * sin(the) * dthe * 2.0 * pi_; //  / Kn_
							}
						}
					}

				}
				S_min[ij] = S_min[ij] / (4.0 * pi_);
			}


			for (int ij = 0; ij < n_S; ij++)
			{
				k->S_m[ij] = S_min[ij];
			}
		}

		////////////////////////////////////////////////////////////////////////////////

		if (k->par[0].F_n > 0)
		{
			for (int i = 0; i < pop_atom; i++)
			{
				if (k->par[0].H_n[i] > 0)
				{
					k->par[0].H_u[i] = k->par[0].H_u[i] / k->par[0].H_n[i];
					k->par[0].H_v[i] = k->par[0].H_v[i] / k->par[0].H_n[i];
					//k->par[0].H_T[i] = (2.0 / 3.0) * (k->par[0].H_T[i] / k->par[0].H_n[i] - kvv(k->par[0].F_u, k->par[0].F_v, 0.0));
					k->par[0].H_T[i] = (2.0 / 3.0) * (k->par[0].H_T[i] / k->par[0].H_n[i] - kvv(k->par[0].H_u[i], k->par[0].H_v[i], 0.0));
					k->par[0].H_uu[i] = k->par[0].H_uu[i] / k->par[0].H_n[i] - k->par[0].H_u[i] * k->par[0].H_u[i];
					k->par[0].H_vv[i] = k->par[0].H_vv[i] / k->par[0].H_n[i] - k->par[0].H_v[i] * k->par[0].H_v[i];
					k->par[0].H_uv[i] = k->par[0].H_uv[i] / k->par[0].H_n[i] - k->par[0].H_u[i] * k->par[0].H_v[i];
					k->par[0].H_uuu[i] = 2.0 * k->par[0].H_v[i] * k->par[0].H_v[i] * k->par[0].H_v[i] + //
						3.0 * k->par[0].H_v[i] * k->par[0].H_vv[i] - k->par[0].H_uuu[i] / k->par[0].H_n[i];
				}
				else
				{
					k->par[0].H_u[i] = 0.0;
					k->par[0].H_v[i] = 0.0;
					k->par[0].H_T[i] = 0.0;
					k->par[0].H_uu[i] = 0.0;
					k->par[0].H_vv[i] = 0.0;
					k->par[0].H_uv[i] = 0.0;
					k->par[0].H_uuu[i] = 0.0;
				}
			}
		}

		if (k->df_s4_bool == true)
		{
			k->df_s4->normir(sum_s / no);
		}

		if (k->df_s3_bool == true)
		{
			k->df_s3->normir(sum_s / no);
		}

		if (k->df_s2_bool == true)
		{
			k->df_s2->normir(sum_s / no);
		}

		if (k->df_s1_bool == true)
		{
			k->df_s1->normir(sum_s / no);
		}

		k->par[0].F_n = sum_s * k->par[0].F_n / no;
		for (int i = 0; i < pop_atom; i++)
		{
			k->par[0].H_n[i] = sum_s * k->par[0].H_n[i] / no;
		}


		k->par[0].I_u = (n_p_LISM_) * (k->par[0].I_u);
		k->par[0].I_v = (n_p_LISM_) * (k->par[0].I_v);
		k->par[0].I_T = (n_p_LISM_) * (k->par[0].I_T);
		k->par[0].II_u = (n_p_LISM_) * (k->par[0].II_u);
		k->par[0].II_v = (n_p_LISM_) * (k->par[0].II_v);
		k->par[0].II_T = (n_p_LISM_) * (k->par[0].II_T);

		// Считаем мультифдюидные источники

		double U_M_H[pop_atom];
		double U_H[pop_atom];
		double sigma_H[pop_atom];
		double nu_H[pop_atom];
		double q2_1, q2_2, q3;
		double u, v, ro, p, Q;

		double u_H[pop_atom];
		double v_H[pop_atom];
		double ro_H[pop_atom];
		double p_H[pop_atom];

		for (int i = 0; i < pop_atom; i++)
		{
			p_H[i] = 0.5 * k->par[0].H_T[i] * k->par[0].H_n[i];
			if (k->par[0].H_n[i] <= 0.0)
			{
				k->par[0].H_n[i] = 0.0000001;
				p_H[i] = 0.0;
			}
		}

		if (k->pui_ == true)
		{
			k->par[0].ro = k->par[0].ro + k->par[0].npui;
			double dd = k->par[0].p;
			k->par[0].p = k->par[0].pp;
			k->par[0].pp = dd;
		}

		u = k->par[0].u;
		v = k->par[0].v;
		ro = k->par[0].ro;
		p = k->par[0].p;

		if (ro <= 0.0)
		{
			ro = 0.0000001;
			p = 0.0;
		}


		for (int i = 0; i < pop_atom; i++)
		{
			U_M_H[i] = sqrt(kv(u - k->par[0].H_u[i]) + kv(v - k->par[0].H_v[i]) + (64.0 / (9.0 * pi_)) //
				* (p / ro + 2.0 * p_H[i] / k->par[0].H_n[i]));

			U_H[i] = sqrt(kv(u - k->par[0].H_u[i]) + kv(v - k->par[0].H_v[i]) + (4.0 / (pi_)) //
				* (p / ro + 2.0 * p_H[i] / k->par[0].H_n[i]));

			sigma_H[i] = kv(1.0 - a_2 * log(U_M_H[i]));

			nu_H[i] = ro * k->par[0].H_n[i] * U_M_H[i] * sigma_H[i];
		}
		
		k->par[0].M_u = 0.0;
		k->par[0].M_v = 0.0;
		k->par[0].M_T = 0.0;

		for (int i = 0; i < pop_atom; i++)
		{
			k->par[0].M_u = k->par[0].M_u + (n_p_LISM_ / Kn_) * nu_H[i] * (k->par[0].H_u[i] - u);
			k->par[0].M_v = k->par[0].M_v + (n_p_LISM_ / Kn_) * nu_H[i] * (k->par[0].H_v[i] - v);
			k->par[0].M_T = k->par[0].M_T + (n_p_LISM_ / Kn_) * nu_H[i] * ((kv(k->par[0].H_u[i]) + kv(k->par[0].H_v[i]) - kv(u) - kv(v)) / 2.0 + //
				(U_H[i] / U_M_H[i]) * (2.0 * p_H[i] / k->par[0].H_n[i] - p / ro));
		}


		if (fabs(k->par[0].M_u) > 0.000001)
		{
			k->par[0].k_u = (k->par[0].I_u) / k->par[0].M_u;
			if (k->par[0].k_u > 10.0 || k->par[0].k_u < 0.05)
			{
				k->par[0].k_u = 1.0;
			}
		}
		else
		{
			k->par[0].k_u = 1.0;
		}

		if (fabs(k->par[0].M_v) > 0.000001)
		{
			k->par[0].k_v = (k->par[0].I_v) / k->par[0].M_v;
			if (k->par[0].k_v > 10.0 || k->par[0].k_v < 0.05)
			{
				k->par[0].k_v = 1.0;
			}
		}
		else
		{
			k->par[0].k_v = 1.0;
		}

		if (fabs(k->par[0].M_T) > 0.000001)
		{
			k->par[0].k_T = (k->par[0].I_T) / k->par[0].M_T;
			if (k->par[0].k_T > 30.0 || k->par[0].k_T < 0.05)
			{
				k->par[0].k_T = 1.0;
			}
		}
		else
		{
			k->par[0].k_T = 1.0;
		}

		double xx, yy;
		k->Get_Center(xx, yy);
		if (false)//(xx > 2.6 && k->axis_ == true)
		{
			k->par[0].F_n = 1.0;
			k->par[0].F_u = Velosity_inf;
			k->par[0].F_v = 0.0;
			k->par[0].F_T = 1.0;
			k->par[0].H_n[0] = 0.0;
			k->par[0].H_u[0] = 0.0;
			k->par[0].H_v[0] = 0.0;
			k->par[0].H_T[0] = 0.0;
			k->par[0].H_n[1] = 0.0;
			k->par[0].H_u[1] = 0.0;
			k->par[0].H_v[1] = 0.0;
			k->par[0].H_T[1] = 0.0;
			k->par[0].H_n[2] = 0.0;
			k->par[0].H_u[2] = 0.0;
			k->par[0].H_v[2] = 0.0;
			k->par[0].H_T[2] = 0.0;
			k->par[0].H_n[3] = 1.0;
			k->par[0].H_u[3] = Velosity_inf;
			k->par[0].H_v[3] = 0.0;
			k->par[0].H_T[3] = 1.0;

			k->par[0].I_u = 0.0;
			k->par[0].I_v = 0.0;
			k->par[0].I_T = 0.0;

			k->par[0].k_u = 0.0;
			k->par[0].k_v = 0.0;
			k->par[0].k_T = 0.0;
		}

	}

	for (auto& k : this->Cell_other)
	{
		k->par[0].F_n = 1.0;
		k->par[0].F_u = Velosity_inf;
		k->par[0].F_v = 0.0;
		k->par[0].F_T = 1.0;
		k->par[0].H_n[0] = 0.0;
		k->par[0].H_u[0] = 0.0;
		k->par[0].H_v[0] = 0.0;
		k->par[0].H_T[0] = 0.0;
		k->par[0].H_n[1] = 0.0;
		k->par[0].H_u[1] = 0.0;
		k->par[0].H_v[1] = 0.0;
		k->par[0].H_T[1] = 0.0;
		k->par[0].H_n[2] = 0.0;
		k->par[0].H_u[2] = 0.0;
		k->par[0].H_v[2] = 0.0;
		k->par[0].H_T[2] = 0.0;
		k->par[0].H_n[3] = 1.0;
		k->par[0].H_u[3] = Velosity_inf;
		k->par[0].H_v[3] = 0.0;
		k->par[0].H_T[3] = 1.0;

		k->par[0].I_u = 0.0;
		k->par[0].I_v = 0.0;
		k->par[0].I_T = 0.0;

		k->par[0].k_u = 0.0;
		k->par[0].k_v = 0.0;
		k->par[0].k_T = 0.0;
	}


	if (func_stat)
	{
		double no = (1.0 * AllNumber);
		fout2 << sum_s / no << endl;

		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < Al_stat; j++)
			{
				fout2 << this->mu_mom[i][j] << " " << this->Vx_mom[i][j] << " " << this->Vy_mom[i][j] //
					<< " " << this->Vxx_mom[i][j] << " " << this->Vyy_mom[i][j] << " " << this->Vxy_mom[i][j] //
					<< " " << this->Vxxx_mom[i][j] << " " << this->T_mom[i][j] << endl;
			}
		}
	}


	for (auto& kk : this->All_Cells)
	{
		double no = (1.0 * AllNumber * kk->Get_Volume_rotate(360.0));

		for (int k = 0; k < 4; k++)
		{
			for (int j = 0; j < pogl_rad_; j++)
			{
				kk->pogloshenie[k][j] = kk->pogloshenie[k][j] * this->sum_s / no;
			}
		}
	}
}

void Setka::MPI_MK_start(int argc, char** argv)
{
#if USEMPI 
	cout << "START  MPI_MK_start" << endl;
	mutex mut_1;
	double Y = fabs(Velosity_inf);
	int st = 1;
	int rank = 0, size = 0;
	int ierr;

	MPI_Init(&argc, &argv);


	MPI_Comm_size(MPI_COMM_WORLD, &size);               // Получить общее число процессов - компов
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);               // Получить номер текущего процесса - компа

	ofstream fout2;
	if (rank == 0)
	{
		fout2.open("stat_moments.txt");
	}


	cout << "START  MPI_MK_start  from  " << rank << endl;

	// omp_set_num_threads(24);

	/*if (omp_get_num_procs() < 28)
	{
		cout << "omp_get_num_procs  =  " << omp_get_num_procs() << "   KILL" << endl;
		exit(-12);
	}*/

#pragma omp parallel for schedule(dynamic) // num_threads(12)
	for (int index = 0; index < 240; index++)
	{
		double A1 = 1.0 + (1.0 + 1.0 / (2.0 * kv(Y))) * erf(Y) + exp(-kv(Y)) / (Y * sqrtpi_);
		double A2 = 0.0; // Нужно задать после нахождения угла theta
		double mu1, mu2, mu3, mu4;
		double Vx, Vy, Vz;
		mu1 = ((this->sqv_1) / this->sum_s) * (1.0 * this->AllNumber / this->Number1);
		mu2 = ((this->sqv_2) / this->sum_s) * (1.0 * this->AllNumber / this->Number2);
		mu3 = ((this->sqv_3) / this->sum_s) * (1.0 * this->AllNumber / this->Number3);
		mu4 = ((this->sqv_4) / this->sum_s) * (1.0 * this->AllNumber / this->Number4);
		int smax = omp_get_num_procs();
		mut_1.lock();
		Sensor* sens1 = Sensors[(28 * rank + omp_get_thread_num()) % 1021];
		Sensor* sens2 = Sensors[(50 * omp_get_thread_num() + rank) % 1021];
		//Sensor* sens1 = Sensors[(rank) % 1021];
		//Sensor* sens2 = Sensors[(rank + 499) % 1021];
		MKmethod MK = MKmethod();
		int s1 = sens2->a1_;
		int s2 = sens2->a2_;
		int s3 = sens2->a3_;
		mut_1.unlock();

		mut_1.lock();
		cout << "Setka.cpp    " << "Mk_start_new   " << st << " potok  is  280;  index = " << index << "  nomer = " << omp_get_thread_num() << "  rank = " << rank << endl;
		st++;
		mut_1.unlock();

		double x, y, z;
		double ksi1, ksi2, ksi3, ksi4, ksi5, ksi6, ksi7;
		double phi, Vphi, Vr;
		double a, b, c;
		double r;

		for (int ii = 0; ii < Number1 / 240; ii++)  //  Вылет с полусферы   Number1 / 135
		{
			if (true)
			{
				vector <double> mu(I_ + 1);
				vector <double> Wt(I_ + 1);
				vector <double> Wp(I_ + 1);
				vector <double> Wr(I_ + 1);
				vector <double> X(I_ + 1);
				double phi, sin_;
				//bool bb = MK.Init_Parametrs(sens1, mu, Wt, Wp, Wr, X);
				//cout << "Setka.cpp    " << "Mk_start_new   " << "A " << endl;
				bool bb = MK.Init_Parametrs(sens1, mu, Wt, Wp, Wr, X);
				//bb = true;
				//mu[I_] = 1.0;
				//cout << "Setka.cpp    " << "Mk_start_new   " << "B " << endl;

				for (int i = 0; i <= I_; i++) // I_  от 0 было
				{
					sin_ = sqrt(1.0 - kv(X[i]));
					x = (Rmax_ - 3.0 / RR_) * X[i];
					phi = 2.0 * pi_ * sens1->MakeRandom();
					y = (Rmax_ - 3.0 / RR_) * sin_ * cos(phi);
					z = (Rmax_ - 3.0 / RR_) * sin_ * sin(phi);
					Cell* Point = Belong_point(1, x, sqrt(kv(z) + kv(y)));
					//dekard_skorost2(Rmax_, acos(X[i]), phi, Wr[i], Wp[i], Wt[i], Vy, Vz, Vx);
					dekard_skorost(y, z, x, Wr[i], Wp[i], Wt[i], Vy, Vz, Vx);

					//Fly_exchenge_Imit_Korol(MK, sens2, x, y, z, Vx, Vy, Vz, Point, mu1 * mu[i], 3, false, mu1 * mu[i], i);

					if (i != I_ || bb == true)
					{
						double time_do_peregel, peregel;
						int ii_z, ii_alp;
						if (Vx * x + Vy * y + Vz * z < 0.0)
						{
							time_do_peregel = (-Vx * x - Vy * y - Vz * z) / kvv(Vx, Vy, Vz);
							peregel = sqrt(kvv(x + Vx * time_do_peregel, y + Vy * time_do_peregel, z + Vz * time_do_peregel));
							ii_z = this->geo_zones(peregel);                     // Номер зоны перегелия атома
							ii_alp = alpha_zones(x + Vx * time_do_peregel, sqrt(kvv(y + Vy * time_do_peregel, z + Vz * time_do_peregel, 0.0)));
						}
						else
						{
							ii_z = I_;
							ii_alp = alpha_zones(x, sqrt(kvv(y, z, 0.0)));
						}

						/*bool** BZ = new bool* [I_];
						for (size_t i2 = 0; i2 < I_; i2++)
						{
							BZ[i2] = new bool[J_];
						}

						for (size_t i2 = 0; i2 < I_; i2++)
						{
							for (size_t j2 = 0; j2 < J_; j2++)
							{
								BZ[i2][j2] = false;
							}
						}*/

						//cout << "A" << endl;
						//Fly_exchenge_Imit_Korol(MK, sens2, BZ, x, y, z, Vx, Vy, Vz, Point, mu1 * mu[i], 3, false, mu1 * mu[i], ii_z, ii_alp, true);

						//Fly_exchenge_Imit_Korol_auto_weight(MK, s1, s2, s3, x, y, z, Vx, Vy, Vz, Point, mu1 * mu[i],
						//	3, false, mu1 * mu[i], ii_z, ii_alp, true);

						//Fly_exchenge_Imit_Korol(MK, s1, s2, s3, x, y, z, Vx, Vy, Vz, Point, mu1 * mu[i], 3, false, mu1 * mu[i], ii_z, ii_alp, true);


						//Fly_exchenge_Imit_Korol_PUI(MK, sens2, x, y, z, Vx, Vy, Vz, Point, mu1 * mu[i], 3, false, mu1, ii_z, ii_alp, true);  // mu1 * mu[i]

						Fly_exchenge_Imit_Korol_2(MK, sens2, x, y, z, Vx, Vy, Vz,//
							Point, mu1 * mu[i], -log(1.0 - sens1->MakeRandom()), 0.0, 3, mu1 * mu[i]);

						//Fly_exchenge_Imit(MK, sens2, x, y, z, Vx, Vy, Vz, Point, mu1 * mu[i], -log(1.0 - sens1->MakeRandom()), 0.0, 3, mu1, i, ii); // i

						/*for (int i = 0; i < I_; ++i) {
							delete[] BZ[i];
						}
						delete[] BZ;*/

					}
				}
			}

		}
		for (int ii = 0; ii < Number2 / 240; ii++)  //  Вылет сверху 
		{
			//cout << "Setka.cpp    " << "Mk_start_new   " << "A" << endl;
			ksi1 = sens1->MakeRandom();
			ksi2 = sens1->MakeRandom();
			ksi3 = sens1->MakeRandom();
			ksi4 = sens1->MakeRandom();
			ksi5 = sens1->MakeRandom();

			double ll = (Left_ + 2.0 / RR_);
			double rr = -0.5 / RR_;
			x = ll + ksi1 * (rr - ll);
			phi = ksi2 * 2.0 * pi_;
			Vphi = cos(2.0 * pi_ * ksi3) * sqrt(-log(1.0 - ksi4));
			Vx = Velosity_inf + sin(2.0 * pi_ * ksi3) * sqrt(-log(1.0 - ksi4));
			Vr = -sqrt(-log(ksi5));
			y = (R5_ - 2.0 / RR_) * cos(phi);
			z = (R5_ - 2.0 / RR_) * sin(phi);

			Cell* Point = Belong_point(2, x, sqrt(kv(z) + kv(y)));  // Находит ячейку, которой принадлежит точка
			//Fly_exchenge_Imit_Korol(sens2, x, y, z, Vx, cos(phi) * Vr - sin(phi) * Vphi, sin(phi) * Vr + cos(phi) * Vphi,//
			//	Point, mu2, 3, false, mu2, I_ - 1);

			//Fly_exchenge_Imit_Korol(MK, sens2, x, y, z, Vx, cos(phi) * Vr - sin(phi) * Vphi, sin(phi) * Vr + cos(phi) * Vphi,//
			//	Point, mu2, 3, false, mu2, I_ - 1, J_ - 1, false);

			Fly_exchenge_Imit_Korol_2(MK, sens2, x, y, z, Vx, cos(phi) * Vr - sin(phi) * Vphi, sin(phi) * Vr + cos(phi) * Vphi,//
				Point, mu2, -log(1.0 - sens1->MakeRandom()), 0.0, 3, mu2);

			//Fly_exchenge_Imit(MK, sens2, x, y, z, Vx, cos(phi) * Vr - sin(phi) * Vphi, sin(phi) * Vr + cos(phi) * Vphi,//
			//	Point, mu2, -log(1.0 - sens1->MakeRandom()), 0.0, 3, mu2, I_ - 1, ii);
		}
		for (int ii = 0; ii < Number3 / 240; ii++)  //  Вылет сзади
		{
			//cout << "Setka.cpp    " << "Mk_start_new   " << "B" << endl;
			Velosity_initial2(sens1, a, b, c);
			ksi1 = sens1->MakeRandom();
			ksi2 = sens1->MakeRandom();

			r = sqrt(ksi1 * (R5_ - 2.0 / RR_) * (R5_ - 2.0 / RR_));
			phi = ksi2 * 2.0 * pi_;
			y = r * cos(phi);
			z = r * sin(phi);

			Cell* Point = Belong_point(3, Left_ + 2.0 / RR_, sqrt(kv(z) + kv(y)));  // Находит ячейку, которой принадлежит точка
			//cout << "Setka.cpp    " << "Mk_start_new   " << "Start" << endl;
			//Fly_exchenge_Imit_Korol(sens2, Left_ + 2.0, y, z, a, b, c, Point, mu3, 3, false, mu3, I_ - 1);

			//Fly_exchenge_Imit_Korol_PUI(MK, sens2, Left_ + 2.0 / RR_, y, z, a, b, c, Point, mu3, 3, false, mu3, I_ - 1, J_ - 1, false);

			//Fly_exchenge_Imit(MK, sens2, Left_ + 2.0/ RR_, y, z, a, b, c, Point, mu3, -log(1.0 - sens1->MakeRandom()), 0.0, //
			//	3, mu3, I_ - 1, ii);

			//Fly_exchenge_Imit_Korol_2(MK, sens2, Left_ + 2.0 / RR_, y, z, a, b, c,//
			//	Point, mu3, -log(1.0 - sens1->MakeRandom()), 0.0, 3, mu3);                // Была эта, потом решили поменять
			 
			Fly_exchenge_Imit_Korol_auto_weight(MK, s1, s2, s3, Left_ + 2.0 / RR_, y, z, a, b, c, Point, mu3,
				3, false, mu3, 0, 0, true);
		}
		for (int ii = 0; ii < Number4 / 240; ii++)  //Number4 / 411
		{
			double a, b, c;
			Velosity_initial(sens1, a, b, c);
			ksi1 = sens1->MakeRandom();
			ksi2 = sens1->MakeRandom();

			r = sqrt(ksi1 * (kv(R5_ - 2.0 / RR_) - kv(Rmax_ - 0.1 / RR_)) + kv(Rmax_ - 0.1 / RR_));  // 0.1
			phi = ksi2 * 2.0 * pi_;
			y = r * cos(phi);
			z = r * sin(phi);

			Cell* Point = Belong_point(4, -0.5 / RR_, sqrt(kv(z) + kv(y)));

			if (Point == nullptr)
			{
				cout << "MIMO  odin  44444" << endl;
				continue;
			}

			//Fly_exchenge_Imit_Korol(MK, sens2, -0.1/RR_, y, z, a, b, c, Point, mu4, 3, false, mu4, I_ - 1, J_ - 1, true);

			Fly_exchenge_Imit_Korol_2(MK, sens2, -0.5 / RR_, y, z, a, b, c, Point, mu4, -log(1.0 - sens1->MakeRandom()), //
				0.0, 3, mu4);
			//Fly_exchenge_Imit(MK, sens2, -0.01 / RR_, y, z, a, b, c, Point, mu4, -log(1.0 - sens1->MakeRandom()), 0.0, //
			//	3, mu4, I_ - 1, ii);
		}
		
		/*mut_1.lock();
		sens2->a1_ = s1;
		sens2->a2_ = s2;
		sens2->a3_ = s3;
		mut_1.unlock();*/
}
	

		MPI_Barrier(MPI_COMM_WORLD);
		cout << "END calculate " << rank << endl;
		MPI_Barrier(MPI_COMM_WORLD);

		double* all_;

		double* all_mom;
		double* one_mom;

		double* one_double_all;
		double* ss_d;

		ss_d = (double*)malloc(sizeof(double));
		one_mom = (double*)malloc(54 * sizeof(double));

		if (rank == 0)
		{
			all_ = (double*)malloc(size * this->All_Cells.size() * sizeof(double));
			one_double_all = (double*)malloc(size * sizeof(double));

			all_mom = (double*)malloc(size * 54 * sizeof(double));
		}

		double* my_send;
		my_send = (double*)malloc(this->All_Cells.size() * sizeof(double));


		// Пересылка данных
		if (true)  
		{
			// Передача моментов на сфере
			for (int num = 0; num < 4; num++)  
			{
				for (int ikld = 0; ikld < 54; ikld++)
				{
					one_mom[ikld] = mu_mom[num][ikld];
				}
				MPI_Gather(one_mom, 54, MPI_DOUBLE, all_mom,
					54, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				if (rank == 0)
				{
					for (int ik = 0; ik < 54; ik++)
					{
						mu_mom[num][ik] = 0.0;
						for (int ikld = 0; ikld < size; ikld++)
						{
							mu_mom[num][ik] = mu_mom[num][ik] + all_mom[ik + ikld * 54];
						}
					}
				}
				///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				for (int ikld = 0; ikld < 54; ikld++)
				{
					one_mom[ikld] = Vx_mom[num][ikld];
				}
				MPI_Gather(one_mom, 54, MPI_DOUBLE, all_mom,
					54, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				if (rank == 0)
				{
					for (int ik = 0; ik < 54; ik++)
					{
						Vx_mom[num][ik] = 0.0;
						for (int ikld = 0; ikld < size; ikld++)
						{
							Vx_mom[num][ik] = Vx_mom[num][ik] + all_mom[ik + ikld * 54];
						}
					}
				}
				///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				for (int ikld = 0; ikld < 54; ikld++)
				{
					one_mom[ikld] = Vy_mom[num][ikld];
				}
				MPI_Gather(one_mom, 54, MPI_DOUBLE, all_mom,
					54, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				if (rank == 0)
				{
					for (int ik = 0; ik < 54; ik++)
					{
						Vy_mom[num][ik] = 0.0;
						for (int ikld = 0; ikld < size; ikld++)
						{
							Vy_mom[num][ik] = Vy_mom[num][ik] + all_mom[ik + ikld * 54];
						}
					}
				}
				///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				for (int ikld = 0; ikld < 54; ikld++)
				{
					one_mom[ikld] = Vxx_mom[num][ikld];
				}
				MPI_Gather(one_mom, 54, MPI_DOUBLE, all_mom,
					54, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				if (rank == 0)
				{
					for (int ik = 0; ik < 54; ik++)
					{
						Vxx_mom[num][ik] = 0.0;
						for (int ikld = 0; ikld < size; ikld++)
						{
							Vxx_mom[num][ik] = Vxx_mom[num][ik] + all_mom[ik + ikld * 54];
						}
					}
				}
				///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				for (int ikld = 0; ikld < 54; ikld++)
				{
					one_mom[ikld] = Vyy_mom[num][ikld];
				}
				MPI_Gather(one_mom, 54, MPI_DOUBLE, all_mom,
					54, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				if (rank == 0)
				{
					for (int ik = 0; ik < 54; ik++)
					{
						Vyy_mom[num][ik] = 0.0;
						for (int ikld = 0; ikld < size; ikld++)
						{
							Vyy_mom[num][ik] = Vyy_mom[num][ik] + all_mom[ik + ikld * 54];
						}
					}
				}
				///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				for (int ikld = 0; ikld < 54; ikld++)
				{
					one_mom[ikld] = Vxy_mom[num][ikld];
				}
				MPI_Gather(one_mom, 54, MPI_DOUBLE, all_mom,
					54, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				if (rank == 0)
				{
					for (int ik = 0; ik < 54; ik++)
					{
						Vxy_mom[num][ik] = 0.0;
						for (int ikld = 0; ikld < size; ikld++)
						{
							Vxy_mom[num][ik] = Vxy_mom[num][ik] + all_mom[ik + ikld * 54];
						}
					}
				}
				///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				for (int ikld = 0; ikld < 54; ikld++)
				{
					one_mom[ikld] = Vxxx_mom[num][ikld];
				}
				MPI_Gather(one_mom, 54, MPI_DOUBLE, all_mom,
					54, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				if (rank == 0)
				{
					for (int ik = 0; ik < 54; ik++)
					{
						Vxxx_mom[num][ik] = 0.0;
						for (int ikld = 0; ikld < size; ikld++)
						{
							Vxxx_mom[num][ik] = Vxxx_mom[num][ik] + all_mom[ik + ikld * 54];
						}
					}
				}
				///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				for (int ikld = 0; ikld < 54; ikld++)
				{
					one_mom[ikld] = T_mom[num][ikld];
				}
				MPI_Gather(one_mom, 54, MPI_DOUBLE, all_mom,
					54, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				if (rank == 0)
				{
					for (int ik = 0; ik < 54; ik++)
					{
						T_mom[num][ik] = 0.0;
						for (int ikld = 0; ikld < size; ikld++)
						{
							T_mom[num][ik] = T_mom[num][ik] + all_mom[ik + ikld * 54];
						}
					}
				}
				///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


			}
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			for (int ikld = 0; ikld < this->All_Cells.size(); ikld++)
			{
				auto k = this->All_Cells[ikld];
				my_send[ikld] = k->par[0].I_u;
			}

			MPI_Gather(my_send, this->All_Cells.size(), MPI_DOUBLE, all_,
				this->All_Cells.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

			if (rank == 0)
			{
				for (int ik = 0; ik < this->All_Cells.size(); ik++)
				{
					auto k = this->All_Cells[ik];
					k->par[0].I_u = 0.0;
					for (int ikld = 0; ikld < size; ikld++)
					{
						k->par[0].I_u = k->par[0].I_u + all_[ik + ikld * this->All_Cells.size()];
					}
				}
			}
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			for (int ikld = 0; ikld < this->All_Cells.size(); ikld++)
			{
				auto k = this->All_Cells[ikld];
				my_send[ikld] = k->par[0].I_v;
			}

			MPI_Gather(my_send, this->All_Cells.size(), MPI_DOUBLE, all_,
				this->All_Cells.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

			if (rank == 0)
			{
				for (int ik = 0; ik < this->All_Cells.size(); ik++)
				{
					auto k = this->All_Cells[ik];
					k->par[0].I_v = 0.0;
					for (int ikld = 0; ikld < size; ikld++)
					{
						k->par[0].I_v = k->par[0].I_v + all_[ik + ikld * this->All_Cells.size()];
					}
				}
			}
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			for (int ikld = 0; ikld < this->All_Cells.size(); ikld++)
			{
				auto k = this->All_Cells[ikld];
				my_send[ikld] = k->par[0].I_T;
			}

			MPI_Gather(my_send, this->All_Cells.size(), MPI_DOUBLE, all_,
				this->All_Cells.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

			if (rank == 0)
			{
				for (int ik = 0; ik < this->All_Cells.size(); ik++)
				{
					auto k = this->All_Cells[ik];
					k->par[0].I_T = 0.0;
					for (int ikld = 0; ikld < size; ikld++)
					{
						k->par[0].I_T = k->par[0].I_T + all_[ik + ikld * this->All_Cells.size()];
					}
				}
			}
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			for (int ikld = 0; ikld < this->All_Cells.size(); ikld++)
			{
				auto k = this->All_Cells[ikld];
				my_send[ikld] = k->par[0].II_u;
			}

			MPI_Gather(my_send, this->All_Cells.size(), MPI_DOUBLE, all_,
				this->All_Cells.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

			if (rank == 0)
			{
				for (int ik = 0; ik < this->All_Cells.size(); ik++)
				{
					auto k = this->All_Cells[ik];
					k->par[0].II_u = 0.0;
					for (int ikld = 0; ikld < size; ikld++)
					{
						k->par[0].II_u = k->par[0].II_u + all_[ik + ikld * this->All_Cells.size()];
					}
				}
			}
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			for (int ikld = 0; ikld < this->All_Cells.size(); ikld++)
			{
				auto k = this->All_Cells[ikld];
				my_send[ikld] = k->par[0].II_v;
			}

			MPI_Gather(my_send, this->All_Cells.size(), MPI_DOUBLE, all_,
				this->All_Cells.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

			if (rank == 0)
			{
				for (int ik = 0; ik < this->All_Cells.size(); ik++)
				{
					auto k = this->All_Cells[ik];
					k->par[0].II_v = 0.0;
					for (int ikld = 0; ikld < size; ikld++)
					{
						k->par[0].II_v = k->par[0].II_v + all_[ik + ikld * this->All_Cells.size()];
					}
				}
			}
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			for (int ikld = 0; ikld < this->All_Cells.size(); ikld++)
			{
				auto k = this->All_Cells[ikld];
				my_send[ikld] = k->par[0].II_T;
			}

			MPI_Gather(my_send, this->All_Cells.size(), MPI_DOUBLE, all_,
				this->All_Cells.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

			if (rank == 0)
			{
				for (int ik = 0; ik < this->All_Cells.size(); ik++)
				{
					auto k = this->All_Cells[ik];
					k->par[0].II_T = 0.0;
					for (int ikld = 0; ikld < size; ikld++)
					{
						k->par[0].II_T = k->par[0].II_T + all_[ik + ikld * this->All_Cells.size()];
					}
				}
			}
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			for (int ikld = 0; ikld < this->All_Cells.size(); ikld++)
			{
				auto k = this->All_Cells[ikld];
				my_send[ikld] = k->par[0].F_n;
			}

			MPI_Gather(my_send, this->All_Cells.size(), MPI_DOUBLE, all_,
				this->All_Cells.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

			if (rank == 0)
			{
				for (int ik = 0; ik < this->All_Cells.size(); ik++)
				{
					auto k = this->All_Cells[ik];
					k->par[0].F_n = 0.0;
					for (int ikld = 0; ikld < size; ikld++)
					{
						k->par[0].F_n = k->par[0].F_n + all_[ik + ikld * this->All_Cells.size()];
					}
				}
			}
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			for (int ikld = 0; ikld < this->All_Cells.size(); ikld++)
			{
				auto k = this->All_Cells[ikld];
				my_send[ikld] = k->par[0].F_u;
			}

			MPI_Gather(my_send, this->All_Cells.size(), MPI_DOUBLE, all_,
				this->All_Cells.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

			if (rank == 0)
			{
				for (int ik = 0; ik < this->All_Cells.size(); ik++)
				{
					auto k = this->All_Cells[ik];
					k->par[0].F_u = 0.0;
					for (int ikld = 0; ikld < size; ikld++)
					{
						k->par[0].F_u = k->par[0].F_u + all_[ik + ikld * this->All_Cells.size()];
					}
				}
			}
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			for (int ikld = 0; ikld < this->All_Cells.size(); ikld++)
			{
				auto k = this->All_Cells[ikld];
				my_send[ikld] = k->par[0].F_v;
			}

			MPI_Gather(my_send, this->All_Cells.size(), MPI_DOUBLE, all_,
				this->All_Cells.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

			if (rank == 0)
			{
				for (int ik = 0; ik < this->All_Cells.size(); ik++)
				{
					auto k = this->All_Cells[ik];
					k->par[0].F_v = 0.0;
					for (int ikld = 0; ikld < size; ikld++)
					{
						k->par[0].F_v = k->par[0].F_v + all_[ik + ikld * this->All_Cells.size()];
					}
				}
			}
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			for (int ikld = 0; ikld < this->All_Cells.size(); ikld++)
			{
				auto k = this->All_Cells[ikld];
				my_send[ikld] = k->par[0].F_T;
			}

			MPI_Gather(my_send, this->All_Cells.size(), MPI_DOUBLE, all_,
				this->All_Cells.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

			if (rank == 0)
			{
				for (int ik = 0; ik < this->All_Cells.size(); ik++)
				{
					auto k = this->All_Cells[ik];
					k->par[0].F_T = 0.0;
					for (int ikld = 0; ikld < size; ikld++)
					{
						k->par[0].F_T = k->par[0].F_T + all_[ik + ikld * this->All_Cells.size()];
					}
				}
			}
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			for (int i = 0; i < pop_atom; i++)
			{
				for (int ikld = 0; ikld < this->All_Cells.size(); ikld++)
				{
					auto k = this->All_Cells[ikld];
					my_send[ikld] = k->par[0].H_n[i];
				}

				MPI_Gather(my_send, this->All_Cells.size(), MPI_DOUBLE, all_,
					this->All_Cells.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

				if (rank == 0)
				{
					for (int ik = 0; ik < this->All_Cells.size(); ik++)
					{
						auto k = this->All_Cells[ik];
						k->par[0].H_n[i] = 0.0;
						for (int ikld = 0; ikld < size; ikld++)
						{
							k->par[0].H_n[i] = k->par[0].H_n[i] + all_[ik + ikld * this->All_Cells.size()];
						}
					}
				}
				////////////////////////////////////////////////////////////////////////////////////////////////////////
				for (int ikld = 0; ikld < this->All_Cells.size(); ikld++)
				{
					auto k = this->All_Cells[ikld];
					my_send[ikld] = k->par[0].H_u[i];
				}

				MPI_Gather(my_send, this->All_Cells.size(), MPI_DOUBLE, all_,
					this->All_Cells.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

				if (rank == 0)
				{
					for (int ik = 0; ik < this->All_Cells.size(); ik++)
					{
						auto k = this->All_Cells[ik];
						k->par[0].H_u[i] = 0.0;
						for (int ikld = 0; ikld < size; ikld++)
						{
							k->par[0].H_u[i] = k->par[0].H_u[i] + all_[ik + ikld * this->All_Cells.size()];
						}
					}
				}
				////////////////////////////////////////////////////////////////////////////////////////////////////////
				for (int ikld = 0; ikld < this->All_Cells.size(); ikld++)
				{
					auto k = this->All_Cells[ikld];
					my_send[ikld] = k->par[0].H_v[i];
				}

				MPI_Gather(my_send, this->All_Cells.size(), MPI_DOUBLE, all_,
					this->All_Cells.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

				if (rank == 0)
				{
					for (int ik = 0; ik < this->All_Cells.size(); ik++)
					{
						auto k = this->All_Cells[ik];
						k->par[0].H_v[i] = 0.0;
						for (int ikld = 0; ikld < size; ikld++)
						{
							k->par[0].H_v[i] = k->par[0].H_v[i] + all_[ik + ikld * this->All_Cells.size()];
						}
					}
				}
				////////////////////////////////////////////////////////////////////////////////////////////////////////
				for (int ikld = 0; ikld < this->All_Cells.size(); ikld++)
				{
					auto k = this->All_Cells[ikld];
					my_send[ikld] = k->par[0].H_T[i];
				}

				MPI_Gather(my_send, this->All_Cells.size(), MPI_DOUBLE, all_,
					this->All_Cells.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

				if (rank == 0)
				{
					for (int ik = 0; ik < this->All_Cells.size(); ik++)
					{
						auto k = this->All_Cells[ik];
						k->par[0].H_T[i] = 0.0;
						for (int ikld = 0; ikld < size; ikld++)
						{
							k->par[0].H_T[i] = k->par[0].H_T[i] + all_[ik + ikld * this->All_Cells.size()];
						}
					}
				}
				////////////////////////////////////////////////////////////////////////////////////////////////////////
				for (int ikld = 0; ikld < this->All_Cells.size(); ikld++)
				{
					auto k = this->All_Cells[ikld];
					my_send[ikld] = k->par[0].H_uu[i];
				}

				MPI_Gather(my_send, this->All_Cells.size(), MPI_DOUBLE, all_,
					this->All_Cells.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

				if (rank == 0)
				{
					for (int ik = 0; ik < this->All_Cells.size(); ik++)
					{
						auto k = this->All_Cells[ik];
						k->par[0].H_uu[i] = 0.0;
						for (int ikld = 0; ikld < size; ikld++)
						{
							k->par[0].H_uu[i] = k->par[0].H_uu[i] + all_[ik + ikld * this->All_Cells.size()];
						}
					}
				}
				////////////////////////////////////////////////////////////////////////////////////////////////////////
				for (int ikld = 0; ikld < this->All_Cells.size(); ikld++)
				{
					auto k = this->All_Cells[ikld];
					my_send[ikld] = k->par[0].H_vv[i];
				}

				MPI_Gather(my_send, this->All_Cells.size(), MPI_DOUBLE, all_,
					this->All_Cells.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

				if (rank == 0)
				{
					for (int ik = 0; ik < this->All_Cells.size(); ik++)
					{
						auto k = this->All_Cells[ik];
						k->par[0].H_vv[i] = 0.0;
						for (int ikld = 0; ikld < size; ikld++)
						{
							k->par[0].H_vv[i] = k->par[0].H_vv[i] + all_[ik + ikld * this->All_Cells.size()];
						}
					}
				}
				////////////////////////////////////////////////////////////////////////////////////////////////////////
				for (int ikld = 0; ikld < this->All_Cells.size(); ikld++)
				{
					auto k = this->All_Cells[ikld];
					my_send[ikld] = k->par[0].H_uv[i];
				}

				MPI_Gather(my_send, this->All_Cells.size(), MPI_DOUBLE, all_,
					this->All_Cells.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

				if (rank == 0)
				{
					for (int ik = 0; ik < this->All_Cells.size(); ik++)
					{
						auto k = this->All_Cells[ik];
						k->par[0].H_uv[i] = 0.0;
						for (int ikld = 0; ikld < size; ikld++)
						{
							k->par[0].H_uv[i] = k->par[0].H_uv[i] + all_[ik + ikld * this->All_Cells.size()];
						}
					}
				}
				////////////////////////////////////////////////////////////////////////////////////////////////////////
				for (int ikld = 0; ikld < this->All_Cells.size(); ikld++)
				{
					auto k = this->All_Cells[ikld];
					my_send[ikld] = k->par[0].H_uuu[i];
				}

				MPI_Gather(my_send, this->All_Cells.size(), MPI_DOUBLE, all_,
					this->All_Cells.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

				if (rank == 0)
				{
					for (int ik = 0; ik < this->All_Cells.size(); ik++)
					{
						auto k = this->All_Cells[ik];
						k->par[0].H_uuu[i] = 0.0;
						for (int ikld = 0; ikld < size; ikld++)
						{
							k->par[0].H_uuu[i] = k->par[0].H_uuu[i] + all_[ik + ikld * this->All_Cells.size()];
						}
					}
				}
			}


		}

		MPI_Barrier(MPI_COMM_WORLD);
		cout << "END peresilka " << rank << endl;
		MPI_Barrier(MPI_COMM_WORLD);
		
		for (int ikld = 0; ikld < this->All_Cells.size(); ikld++)
		{
			auto k = this->All_Cells[ikld];

			// Соберём данные со всех компов
			if (true)
			{
				if (k->df_s4_bool == true)
				{
					for (int i = 0; i < k->df_s4->n1; i++)
					{
						for (int j = 0; j < k->df_s4->n2; j++)
						{
							for (int jj = 0; jj < k->df_s4->n3; jj++)
							{
								ss_d[0] = k->df_s4->V[i][j][jj];
								MPI_Gather(ss_d, 1, MPI_DOUBLE, one_double_all,
									1, MPI_DOUBLE, 0, MPI_COMM_WORLD);         // Собираем число a со всех процессов в Главный
								MPI_Barrier(MPI_COMM_WORLD);
								if (rank == 0)
								{
									k->df_s4->V[i][j][jj] = 0.0;
									for (int ii = 0; ii < size; ii++)
									{
										k->df_s4->V[i][j][jj] = k->df_s4->V[i][j][jj] + one_double_all[ii];
									}
								}
							}
						}
					}
				}

				if (k->df_s3_bool == true)
				{
					for (int i = 0; i < k->df_s3->n1; i++)
					{
						for (int j = 0; j < k->df_s3->n2; j++)
						{
							for (int jj = 0; jj < k->df_s3->n3; jj++)
							{
								ss_d[0] = k->df_s3->V[i][j][jj];
								MPI_Gather(ss_d, 1, MPI_DOUBLE, one_double_all,
									1, MPI_DOUBLE, 0, MPI_COMM_WORLD);         // Собираем число a со всех процессов в Главный
								MPI_Barrier(MPI_COMM_WORLD);
								if (rank == 0)
								{
									k->df_s3->V[i][j][jj] = 0.0;
									for (int ii = 0; ii < size; ii++)
									{
										k->df_s3->V[i][j][jj] = k->df_s3->V[i][j][jj] + one_double_all[ii];
									}
								}
							}
						}
					}
				}
			}

		
			if (rank == 0)
			{
			//cout << "Hellow" <<  endl;
			double no = (1.0 * AllNumber * size * k->Get_Volume_rotate(360.0));
			//cout << "no = " << no << endl;

			k->par[0].I_u = sum_s * k->par[0].I_u / no;
			k->par[0].I_v = sum_s * k->par[0].I_v / no;
			k->par[0].I_T = sum_s * k->par[0].I_T / no;

			k->par[0].II_u = sum_s * k->par[0].II_u / no;
			k->par[0].II_v = sum_s * k->par[0].II_v / no;
			k->par[0].II_T = sum_s * k->par[0].II_T / no;



			if (k->par[0].F_n > 0.00000001)
			{
				k->par[0].F_u = k->par[0].F_u / k->par[0].F_n;
				k->par[0].F_v = k->par[0].F_v / k->par[0].F_n;
				k->par[0].F_T = (2.0 / 3.0) * (k->par[0].F_T / k->par[0].F_n - kvv(k->par[0].F_u, k->par[0].F_v, 0.0));
			}
			else
			{
				k->par[0].F_n = 0.0;
				k->par[0].F_u = 0.0;
				k->par[0].F_v = 0.0;
				k->par[0].F_T = 0.0;
			}



			if (k->par[0].F_n > 0.00000001)
			{
				for (int i = 0; i < pop_atom; i++)
				{
					if (k->par[0].H_n[i] > 0.00000001)
					{
						k->par[0].H_u[i] = k->par[0].H_u[i] / k->par[0].H_n[i];
						k->par[0].H_v[i] = k->par[0].H_v[i] / k->par[0].H_n[i];
						k->par[0].H_T[i] = (2.0 / 3.0) * (k->par[0].H_T[i] / k->par[0].H_n[i] - kvv(k->par[0].H_u[i], k->par[0].H_v[i], 0.0));
						k->par[0].H_uu[i] = k->par[0].H_uu[i] / k->par[0].H_n[i] - k->par[0].H_u[i] * k->par[0].H_u[i];
						k->par[0].H_vv[i] = k->par[0].H_vv[i] / k->par[0].H_n[i] - k->par[0].H_v[i] * k->par[0].H_v[i];
						k->par[0].H_uv[i] = k->par[0].H_uv[i] / k->par[0].H_n[i] - k->par[0].H_u[i] * k->par[0].H_v[i];
						k->par[0].H_uuu[i] = 2.0 * k->par[0].H_v[i] * k->par[0].H_v[i] * k->par[0].H_v[i] + //
							3.0 * k->par[0].H_v[i] * k->par[0].H_vv[i] - k->par[0].H_uuu[i] / k->par[0].H_n[i]; !ERROR
					}
					else
					{
						k->par[0].H_u[i] = 0.0;
						k->par[0].H_v[i] = 0.0;
						k->par[0].H_T[i] = 0.0;
						k->par[0].H_uu[i] = 0.0;
						k->par[0].H_vv[i] = 0.0;
						k->par[0].H_uv[i] = 0.0;
						k->par[0].H_uuu[i] = 0.0;
					}
				}
			}


			if (k->df_s4_bool == true)
			{

				if (rank == 0)
				{
					k->df_s4->normir(sum_s / no);
				}
			}

			if (k->df_s3_bool == true)
			{

				if (rank == 0)
				{
					k->df_s3->normir(sum_s / no);
				}
			}

			


				k->par[0].F_n = sum_s * k->par[0].F_n / no;
				for (int i = 0; i < pop_atom; i++)
				{
					k->par[0].H_n[i] = sum_s * k->par[0].H_n[i] / no;
				}


				k->par[0].I_u = (n_p_LISM_) * (k->par[0].I_u);
				k->par[0].I_v = (n_p_LISM_) * (k->par[0].I_v);
				k->par[0].I_T = (n_p_LISM_) * (k->par[0].I_T);
				k->par[0].II_u = (n_p_LISM_) * (k->par[0].II_u);
				k->par[0].II_v = (n_p_LISM_) * (k->par[0].II_v);
				k->par[0].II_T = (n_p_LISM_) * (k->par[0].II_T);

				// Считаем мультифдюидные источники

				double U_M_H[pop_atom];
				double U_H[pop_atom];
				double sigma_H[pop_atom];
				double nu_H[pop_atom];
				double q2_1, q2_2, q3;
				double u, v, ro, p, Q;

				double u_H[pop_atom];
				double v_H[pop_atom];
				double ro_H[pop_atom];
				double p_H[pop_atom];

				for (int i = 0; i < pop_atom; i++)
				{
					p_H[i] = 0.5 * k->par[0].H_T[i] * k->par[0].H_n[i];
					if (k->par[0].H_n[i] <= 0.0)
					{
						k->par[0].H_n[i] = 0.0000001;
						p_H[i] = 0.0;
					}
				}

				if (k->pui_ == true)
				{
					k->par[0].ro = k->par[0].ro + k->par[0].npui;
					double dd = k->par[0].p;
					k->par[0].p = k->par[0].pp;
					k->par[0].pp = dd;
				}

				u = k->par[0].u;
				v = k->par[0].v;
				ro = k->par[0].ro;
				p = k->par[0].p;

				if (ro <= 0.0)
				{
					ro = 0.0000001;
					p = 0.0;
				}


				for (int i = 0; i < pop_atom; i++)
				{
					U_M_H[i] = sqrt(kv(u - k->par[0].H_u[i]) + kv(v - k->par[0].H_v[i]) + (64.0 / (9.0 * pi_)) //
						* (p / ro + 2.0 * p_H[i] / k->par[0].H_n[i]));

					U_H[i] = sqrt(kv(u - k->par[0].H_u[i]) + kv(v - k->par[0].H_v[i]) + (4.0 / (pi_)) //
						* (p / ro + 2.0 * p_H[i] / k->par[0].H_n[i]));

					sigma_H[i] = kv(1.0 - a_2 * log(U_M_H[i]));

					nu_H[i] = ro * k->par[0].H_n[i] * U_M_H[i] * sigma_H[i];
				}

				k->par[0].M_u = 0.0;
				k->par[0].M_v = 0.0;
				k->par[0].M_T = 0.0;

				for (int i = 0; i < pop_atom; i++)
				{
					k->par[0].M_u = k->par[0].M_u + (n_p_LISM_ / Kn_) * nu_H[i] * (k->par[0].H_u[i] - u);
					k->par[0].M_v = k->par[0].M_v + (n_p_LISM_ / Kn_) * nu_H[i] * (k->par[0].H_v[i] - v);
					k->par[0].M_T = k->par[0].M_T + (n_p_LISM_ / Kn_) * nu_H[i] * ((kv(k->par[0].H_u[i]) + kv(k->par[0].H_v[i]) - kv(u) - kv(v)) / 2.0 + //
						(U_H[i] / U_M_H[i]) * (2.0 * p_H[i] / k->par[0].H_n[i] - p / ro));
				}


				if (fabs(k->par[0].M_u) > 0.000001)
				{
					k->par[0].k_u = (k->par[0].I_u) / k->par[0].M_u;
					if (k->par[0].k_u > 10.0 || k->par[0].k_u < 0.1)
					{
						k->par[0].k_u = 1.0;
					}
				}
				else
				{
					k->par[0].k_u = 1.0;
				}

				if (fabs(k->par[0].M_v) > 0.000001)
				{
					k->par[0].k_v = (k->par[0].I_v) / k->par[0].M_v;
					if (k->par[0].k_v > 10.0 || k->par[0].k_v < 0.1)
					{
						k->par[0].k_v = 1.0;
					}
				}
				else
				{
					k->par[0].k_v = 1.0;
				}

				if (fabs(k->par[0].M_T) > 0.000001)
				{
					k->par[0].k_T = (k->par[0].I_T) / k->par[0].M_T;
					if (k->par[0].k_T > 30.0 || k->par[0].k_T < 0.1)
					{
						k->par[0].k_T = 1.0;
					}
				}
				else
				{
					k->par[0].k_T = 1.0;
				}
			}

		}

		MPI_Barrier(MPI_COMM_WORLD);
		cout << "END rank = 0  ___  from   " << rank << endl;
		MPI_Barrier(MPI_COMM_WORLD);

		if (rank == 0)
		{
			free(all_);
		}
		free(my_send);

		if (rank == 0)  // Это нужно только для основного компа
		{
			for (auto& k : this->Cell_other)
			{
				k->par[0].F_n = 1.0;
				k->par[0].F_u = Velosity_inf;
				k->par[0].F_v = 0.0;
				k->par[0].F_T = 1.0;
				k->par[0].H_n[0] = 0.0;
				k->par[0].H_u[0] = 0.0;
				k->par[0].H_v[0] = 0.0;
				k->par[0].H_T[0] = 0.0;
				k->par[0].H_n[1] = 0.0;
				k->par[0].H_u[1] = 0.0;
				k->par[0].H_v[1] = 0.0;
				k->par[0].H_T[1] = 0.0;
				k->par[0].H_n[2] = 0.0;
				k->par[0].H_u[2] = 0.0;
				k->par[0].H_v[2] = 0.0;
				k->par[0].H_T[2] = 0.0;
				k->par[0].H_n[3] = 1.0;
				k->par[0].H_u[3] = Velosity_inf;
				k->par[0].H_v[3] = 0.0;
				k->par[0].H_T[3] = 1.0;

				k->par[0].I_u = 0.0;
				k->par[0].I_v = 0.0;
				k->par[0].I_T = 0.0;

				k->par[0].k_u = 0.0;
				k->par[0].k_v = 0.0;
				k->par[0].k_T = 0.0;
			}
		}


		if (func_stat && rank == 0)
		{
			double no;// = (1.0 * AllNumber * size);
			//fout2 << sum_s / no << endl;
			// (int)(the_ / (pi_ / Al_stat))
			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < Al_stat; j++)
				{
					no = (1.0 * AllNumber * size * kv(R_stat/ RR_) * 2.0 * pi_ * (cos(j * (pi_ / Al_stat)) - cos((j + 1.0) * (pi_ / Al_stat))) );

					fout2 << sum_s / no << " " <<  this->mu_mom[i][j] << " " << this->Vx_mom[i][j] << " " << this->Vy_mom[i][j] //
						<< " " << this->Vxx_mom[i][j] << " " << this->Vyy_mom[i][j] << " " << this->Vxy_mom[i][j] //
						<< " " << this->Vxxx_mom[i][j] << " " << this->T_mom[i][j] << endl;
				}
			}

			fout2.close();
		}

		MPI_Barrier(MPI_COMM_WORLD);
		cout << "END MPI_MK_start " << rank << endl;
		MPI_Barrier(MPI_COMM_WORLD);
#endif
}

void Setka::Fly_exchenge(Sensor* sens, double x_0, double y_0, double z_0, double Vx, double Vy, double Vz, Cell* now, double mu,//
	const double& mu_0, bool ExCh)
{
	//cout << "Setka.cpp    " << "NEW" << endl;
	//cout << "Setka.cpp    "  << x_0 << " " << sqrt(kvv(y_0, z_0, 0.0)) << " " << mu << endl;
	//double f1, f2;
	//now->Get_Center(f1, f2);
	//cout << "Setka.cpp    " << "Cent = " << f1 << " " << f2 << endl;
	double normV = sqrt(kvv(Vx, Vy, Vz));
	double dt = geo_accur * now->L / normV;
	//cout << "Setka.cpp    " << "dt = " << dt << endl;


	// Нахождение времени до попадания в следующую ячейку
	//cout << "Setka.cpp    " << "V = " << Vx << " " << Vy << " " << Vz << endl;
	while (now->belong(x_0 + 1000.0 * dt * Vx, sqrt(kv(y_0 + 1000.0 * dt * Vy) + kv(z_0 + 1000.0 * dt * Vz))) == false) // Слишком большой шаг
	{
		dt = dt / 2.0;
		//cout << "Setka.cpp    " << "dt = " << dt << endl;
	}

	int a = 1000;
	int b = 10000;
	while (now->belong(x_0 + (1.0 * b) * dt * Vx, sqrt(kv(y_0 + (1.0 * b) * dt * Vy) + kv(z_0 + (1.0 * b) * dt * Vz))) == true)
	{
		b = b + 500;
		//cout << "Setka.cpp    " << "b = " <<  x_0 + (1.0 * b) * dt * Vx << " " << sqrt(kv(y_0 + (1.0 * b) * dt * Vy) + kv(z_0 + (1.0 * b) * dt * Vz)) << " " << b <<  endl;
	}

	int k = (a + b) / 2;
	while (k != a)
	{
		if (now->belong(x_0 + (1.0 * k) * dt * Vx, sqrt(kv(y_0 + (1.0 * k) * dt * Vy) + kv(z_0 + (1.0 * k) * dt * Vz))) == true)
		{
			a = k;
		}
		else
		{
			b = k;
		}
		k = (a + b) / 2;
		//cout << "Setka.cpp    " << a << " " << b << endl;
	}

	if (b != a + 1)  // Можно потом удалить проверку
	{
		cout << "Setka.cpp    " << "ERRORIUEHUYCEUVDCC W" << endl;
		exit(-1);
	}

	double time = dt * (b); // Время нахождения в ячейке

	double uz, uz_M, uz_E, t1;// y, z, r;
	double cp = sqrt(now->par[0].p / now->par[0].ro);
	bool change_was = false;
	double vx = now->par[0].u;
	double vy = now->par[0].v;
	double ro = now->par[0].ro;

	double t_ex = 0.0;
	double mu2 = mu;
	double nu_ex, u, alpha, mu_ex;
	double t2 = time;

	double x_ex = x_0, y_ex = y_0, z_ex = z_0;
	double u1, u2, u3;

	if (ExCh == false) // Если перезарядки в этой ячейке ещё не было
	{
		alpha = polar_angle(y_0, z_0);
		u = sqrt(kvv(Vx - vx, Vy - vy * cos(alpha), Vz - vy * sin(alpha)));
		uz = Velosity_1(u, cp);
		nu_ex = (ro * uz * sigma2(uz, cp));
		double kappa = (nu_ex * time) / Kn_;
		t_ex = -(time / kappa) * log(1.0 - sens->MakeRandom() * (1.0 - exp(-kappa))); // Время до перезарядки
		t2 = t2 - t_ex;
		mu_ex = mu * (1.0 - exp(-kappa));
		mu2 = mu - mu_ex;
		alpha = polar_angle(y_0 + 0.5 * (t_ex)*Vy, z_0 + 0.5 * (t_ex)*Vz);


		x_ex = x_0 + t_ex * Vx;  // Координаты перезарядки
		y_ex = y_0 + t_ex * Vy;
		z_ex = z_0 + t_ex * Vz;
		if (now->belong(x_ex, sqrt(kv(y_ex) + kv(z_ex))) == false)
		{
			//cout << "Setka.cpp    " << "TUT" << endl;
			x_ex = x_0 + (1.0 * a * dt) * Vx;  // Координаты перезарядки
			y_ex = y_0 + (1.0 * a * dt) * Vy;
			z_ex = z_0 + (1.0 * a * dt) * Vz;
			t_ex = (1.0 * a * dt);
			t2 = time - t_ex;
			alpha = polar_angle(y_0 + 0.5 * (t_ex)*Vy, z_0 + 0.5 * (t_ex)*Vz);
		}

		// Считаем источники для частицы до перезарядки
		u = sqrt(kvv(Vx - vx, Vy - vy * cos(alpha), Vz - vy * sin(alpha)));
		uz = Velosity_1(u, cp);
		uz_M = Velosity_2(u, cp) / (uz * kv(cp) * cp * pi * sqrtpi_);
		uz_E = Velosity_3(u, cp);
		u1 = vx - Vx;
		u2 = vy * cos(alpha) - Vy;
		u3 = vy * sin(alpha) - Vz;
		double skalar = Vx * u1 + Vy * u2 + Vz * u3;


		now->mut.lock();
		now->par[0].F_n += t_ex * mu;
		now->par[0].F_u += t_ex * Vx * mu;
		now->par[0].F_v += t_ex * (Vy * cos(alpha) + Vz * sin(alpha)) * mu;
		now->par[0].F_T += t_ex * kvv(Vx, Vy, Vz) * mu;
		now->par[0].I_u += mu * t_ex * uz_M * uz * sigma(uz_M) * u1 / u;
		now->par[0].I_v += mu * t_ex * uz_M * uz * sigma(uz_M) * (u2 * cos(alpha) + u3 * sin(alpha)) / u;
		now->par[0].I_T += mu * t_ex * (0.5 * (3.0 * kv(cp) + 2.0 * kv(u)) * uz_E * sigma(uz_E) + 2.0 * uz_M * uz * sigma(uz_M) * skalar / u);
		now->mut.unlock();

		bool bb = false;
		if (mu_ex >= weight_ * mu_0)
		{
			//cout << "Setka.cpp    " << mu_ex << " podhodit" << endl;
			bb = true;
		}
		else
		{
			if (mu_ex >= sens->MakeRandom() * weight_ * mu_0)
			{
				mu_ex = weight_ * mu_0;
				//cout << "Setka.cpp    " << mu_ex << " ne podhodit = 0.25" << endl;
				bb = true;
			}
			else
			{
				bb = false;
				//cout << "Setka.cpp    " << mu_ex << " zabili" << endl;
			}
		}

		if (bb)  // Иначе нужно забить на эти траектории
		{

			double alpha = polar_angle(y_ex, z_ex);

			double Vx2, Vy2, Vz2;

			Change_Velosity(sens, vx / cp, vy * cos(alpha) / cp, vy * sin(alpha) / cp,//
				Vx / cp, Vy / cp, Vz / cp, Vx2, Vy2, Vz2, cp);  // здесь подаются  r, theta, phi
			Vx2 = Vx2 * cp;
			Vy2 = Vy2 * cp;
			Vz2 = Vz2 * cp;

			//cout << "Setka.cpp    " << "Change  " << endl;
			Fly_exchenge(sens, x_ex, y_ex, z_ex, Vx2, Vy2, Vz2, now, mu_ex, mu_0, true);
			//cout << "Setka.cpp    " << "Vozvrat change" << endl;
		}
	}

	//cout << "Setka.cpp    " << mu2 << " mu2" << endl;
	if (mu2 < weight_ * mu_0 && ExCh == false)
	{
		if (mu2 >= sens->MakeRandom() * weight_ * mu_0)
		{
			//cout << "Setka.cpp    " << mu2 << " mu2 = 0.25" << endl;
			mu2 = weight_ * mu_0;
		}
		else
		{
			//cout << "Setka.cpp    " << mu2 << " zabili" << endl;
			return;
		}
	}

	alpha = polar_angle(y_ex + 0.5 * t2 * Vy, z_ex + 0.5 * t2 * Vz);
	u = sqrt(kvv(Vx - vx, Vy - vy * cos(alpha), Vz - vy * sin(alpha)));
	uz = Velosity_1(u, cp);
	uz_M = Velosity_2(u, cp) / (uz * kv(cp) * cp * pi * sqrtpi_);
	uz_E = Velosity_3(u, cp);
	u1 = vx - Vx;
	u2 = vy * cos(alpha) - Vy;
	u3 = vy * sin(alpha) - Vz;
	double skalar = Vx * u1 + Vy * u2 + Vz * u3;

	now->mut.lock();
	now->par[0].F_n += t2 * mu2;
	now->par[0].F_u += t2 * Vx * mu2;
	now->par[0].F_v += t2 * (Vy * cos(alpha) + Vz * sin(alpha)) * mu2;
	now->par[0].F_T += t2 * kvv(Vx, Vy, Vz) * mu2;
	now->par[0].I_u += mu2 * t2 * uz_M * uz * sigma(uz_M) * u1 / u;
	now->par[0].I_v += mu2 * t2 * uz_M * uz * sigma(uz_M) * (u2 * cos(alpha) + u3 * sin(alpha)) / u;
	now->par[0].I_T += mu2 * t2 * (0.5 * (3.0 * kv(cp) + 2.0 * kv(u)) * uz_E * sigma(uz_E) + 2.0 * uz_M * uz * sigma(uz_M) * skalar / u);
	now->mut.unlock();

	double X = x_ex + t2 * Vx;
	double Y = y_ex + t2 * Vy;
	double Z = z_ex + t2 * Vz;


	Cell* next = nullptr;             // В какую ячейку попадаем следующую
	double xk, yk;
	xk = x_0 + (1.0 * b) * dt * Vx;
	yk = sqrt(kv(y_0 + (1.0 * b) * dt * Vy) + kv(z_0 + (1.0 * b) * dt * Vz));
	if (xk < Left_ || (sqrt(kv(xk) + kv(yk)) > R5_ && xk > 0) || (xk < 0 && yk > R5_))
	{
		next = nullptr;
	}
	else
	{
		for (auto& i : now->Grans)
		{
			if (i->Sosed != nullptr)
			{
				if (i->Sosed->belong(xk, yk) == true)
				{
					next = i->Sosed;
				}
			}
		}
		if (next == nullptr)
		{
			for (auto& k : now->Grans)
			{
				if (k->Sosed != nullptr)
				{
					for (auto& i : k->Sosed->Grans)
					{
						if (i->Sosed != nullptr)
						{
							if (i->Sosed->belong(xk, yk) == true)
							{
								next = i->Sosed;
							}
						}
					}
				}
			}
		}
		/*if (next == nullptr)
		{
			cout << "Setka.cpp    " << xk << " " << yk << " " << "no sosed" << endl;
		}*/
	}

	if (next != nullptr)
	{
		//cout << "Setka.cpp    " << "Next cell" << endl;
		Fly_exchenge(sens, X, Y, Z, Vx, Vy, Vz, next, mu2, mu_0, false);
		//cout << "Setka.cpp    " << "Vozvrat next" << endl;
	}
	return;

	// Попытка сделать честно, удалять не стал пока
	//for (auto& i : now->Grans)
	//{
	//	a = i->a;
	//	b = i->b;
	//	A = -2.0 * a * b * Vx - 2.0 * kv(a) * Vx * x_0 + 2.0 * Vy * y_0 + 2.0 * Vz * z_0;
	//	C = 2.0 * (kv(a) * kv(Vx) - kv(Vy) - kv(Vz));
	//	B = kv(2.0 * a * b * Vx + 2.0 * kv(a) * Vx * x_0 - 2.0 * Vy * y_0 - 2.0 * Vz * z_0) - 4.0 * (kv(a * Vx) - kv(Vy) - kv(Vz)) * //
	//		(kv(b) + 2.0 * a * b * x_0 + kv(a) * kv(x_0) - kv(y_0) - kv(z_0));
	//	if (fabs(C) > 0.000001)
	//	{
	//		if (B >= 0)
	//		{
	//			t1 = (A - sqrt(B)) / C;
	//			if (t1 > 0 && i->belong_gran(x_0 + t1 * Vx, sqrt(kv(y_0 + t1 * Vy) + kv(z_0 + t1 * Vz))))
	//			{
	//				if (t1 < t)
	//				{
	//					t = t1;

	//				}
	//			}
	//			t2 = (A + sqrt(B)) / C;
	//		}
	//	}
	//}
}

void Setka::Fly_exchenge_Imit_Korol(MKmethod& MK, Sensor* sens, double x_0, double y_0, double z_0, double Vx, double Vy, double Vz, Cell* now, double mu, //
	int area, bool ExCh, const double& mu_start, int to_I, int to_J, bool georaschep, int zon_stat)
	// Новая версия полёта нейтрального атома. Функция без лишнего кода и насыщена комментариями
	// Функция использует и геометрическое и физическое расщепление!!!
	// Описание переменных:
	// sens - датчик случайных чисел
	// x_0 ... -  координаты начала полёта нейтрального атома в данной ячейке
	// Vx ... -   скорость летящего нейтрального атома
	// now -      текущая ячейка
	// mu -       вес нейтрального атома
	// area -     какой атом летит
	// mu_start - начальный вес атома при инициализации (чтобы понимать, сколько от исходного осталось)
	// to_I - в какую область предназначался атом
{
	Cell* next_mb = nullptr;
	Sensor& adv = *sens;
	// Алгоритм вырубания лишних траектори (траектории вырубаются при пересечении новой зоны, далее 
	// вес лмбо увеличивается, либо частица уничтожается)
	// Не уверен, что to_I бывает больше, чем now->zona, кажется не должно такого быть

	if (false)//(area > 1 && to_I < now->zona)// (true && to_I < now->zona)
	{
		if (x_0 * Vx + y_0 * Vy + z_0 * Vz > 0.0 && to_I < now->zona) // Если летит вне солнца
		{
			double all = polar_angle(x_0, sqrt(kv(y_0) + kv(z_0)));
			to_I = now->zona;
			if (fabs(mu) < Mu[area][now->zona] * mu_start * 0.3 * sin(all))
			{
				if (fabs(mu) >= adv.MakeRandom() * Mu[area][now->zona] * mu_start * 0.3 * sin(all))
				{
					mu = Mu[area][now->zona] * mu_start * sign(mu) * 0.3 * sin(all);
				}
				else
				{
					return;
				}
			}
		}
	}

	/*if (now->number == 34)
	{
		if (fabs(mu) > Mu[area][now->zona] * mu_start * 0.03)
		{
			this->mut_1.lock();
			this->k_1++;
			this->mut_1.unlock();
		}
	}*/

	double normV = sqrt(kvv(Vx, Vy, Vz));
	double dt = now->L / normV;  // Время на один мини-шаг в ячейке (планируется сделать много шагов в ячейке)
	double ddt = 0.001 * dt;
	// Нахождение времени до попадания в следующую ячейку
	//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "V = " << Vx << " " << Vy << " " << Vz << endl;
	//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "A" << endl;
	double a = 0.0, b = dt;
	double t_per = 100000000000.0;
	// Блок проверки перигелия к оси
	/*double t_per = (-Vy * y_0 - Vz * z_0) / (kv(Vy) + kv(Vz));
	if (t_per > 0.0000001)
	{
		if (now->belong(x_0 + t_per * Vx, sqrt(kv(y_0 + t_per * Vy) + kv(z_0 + t_per * Vz))) == true)
		{
			t_per = -1.0;
		}
	}
	else
	{
		t_per = -1.0;
	}*/

	//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "A" << endl;
	double tt1, tt2;
	double a1, a2, a3;

	// Пытаемся определить время выхода точно
	if (true)
	{
		double A, B, C;
		for (auto& i : now->Grans)
		{
			if (i->type == Axis)
			{
				continue;
			}
			A = i->aa;
			B = i->bb;
			C = i->cc;

			/*if (fabs(B) < 0.001)
			{
				tt1 = (-C / A - x_0) / Vx;
				if (std::fpclassify(tt1) == FP_NORMAL && tt1 > 0.0)
				{
					if (t_per > tt1)
					{
						if (now->belong(x_0 + (tt1 + min(ddt, 0.03 * tt1)) * Vx, sqrt(kv(y_0 + (tt1 + min(ddt, 0.03 * tt1)) * Vy) + kv(z_0 + (tt1 + min(ddt, 0.03 * tt1)) * Vz))) == false)
						{
							t_per = tt1;
							next_mb = i->Sosed;
						}
					}
				}
				continue;
			}*/
			/*cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "1 - " << A << " " << B << " " << C << endl;
			cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "2 - " << x_0 << " " << y_0 << " " << z_0 << endl;
			cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "3 - " << Vx << " " << Vy << " " << Vz << endl;*/
			a1 = 4.0 * kv(A * Vx * (C + A * x_0) - kv(B) * (Vy * y_0 + Vz * z_0)) - 4.0 * (kv(A * Vx) - //
				B * B * (kv(Vy) + kv(Vz))) * (kv(C + A * x_0) - B * B * (kv(y_0) + kv(z_0)));
			a2 = -A * Vx * (C + A * x_0) + B * B * (Vy * y_0 + Vz * z_0);
			a3 = A * A * Vx * Vx - B * B * (Vy * Vy + Vz * Vz);
			//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "4 - " << a1 << " " << a2 << " " << a3 << endl;
			if (fabs(a1) < 0.0000001)
			{
				a1 = 0.0;
			}
			if (a1 >= 0.0)
			{
				tt1 = (a2 + 0.5 * sqrt(a1)) / a3;
				tt2 = (a2 - 0.5 * sqrt(a1)) / a3;
				//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< A << " " << B << " " << C << endl;
				if (std::fpclassify(tt1) == FP_NORMAL && tt1 > 0.0)
				{
					if (t_per > tt1)
					{
						if (now->belong(x_0 + (tt1 + min(ddt, 0.005 * tt1)) * Vx, sqrt(kv(y_0 + (tt1 + min(ddt, 0.005 * tt1)) * Vy) + kv(z_0 + (tt1 + min(ddt, 0.01 * tt1)) * Vz))) == false)
						{
							t_per = tt1;
							next_mb = i->Sosed;
						}
					}
				}
				if (std::fpclassify(tt2) == FP_NORMAL && tt2 > 0.0)
				{
					if (t_per > tt2)
					{
						if (now->belong(x_0 + (tt2 + min(ddt, 0.005 * tt2)) * Vx, sqrt(kv(y_0 + (tt2 + min(ddt, 0.005 * tt2)) * Vy) + kv(z_0 + (tt2 + min(ddt, 0.01 * tt2)) * Vz))) == false)
						{
							t_per = tt2;
							next_mb = i->Sosed;
						}
					}
				}
			}
		}

		a = t_per - min(ddt, 0.005 * t_per);
		b = t_per + min(ddt, 0.005 * t_per);
	}

	bool hand = false;
	if (now->belong(x_0 + a * Vx, sqrt(kv(y_0 + a * Vy) + kv(z_0 + a * Vz))) == false)
	{
		hand = true;
	}
	else if (now->belong(x_0 + b * Vx, sqrt(kv(y_0 + b * Vy) + kv(z_0 + b * Vz))) == true)
	{
		hand = true;
	}
	else if (t_per > 100000000.0)
	{
		//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "ERROR  8568  " << t_per << endl;
		hand = true;
	}

	double k;

	// Находим время  в ячейке
	int lk = 0;
	if (hand)
	{
		a = 0.0;
		b = dt;
		t_per = 100000000000.0;
		for (auto& i : now->Grans)
		{
			a1 = -kv(i->a * Vx) * (kv(i->a * Vx) - kv(Vy) - kv(Vz)) * kv(Vz * y_0 - Vy * z_0);
			a2 = (kv(i->a * Vx) - kv(Vy) - kv(Vz)) * (Vz * z_0 + Vy * y_0);
			a3 = (kv(Vy) + kv(Vz)) * (-kv(i->a * Vx) + kv(Vy) + kv(Vz));
			if (a1 >= 0)
			{
				tt1 = (sqrt(a1) + a2) / a3;
				tt2 = (-sqrt(a1) + a2) / a3;
				if (std::fpclassify(tt1) == FP_NORMAL && tt1 > 0.0)
				{
					if (now->belong(x_0 + tt1 * Vx, sqrt(kv(y_0 + tt1 * Vy) + kv(z_0 + tt1 * Vz))) == false)
					{
						if (t_per > tt1)
						{
							t_per = tt1;
						}
					}
				}
				if (std::fpclassify(tt2) == FP_NORMAL && tt2 > 0.0)
				{
					if (now->belong(x_0 + tt2 * Vx, sqrt(kv(y_0 + tt2 * Vy) + kv(z_0 + tt2 * Vz))) == false)
					{
						if (t_per > tt2)
						{
							t_per = tt2;
						}
					}
				}
			}
		}

		while (now->belong(x_0 + b * Vx, sqrt(kv(y_0 + b * Vy) + kv(z_0 + b * Vz))) == true) // Слишком маленький шаг
		{
			lk++;
			if (lk > 1000)
			{
				cout << "Setka.cpp    " << "Fly_ex_Korol   " << "8055 dt = " << dt << " " << x_0 << " " << sqrt(kv(y_0) + kv(z_0)) << endl;
				cout << "Setka.cpp    " << "Fly_ex_Korol   " << normV << " " << Vx << " " << Vy << " " << Vz << endl;
				exit(-1);
			}
			b = b * 1.1;
		}

		if (t_per < b)
		{
			b = t_per;
		}

		k = (a + b) / 2.0;
		while (fabs(a - b) * normV > now->L / 300.0)
		{
			if (now->belong(x_0 + k * Vx, sqrt(kv(y_0 + k * Vy) + kv(z_0 + k * Vz))) == true)
			{
				a = k;
			}
			else
			{
				b = k;
			}
			k = (a + b) / 2.0;
		}
	}

	double time = b; // Время нахождения в ячейке


	bool slay = false;
	double t1, t2;
	//Нужно проверить пересечение со сферой
	double D = kv(2.0 * Vx * x_0 + 2.0 * Vy * y_0 + 2.0 * Vz * z_0) - 4.0 * kvv(Vx, Vy, Vz) * (kvv(x_0, y_0, z_0) - kv(Rmax_));
	if (D >= 0.0)
	{
		t1 = (-2.0 * (Vx * x_0 + Vy * y_0 + Vz * z_0) - sqrt(D)) / (2.0 * kv(normV));
		t2 = (-2.0 * (Vx * x_0 + Vy * y_0 + Vz * z_0) + sqrt(D)) / (2.0 * kv(normV));
		if (t1 > 0.0 && t1 < time)
		{
			if (x_0 + Vx * t1 >= 0.0)
			{
				time = t1;
				a = time * 0.999;
				b = time * 1.001;
				slay = true;
			}
		}
		if (t2 > 0.0 && t2 < time)
		{
			if (x_0 + Vx * t2 >= 0.0)  // Для правой полусферы
			{
				time = t2;
				a = time * 0.999;
				b = time * 1.001;
				slay = true;
			}
		}
	}
	// Теперь пересечение с правой границей
	if (fabs(Vx) > 0.000001)
	{
		t1 = -x_0 / Vx;
		if (t1 > 0.0 && t1 <= time && kvv(0.0, y_0 + t1 * Vy, z_0 + t1 * Vz) >= kv(Rmax_))
		{
			time = t1;
			a = time * 0.999;
			b = time * 1.001;
			slay = true;
		}
	}

	//Мини - проверка
	//for (double kk = 0.0; kk < a ; kk += time / 50.0)
	//{
	//	if (now->belong(x_0 + kk * Vx, sqrt(kv(y_0 + kk * Vy) + kv(z_0 + kk * Vz))) == false)
	//	{
	//		cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "ERROR  8280 iruhfuerfefr  " << x_0 + kk * Vx << " " << sqrt(kv(y_0 + kk * Vy) + kv(z_0 + kk * Vz)) << " " << now->contour[0]->x << " " << now->contour[0]->y << " " << kk << " " << time << endl;
	//		cout << "Setka.cpp    " << "Fly_ex_Korol   "<< a << " " << b << endl;
	//		cout << "Setka.cpp    " << "Fly_ex_Korol   "<< x_0 << " " << sqrt(kv(y_0) + kv(z_0)) << endl;
	//		cout << "Setka.cpp    " << "Fly_ex_Korol   "<< x_0 + 0.01 * Vx << " " << sqrt(kv(y_0 + 0.01 * Vy) + kv(z_0 + 0.01 * Vz)) << endl;
	//		double xx, yy;
	//		now->Get_Center(xx, yy);
	//		cout << "Setka.cpp    " << "Fly_ex_Korol   "<< xx << " " << yy << endl;
	//		cout << "Setka.cpp    " << "Fly_ex_Korol   "<< " FFF " << endl;
	//		for (auto& i : now->Grans)
	//		{
	//			if (i->type == Axis)
	//			{
	//				continue;
	//			}
	//			double A = i->aa;
	//			double B = i->bb;
	//			double C = i->cc;
	//			cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "1 - " << A << " " << B << " " << C << endl;
	//			a1 = 4.0 * kv(A * Vx * (C + A * x_0) - kv(B) * (Vy * y_0 + Vz * z_0)) - 4.0 * (kv(A * Vx) - //
	//				B * B * (kv(Vy) + kv(Vz))) * (kv(C + A * x_0) - B * B * (kv(y_0) + kv(z_0)));
	//			a2 = -A * Vx * (C + A * x_0) + B * B * (Vy * y_0 + Vz * z_0);
	//			a3 = A * A * Vx * Vx - B * B * (Vy * Vy + Vz * Vz);
	//			cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "4 - " << a1 << " " << a2 << " " << a3 << endl;
	//			if (fabs(a1) < 0.000000001)
	//			{
	//				a1 = 0.0;
	//			}
	//			if (a1 >= 0.0)
	//			{
	//				tt1 = (a2 + 0.5 * sqrt(a1)) / a3;
	//				tt2 = (a2 - 0.5 * sqrt(a1)) / a3;
	//				cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "tt1 = " << tt1 << endl;
	//				cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "tt2 = " << tt2 << endl;
	//			}
	//		}
	//		exit(-1);
	//	}
	//}

	// теперь время time, a, b уже определены
	// a - последний момент времени в ячейке, b - первый момент за пределами ячейки


	//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< a << " " << b << endl;
	double uz, uz_M, uz_E;								// средние скорости в интегралах
	double cp = sqrt(now->par[0].p / now->par[0].ro);
	double vx = now->par[0].u;							// Скорости плазмы в ячейке
	double vy = now->par[0].v;
	double ro = now->par[0].ro;


	double t_ex = 0.0;									// время до перезарядки
	// time - время нахождения атома в ячейке
	t2 = time;									// время после перезарядки (будет ниже)

	double x_ex = x_0, y_ex = y_0, z_ex = z_0;		    // координаты перезарядки
	double u1, u2, u3;

	double y_start = sqrt(kv(y_0) + kv(z_0));

	double l = sqrt(kvv(time * Vx, time * Vy, time * Vz));                    // Расстояние до выхода из ячейки

	double alpha = 0.0;								    // угол 
	double u = 0.0;                                     // модуль относительной скорости атома и плазмы
	double nu_ex = 0.0;								    // частота перезарядки
	double mu_ex = 0.0;									// вес перезаряженного атома
	double mu2 = mu;									// вес не-перезаряженного атома

	double drob = 5.0;     // Для того чтобы точнее учесть проворот скорости во время полёта

	if (true)//(ExCh == false) // Если перезарядки в этой ячейке ещё не было
	{
		double kappa = 0.0;
		//drob = min(fabs(polar_angle(y_0, z_0) - polar_angle(y_0 + time * Vy, z_0 + time * Vz)) / 0.017 + 5.0, 80.0);

		//if (y_start < 40.0/RR_)
		//{
		//	drob = 30.0;  // 30.0
		//}
		//else if (y_start < 60.0 / RR_)
		//{
		//	drob = 20.0;  // 20.0
		//}
		//else if (y_start < 80.0 / RR_)
		//{
		//	drob = 15.0;
		//}
		//else if (y_start < 120.0 / RR_)
		//{
		//	drob = 10.0;
		//}


		if (true)
		{

			for (double tt = 0.0; tt < time * 0.98; tt += time / drob)
			{
				double yy = sqrt(kv(y_0 + (tt + time / (2.0 * drob)) * Vy) + kv(z_0 + (tt + time / (2.0 * drob)) * Vz));

				alpha = polar_angle(y_0 + (tt + time / (2.0 * drob)) * Vy, z_0 + (tt + time / (2.0 * drob)) * Vz);
				u = sqrt(kvv(Vx - vx, Vy - vy * cos(alpha), Vz - vy * sin(alpha)));


				if (u / cp > 7.0)
				{
					uz = Velosity_1(u, cp);
					nu_ex = (ro * uz * sigma(uz)) / Kn_;
					//cout << x_0 << " " << y_start << " " << u / cp << " " << u << " " << cp <<  endl;
				}
				else
				{
					nu_ex = (ro * MK.int_1(u, cp)) / Kn_;        // Пробуем вычислять интеграллы численно

					//uz = Velosity_1(u, cp);
					//cout << (uz * sigma(uz)) << " " << MK.int_1(u, cp) << " " << u << " " << cp << endl;
				}
				kappa += (nu_ex * time / drob);
			}
		}
		else
		{
			alpha = polar_angle(y_0, z_0);
			u = sqrt(kvv(Vx - vx, Vy - vy * cos(alpha), Vz - vy * sin(alpha)));
			if (u / cp > 7.0)
			{
				uz = Velosity_1(u, cp);
				nu_ex = (ro * uz * sigma(uz)) / Kn_;
			}
			else
			{
				nu_ex = (ro * MK.int_1(u, cp)) / Kn_;        // Пробуем вычислять интеграллы численно
			}
			kappa = (nu_ex * time);
		}


		t_ex = -(time / kappa) * log(1.0 - adv.MakeRandom() * (1.0 - exp(-kappa))); // Время до перезарядки
		t2 = time - t_ex;   // Время сколько лететь после того, как атом перезарядился
		mu_ex = mu * (1.0 - exp(-kappa)); // вес перезаряженного атома
		mu2 = mu - mu_ex;                 // вес оставшегося неперезаряженного атома

		x_ex = x_0 + t_ex * Vx;  // Координаты перезарядки
		y_ex = y_0 + t_ex * Vy;
		z_ex = z_0 + t_ex * Vz;

		alpha = polar_angle(y_ex, z_ex);

		// проверка, если перезарядка произошла за пределами ячейки (нужно немного сдвигать её в этом случае)
		if (now->belong(x_ex, sqrt(kv(y_ex) + kv(z_ex))) == false)
		{
			x_ex = x_0 + a * Vx;  // Координаты перезарядки
			y_ex = y_0 + a * Vy;
			z_ex = z_0 + a * Vz;
			t_ex = a;
			t2 = time - t_ex;
			alpha = polar_angle(y_ex, z_ex);
		}

		double r = sqrt(kvv(x_ex, y_ex, z_ex));
		double all = polar_angle(x_ex, sqrt(kv(y_ex) + kv(z_ex)));

		// для сбора статистики
		if (func_stat)
		{
			if (r < 100.0 / RR_ && r > 60.0 / RR_)
			{
				bool sec = false;
				bool sec2 = false;
				double tt;
				double a1 = Vx * x_0 + Vy * y_0 + Vz * z_0;
				double a2 = kvv(Vx, Vy, Vz);
				double a3 = kvv(x_0, y_0, z_0) - kv(R_stat / RR_);
				if (4.0 * a1 * a1 - 4.0 * a2 * a3 >= 0.0)
				{
					double t1 = (-2.0 * a1 + sqrt(4.0 * a1 * a1 - 4.0 * a2 * a3)) / (2.0 * a2);
					double t2 = (-2.0 * a1 - sqrt(4.0 * a1 * a1 - 4.0 * a2 * a3)) / (2.0 * a2);
					if (t1 > 0.000000001 && t1 < time)
					{
						sec = true;
						tt = t1;
						if (t1 < t_ex)
						{
							sec2 = true; // пролетел до перезярядки
						}
						else
						{
							sec2 = false;
						}
					}
					if (t2 > 0.000000001 && t2 < time)
					{
						sec = true;
						tt = t2;
						if (t2 < t_ex)
						{
							sec2 = true; // пролетел до перезярядки
						}
						else
						{
							sec2 = false;
						}
					}
				}

				if (sec == true)
				{
					double xx, yy, zz;
					xx = x_0 + tt * Vx;
					yy = y_0 + tt * Vy;
					zz = z_0 + tt * Vz;
					double Vr, Vphi, Vthe;
					double phi_ = polar_angle(yy, zz);  //polar_angle(xx, sqrt(kv(yy) + kv(zz)));
					spherical_skorost(yy, zz, xx, Vy, Vz, Vx, Vr, Vphi, Vthe);

					double Vn = fabs((Vx * xx + Vy * yy + Vz * zz) / (R_stat / RR_));
					//double Vn = sqrt(kvv(Vx, Vy, Vz));
					double Vrr = (Vy * cos(phi_) + Vz * sin(phi_));
					double mumu = 0.0;
					if (sec2)
					{
						mumu = mu / Vn;
					}
					else
					{
						mumu = mu2 / Vn;
					}

					mut_mom.lock();
					this->mu_mom[area][(int)(phi_ / (pi_ / Al_stat))] += mumu;
					this->Vx_mom[area][(int)(phi_ / (pi_ / Al_stat))] += mumu * Vx;
					this->Vy_mom[area][(int)(phi_ / (pi_ / Al_stat))] += mumu * Vrr;
					this->Vxx_mom[area][(int)(phi_ / (pi_ / Al_stat))] += mumu * Vx * Vx;
					this->Vyy_mom[area][(int)(phi_ / (pi_ / Al_stat))] += mumu * Vrr * Vrr;
					this->Vxy_mom[area][(int)(phi_ / (pi_ / Al_stat))] += mumu * Vx * Vrr;
					this->Vxxx_mom[area][(int)(phi_ / (pi_ / Al_stat))] += mumu * Vx * Vx * Vx;
					mut_mom.unlock();

					//if (fabs(mu) > pred(area, now->zona, all) * 0.3)
					//{
					//	if (phi_ < 0.17453293)  // Угол в 10 градусов
					//	{
					//		mut_stat.lock();
					//		this->V_r_stat[this->number_stat] = Vr;
					//		this->V_t_stat[this->number_stat] = Vthe;
					//		this->V_p_stat[this->number_stat] = Vphi;
					//		this->phi_stat[this->number_stat] = phi_;
					//		this->num_stat[this->number_stat] = area;
					//		if (sec2)
					//		{
					//			this->mu_stat[this->number_stat] = mu / Vn;
					//		}
					//		else
					//		{
					//			this->mu_stat[this->number_stat] = mu2 / Vn;
					//		}
					//		this->number_stat++;
					//		mut_stat.unlock();
					//	}
					//}
				}
			}
		}

		double rr = sqrt(kvv(0.0, y_ex, z_ex));


		// Считаем источники для частицы до перезарядки

		u = sqrt(kvv(Vx - vx, Vy - vy * cos(alpha), Vz - vy * sin(alpha)));
		u1 = vx - Vx;
		u2 = vy * cos(alpha) - Vy;
		u3 = vy * sin(alpha) - Vz;
		double skalar = Vx * u1 + Vy * u2 + Vz * u3;

		// определим область перезарядки (по геометрическим зонам)
		int i_z = 0, i_alp = 0;
		i_z = geo_zones(r);
		i_alp = alpha_zones(x_ex, rr);

		int area2 = 0;
		// определим область в которой находится атом сейчас (это параметр самой ячейки)
		if (now->type == C_5)
		{
			area2 = 3;
		}
		else if (now->type == C_4)
		{
			area2 = 2;
		}
		else if (now->type == C_1 || now->type == C_centr)
		{
			area2 = 0;
		}
		else
		{
			area2 = 1;
		}


		now->mut.lock();

		// Блок расчёта функции распределения, если это нужно
		if (area == 3)
		{
			if (now->df_s4_bool == true)
			{
				now->df_s4->Add_point(Vx, Vy, Vz, y_ex, z_ex, t_ex* mu);
			}
		}
		else if (area == 2)
		{
			if (now->df_s3_bool == true)
			{
				now->df_s3->Add_point(Vx, Vy, Vz, y_ex, z_ex, t_ex * mu);
			}
		}

		//now->par[0].w_m[area] += mu / max(sin(alpha), 0.3 * Sinus[i_alp]);
		now->par[0].F_n += t_ex * mu;
		now->par[0].F_u += t_ex * Vx * mu;
		now->par[0].F_T += t_ex * kvv(Vx, Vy, Vz) * mu;
		now->par[0].H_n[area] += t_ex * mu;
		now->par[0].H_u[area] += t_ex * Vx * mu;
		now->par[0].H_T[area] += t_ex * kvv(Vx, Vy, Vz) * mu;

		double Vrr = (Vy * cos(alpha) + Vz * sin(alpha));
		now->par[0].H_uu[area] += t_ex * Vx * Vx * mu;
		now->par[0].H_uv[area] += t_ex * Vrr * Vx * mu;
		now->par[0].H_vv[area] += t_ex * Vrr * Vrr * mu;
		now->par[0].H_uuu[area] += t_ex * Vx * Vx * Vx * mu;

		if (u / cp > 7.0)
		{
			uz = Velosity_1(u, cp);
			uz_M = Velosity_2(u, cp) / (uz * kv(cp) * cp * pi * sqrtpi_);
			uz_E = Velosity_3(u, cp);

			now->par[0].I_u += -mu_ex * uz_M * u1 / u;
			now->par[0].I_v += -mu_ex * uz_M * (u2 * cos(alpha) + u3 * sin(alpha)) / u;
			now->par[0].I_T += mu_ex * (-0.25 * (3.0 * kv(cp) + 2.0 * kv(u)) * (uz_E / uz) - //
				uz_M * skalar / u);
		}
		else
		{
			double k1 = MK.int_1(u, cp);
			double k2 = MK.int_2(u, cp);
			double k3 = MK.int_3(u, cp);
			now->par[0].I_u += mu_ex * (k2 / k1) * u1 / u;
			now->par[0].I_v += mu_ex * (k2 / k1) * (u2 * cos(alpha) + u3 * sin(alpha)) / u;
			now->par[0].I_T += mu_ex * (-0.5 * k3 / k1 + k2 / k1 * skalar / u);
		}

		//drob = min(fabs(polar_angle(y_0, z_0) - polar_angle(y_0 + t_ex * Vy, z_0 + t_ex * Vz)) / 0.017 + 5.0, 80.0);

		if (true)
		{
			for (double tt = 0.0; tt < t_ex * 0.99; tt += t_ex / drob)
			{
				alpha = polar_angle(y_0 + (tt + t_ex / (2.0 * drob)) * Vy, z_0 + (tt + t_ex / (2.0 * drob)) * Vz);

				now->par[0].F_v += (t_ex / drob) * (Vy * cos(alpha) + Vz * sin(alpha)) * mu;
				now->par[0].H_v[area] += (t_ex / drob) * (Vy * cos(alpha) + Vz * sin(alpha)) * mu;
			}
		}
		else
		{
			now->par[0].F_v += t_ex * (Vy * cos(alpha) + Vz * sin(alpha)) * mu;
			now->par[0].H_v[area] += t_ex * (Vy * cos(alpha) + Vz * sin(alpha)) * mu;
		}

		//now->par[0].I_u += -mu_ex * uz_M * (sigma2(uz_M, cp) / sigma2(uz, cp)) * u1 / u;
		//now->par[0].I_v += -mu_ex * uz_M * (sigma2(uz_M, cp) / sigma2(uz, cp)) * (u2 * cos(alpha) + u3 * sin(alpha)) / u;
		//now->par[0].I_T += mu_ex * (-0.25 * (3.0 * kv(cp) + 2.0 * kv(u)) * (uz_E / uz) * (sigma2(uz_E, cp) / sigma2(uz, cp)) - //
		//	uz_M * (sigma2(uz_M, cp) / sigma2(uz, cp)) * skalar / u);

		/*now->par[0].II_u += -mu_ex * u1;
		now->par[0].II_v += -mu_ex * (u2 * cos(alpha) + u3 * sin(alpha));
		now->par[0].II_T += mu_ex * 0.5 * (kvv(Vx, Vy, Vz) - kvv(vx, vy, 0.0));*/

		now->mut.unlock();


		double Ur, Uthe, Uphi, Vr, Vthe, Vphi;
		spherical_skorost(y_ex, z_ex, x_ex, vy * cos(alpha), vy * sin(alpha), vx, Ur, Uphi, Uthe);
		spherical_skorost(y_ex, z_ex, x_ex, Vy, Vz, Vx, Vr, Vphi, Vthe);
		double mu3;


		int stat_zone_ = now->zona;

		if (area != area2)
		{
			stat_zone_ = -1;
		}


		if (area2 == 0 || Ur / cp > 3.0)   // Без геометрического расщепления
		{
		aa:
			bool kj = true;
			mu3 = mu_ex;
			double alpha = polar_angle(y_ex, z_ex);
			double aa, bb, cc;  // Декартовы скорости после перезарядки
			vector <double> Wr(1);
			vector <double> Wthe(1);
			vector <double> Wphi(1);
			vector <double> mu_(1);
			MK.Change_Velosity4(sens, Ur / cp, Uthe / cp, Uphi / cp, Vr / cp, Vthe / cp, Vphi / cp, Wr, Wthe, Wphi, mu_, cp, r, 0);
			Wr[0] = Wr[0] * cp;
			Wthe[0] = Wthe[0] * cp;
			Wphi[0] = Wphi[0] * cp;
			dekard_skorost(y_ex, z_ex, x_ex, Wr[0], Wphi[0], Wthe[0], bb, cc, aa);
			double peregel, time_do_peregel;
			int ii_z;
			int ii_alp;

			if (Wr[0] < 0.0)
			{
				//wwt = sqrt(kv(Wphi[0]) + kv(Wthe[0]));
				//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "Letit vniz  " << x_ex * RR_ << " " << rr * RR_ << " " << mu3 << " " << area << " " << area2 << endl;
				//goto aa;
				time_do_peregel = (-aa * x_ex - bb * y_ex - cc * z_ex) / kvv(aa, bb, cc);
				peregel = sqrt(kvv(x_ex + aa * time_do_peregel, y_ex + bb * time_do_peregel, z_ex + cc * time_do_peregel));
				ii_z = this->geo_zones(peregel);                     // Номер зоны перегелия атома
				ii_alp = alpha_zones(x_ex + aa * time_do_peregel, sqrt(kvv(y_ex + bb * time_do_peregel, z_ex + cc * time_do_peregel, 0.0)));
				//kj = false;   // Убиваем траекторию, иначе может испортить статистику
			}
			else
			{
				ii_z = i_z;
				ii_alp = i_alp;
			}

			if (kj == true)
			{
				if (fabs(mu3) >= pred(area2, ii_z, all))
				{
					kj = true;
				}
				else
				{
					if (fabs(mu3) >= adv.MakeRandom() * pred(area2, ii_z, all))
					{
						mu3 = pred(area2, ii_z, all) * sign(mu3);
						kj = true;
					}
					else
					{
						kj = false;
					}
				}
			}

			if (kj == true)
			{
				now->mut.lock();
				now->par[0].II_u += mu3 * (Vx - aa);
				now->par[0].II_v += mu3 * ((Vy - bb) * cos(alpha) + (Vz - cc) * sin(alpha));
				now->par[0].II_T += mu3 * 0.5 * (kvv(Vx, Vy, Vz) - kvv(aa, bb, cc));
				now->mut.unlock();
				Fly_exchenge_Imit_Korol(MK, sens, x_ex, y_ex, z_ex, aa, bb, cc, now, mu3, area2, false, mu_start, ii_z, ii_alp, true, stat_zone_);
			}
			//this->Fly_exchenge_Imit(MK, sens, x_ex, y_ex, z_ex, aa, bb, cc, now, mu3, -log(1.0 - sens->MakeRandom()),//
			//	0.0, area2, mu_start, to_I, 0);

		}
		else  // делаем геометрическое расщепление
		{
			int I = geo_zones(r, 1.2);  // Число дополнительных траекторий
			if (I > to_I)  // Для того чтобы не расщеплять атом, который летел в to_I область на атомы из более крупных областей
			{
				I = to_I;
			}

			vector <double> Wr(I + 1);
			vector <double> Wthe(I + 1);
			vector <double> Wphi(I + 1);
			vector <double> mu_(I + 1);

			//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< " C " << endl;
			bool bbb = MK.Change_Velosity4(sens, Ur / cp, Uthe / cp, Uphi / cp, Vr / cp, Vthe / cp, Vphi / cp, Wr, Wthe, Wphi, mu_, cp, r, I);
			//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< " D " << endl;

			for (int i = 0; i <= I; i++)
			{
				Wr[i] = Wr[i] * cp;
				Wthe[i] = Wthe[i] * cp;
				Wphi[i] = Wphi[i] * cp;
			}

			double aa, bb, cc, time_do_peregel, peregel;
			int ii_z, ii_alp;
			bool kj = true;
			//now->par[0].num_atoms++;
			for (int i = 0; i < I; i++)
			{
				dekard_skorost(y_ex, z_ex, x_ex, Wr[i], Wphi[i], Wthe[i], bb, cc, aa);
				mu3 = mu_[i] * mu_ex;
				time_do_peregel = (-aa * x_ex - bb * y_ex - cc * z_ex) / kvv(aa, bb, cc);
				peregel = sqrt(kvv(x_ex + aa * time_do_peregel, y_ex + bb * time_do_peregel, z_ex + cc * time_do_peregel));
				ii_z = this->geo_zones(peregel);                     // Номер зоны перегелия атома
				ii_alp = alpha_zones(x_ex + aa * time_do_peregel, sqrt(kvv(y_ex + bb * time_do_peregel, z_ex + cc * time_do_peregel, 0.0)));

				if (fabs(mu3) >= pred(area2, ii_z, all))
				{
					kj = true;
				}
				else
				{
					if (fabs(mu3) >= adv.MakeRandom() * pred(area2, ii_z, all))
					{
						mu3 = pred(area2, ii_z, all);
						kj = true;
					}
					else
					{
						kj = false;
					}
				}
				if (kj == true)
				{
					now->mut.lock();
					now->par[0].II_u += mu3 * (Vx - aa);
					now->par[0].II_v += mu3 * ((Vy - bb) * cos(alpha) + (Vz - cc) * sin(alpha));
					now->par[0].II_T += mu3 * 0.5 * (kvv(Vx, Vy, Vz) - kvv(aa, bb, cc));
					now->mut.unlock();
					Fly_exchenge_Imit_Korol(MK, sens, x_ex, y_ex, z_ex, aa, bb, cc, now, mu3, area2, false, mu_start, ii_z, ii_alp, true, stat_zone_);
				}
			}

			if (bbb == true)  // Если вообще нужно запускать основной атом
			{
				mu3 = mu_[I] * mu_ex;

				// Нужно определить номер основного атома (от его перегелия)
				dekard_skorost(y_ex, z_ex, x_ex, Wr[I], Wphi[I], Wthe[I], bb, cc, aa);
				if (Wr[I] < 0.0)
				{
					time_do_peregel = (-aa * x_ex - bb * y_ex - cc * z_ex) / kvv(aa, bb, cc);
					peregel = sqrt(kvv(x_ex + aa * time_do_peregel, y_ex + bb * time_do_peregel, z_ex + cc * time_do_peregel));
					ii_z = this->geo_zones(peregel);                     // Номер зоны перегелия атома
					ii_alp = alpha_zones(x_ex + aa * time_do_peregel, sqrt(kvv(y_ex + bb * time_do_peregel, z_ex + cc * time_do_peregel, 0.0)));
				}
				else
				{
					ii_z = i_z;
					ii_alp = i_alp;
				}

				if (fabs(mu3) >= pred(area2, ii_z, all))
				{
					kj = true;
				}
				else
				{
					if (fabs(mu3) >= adv.MakeRandom() * pred(area2, ii_z, all))
					{
						mu3 = pred(area2, ii_z, all) * sign(mu3);
						kj = true;
					}
					else
					{
						kj = false;
					}
				}

				if (kj == true)
				{
					now->mut.lock();
					now->par[0].II_u += mu3 * (Vx - aa);
					now->par[0].II_v += mu3 * ((Vy - bb) * cos(alpha) + (Vz - cc) * sin(alpha));
					now->par[0].II_T += mu3 * 0.5 * (kvv(Vx, Vy, Vz) - kvv(aa, bb, cc));
					now->mut.unlock();
					Fly_exchenge_Imit_Korol(MK, sens, x_ex, y_ex, z_ex, aa, bb, cc, now, mu3, area2, true, mu_start, ii_z, ii_alp, false, stat_zone_);
				}
			}

		}


		// Здесь можно добавить рулетку для оставшегося неперезаряженного атома
		if (fabs(mu2) < pred(area, to_I, all))
		{
			if (fabs(mu2) >= adv.MakeRandom() * pred(area, to_I, all))
			{
				mu2 = pred(area, to_I, all) * sign(mu2);
			}
			else
			{
				return;
			}
		}

	}


	// этот блок кода будет работать и для частиц, которые не попали в предыдущий цикл
	// здесь t2, _ex будут иметь исходные значения, как при инициализации
	alpha = polar_angle(y_ex + 0.5 * t2 * Vy, z_ex + 0.5 * t2 * Vz);
	now->mut.lock();

	// Блок расчёта функции распределения, если это нужно
	if (area == 3)
	{
		if (now->df_s4_bool == true)
		{
			now->df_s4->Add_point(Vx, Vy, Vz, y_ex, z_ex, t2 * mu2);
		}
	}
	else if (area == 2)
	{
		if (now->df_s3_bool == true)
		{
			now->df_s3->Add_point(Vx, Vy, Vz, y_ex, z_ex, t2 * mu2);
		}
	}

	now->par[0].F_n += t2 * mu2;
	now->par[0].F_u += t2 * Vx * mu2;
	now->par[0].F_T += t2 * kvv(Vx, Vy, Vz) * mu2;
	now->par[0].H_n[area] += t2 * mu2;
	now->par[0].H_u[area] += t2 * Vx * mu2;
	now->par[0].H_T[area] += t2 * kvv(Vx, Vy, Vz) * mu2;

	double Vrr = (Vy * cos(alpha) + Vz * sin(alpha));
	now->par[0].H_uu[area] += t2 * Vx * Vx * mu2;
	now->par[0].H_uv[area] += t2 * Vrr * Vx * mu2;
	now->par[0].H_vv[area] += t2 * Vrr * Vrr * mu2;
	now->par[0].H_uuu[area] += t2 * Vx * Vx * Vx * mu2;

	//drob = min(fabs(polar_angle(y_ex, z_ex) - polar_angle(y_ex + t2 * Vy, z_ex + t2 * Vz)) / 0.017 + 5.0, 80.0);

	if (true)
	{
		for (double tt = 0.0; tt < t2 * 0.99; tt += t2 / drob)  // было 20 и норм работало
		{
			alpha = polar_angle(y_ex + (tt + t2 / (2.0 * drob)) * Vy, z_ex + (tt + t2 / (2.0 * drob)) * Vz);

			now->par[0].F_v += (t2 / drob) * (Vy * cos(alpha) + Vz * sin(alpha)) * mu2;
			now->par[0].H_v[area] += (t2 / drob) * (Vy * cos(alpha) + Vz * sin(alpha)) * mu2;
		}
	}
	else
	{
		now->par[0].F_v += t2 * (Vy * cos(alpha) + Vz * sin(alpha)) * mu2;
		now->par[0].H_v[area] += t2 * (Vy * cos(alpha) + Vz * sin(alpha)) * mu2;
	}
	now->mut.unlock();

	double mu3 = mu2;
	/*double X = x_ex + t2 * Vx;
	double Y = y_ex + t2 * Vy;
	double Z = z_ex + t2 * Vz;*/
	double X = x_0 + b * Vx;
	double Y = y_0 + b * Vy;
	double Z = z_0 + b * Vz;

	Cell* next = nullptr;             // В какую ячейку попадаем следующую
	double xk, yk;
	xk = x_0 + b * Vx;
	yk = sqrt(kv(y_0 + b * Vy) + kv(z_0 + b * Vz));
	double r_k = sqrt(kv(xk) + kv(yk));
	double r_0 = sqrt(kv(x_0) + kv(y_0) + kv(z_0));

	/*if (yk < 5.0)
	{
		cout << "Setka.cpp    " << "Fly_ex_Korol   "<< b << "  " << xk << endl;
	}*/

	if (slay == true)  //Если надо убить траекторию
	{
		return;
	}

	if (xk < Left_ + 2.0 / RR_)  // ЭТОГО НЕ БЫЛО
	{
		return;
	}

	if (yk > R5_ - 2.0 / RR_)  // ЭТОГО НЕ БЫЛО
	{
		return;
	}

	if (r_k >= Rmax_ - 3.0 / RR_ && xk >= -0.5 / RR_)
	{
		return;
		next = nullptr;
	}

	if (next_mb != nullptr)
	{
		if (next_mb->belong(xk, yk) == true)
		{
			next = next_mb;
		}
	}

	if (next == nullptr)
	{
		for (auto& i : now->Grans)
		{
			if (i->Sosed != nullptr)
			{
				if (i->Sosed->belong(xk, yk) == true)
				{
					next = i->Sosed;
				}
			}
		}
		if (next == nullptr)
		{
			for (auto& k : now->Grans)
			{
				if (k->Sosed != nullptr)
				{
					for (auto& i : k->Sosed->Grans)
					{
						if (i->Sosed != nullptr)
						{
							if (i->Sosed->belong(xk, yk) == true)
							{
								next = i->Sosed;
							}
						}
					}
				}
			}
		}
	}

	if (next == nullptr)
	{
		for (auto& i : this->All_Cells_zero)
		{
			if (i->belong(xk, yk) == true)
			{
				next = i;
			}
		}
	}

	if (next == nullptr)
	{
		this->Smut.lock();
		cout << "Setka.cpp    " << "Fly_ex_Korol   " << "Long   find  " << endl;
		cout << "Setka.cpp    " << "Fly_ex_Korol  end " << xk << " " << yk << endl;
		cout << "Setka.cpp    " << "Fly_ex_Korol  start " << x_0 << " " << y_start << endl;
		this->Smut.unlock();
		//exit(-10);
		//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< x_0 << " " << sqrt(kvv(0.0, y_0, z_0)) << endl;
		//t_per = 100000000.0;
		//if (true)
		//{
		//	double A, B, C;
		//	for (auto& i : now->Grans)
		//	{
		//		double xx, yy;
		//		i->Get_Center(xx, yy);
		//		cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "centr = " << xx << " " << yy << endl;
		//		cout << "Setka.cpp    " << "Fly_ex_Korol   "<< t_per << endl;
		//		if (i->type == Axis)
		//		{
		//			continue;
		//		}
		//		A = i->aa;
		//		B = i->bb;
		//		C = i->cc;
		//		if (fabs(B) < 0.00001)
		//		{
		//			tt1 = (-C / A - x_0) / Vx;
		//			cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "TT1 = " << tt1 << endl;
		//			if (std::fpclassify(tt1) == FP_NORMAL && tt1 > 0.0)
		//			{
		//				if (t_per > tt1)
		//				{
		//					if (now->belong(x_0 + 1.02 * tt1 * Vx, sqrt(kv(y_0 + 1.02 * tt1 * Vy) + kv(z_0 + 1.02 * tt1 * Vz))) == false)
		//					{
		//						t_per = tt1;
		//						next_mb = i->Sosed;
		//					}
		//				}
		//			}
		//			continue;
		//		}
		//		cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "1 - " << A << " " << B << " " << C << endl;
		//		//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "2 - " << x_0 << " " << y_0 << " " << z_0 << endl;
		//		//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "3 - " << Vx << " " << Vy << " " << Vz << endl;
		//		a1 = 4.0 * kv(A * Vx * (C + A * x_0) - kv(B) * (Vy * y_0 + Vz * z_0)) - 4.0 * (kv(A * Vx) - //
		//			B * B * (kv(Vy) + kv(Vz))) * (kv(C + A * x_0) - B * B * (kv(y_0) + kv(z_0)));
		//		a2 = -A * Vx * (C + A * x_0) + B * B * (Vy * y_0 + Vz * z_0);
		//		a3 = A * A * Vx * Vx - B * B * (Vy * Vy + Vz * Vz);
		//		cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "4 - " << a1 << " " << a2 << " " << a3 << endl;
		//		if (fabs(a1) < 0.000001)
		//		{
		//			a1 = 0.0;
		//		}
		//		if (a1 >= 0.0)
		//		{
		//			tt1 = (a2 + 0.5 * sqrt(a1)) / a3;
		//			tt2 = (a2 - 0.5 * sqrt(a1)) / a3;
		//			//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< A << " " << B << " " << C << endl;
		//			cout << "Setka.cpp    " << "Fly_ex_Korol   "<< tt1 << " " << tt2 << endl;
		//			if (std::fpclassify(tt1) == FP_NORMAL && tt1 > 0.0)
		//			{
		//				if (t_per > tt1)
		//				{
		//					if (now->belong(x_0 + 1.03 * tt1 * Vx, sqrt(kv(y_0 + 1.03 * tt1 * Vy) + kv(z_0 + 1.03 * tt1 * Vz))) == false)
		//					{
		//						cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "A " << endl;
		//						t_per = tt1;
		//					}
		//				}
		//			}
		//			if (std::fpclassify(tt2) == FP_NORMAL && tt2 > 0.0)
		//			{
		//				if (t_per > tt2)
		//				{
		//					if (now->belong(x_0 + 1.03 * tt2 * Vx, sqrt(kv(y_0 + 1.03 * tt2 * Vy) + kv(z_0 + 1.03 * tt2 * Vz))) == false)
		//					{
		//						cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "B " << endl;
		//						t_per = tt2;
		//					}
		//				}
		//			}
		//		}
		//	}
		//	cout << "Setka.cpp    " << "Fly_ex_Korol   "<< t_per << endl;
		//	xk = x_0 + t_per * Vx;
		//	yk = sqrt(kv(y_0 + t_per * Vy) + kv(z_0 + t_per * Vz));
		//	cout << "Setka.cpp    " << "Fly_ex_Korol   "<< xk << " " << yk << endl;

		//	a = t_per * 0.99;
		//	b = t_per * 1.01;
		//}




		//exit(-1);
		for (auto& i : this->All_Cells)
		{
			if (i->belong(xk, yk) == true)
			{
				next = i;
				break;
			}
		}
	}

	if (next != nullptr)
	{
		// Для подсчёта статистики весов
		if (mu_statistic)
		{
			if (now->zona != zon_stat)
			{
				Mu_stat[area][now->zona] += mu;
				I_stat[area][now->zona] ++;
				zon_stat = now->zona;
			}
		}

		Fly_exchenge_Imit_Korol(MK, sens, X, Y, Z, Vx, Vy, Vz, next, mu3, area, false, mu_start, to_I, to_J, georaschep, zon_stat);
	}
	else
	{

		if (xk > -1000.0 / RR_ && xk < 1000.0 / RR_ && yk < 1600.0 / RR_)
		{
			cout << "Setka.cpp    " << "Fly_ex_Korol   " << "ERROR  8670  poteryal atom" << endl;
			cout << "Setka.cpp    " << "Fly_ex_Korol   " << "SDFSDF    " << xk << " " << yk << endl;
		}
	}

	return;
}

void Setka::Fly_exchenge_Imit_Korol(MKmethod& MK, int& s1, int& s2, int& s3, double x_0, double y_0, double z_0, double Vx, double Vy, double Vz, Cell* now, double mu, //
	int area, bool ExCh, const double& mu_start, int to_I, int to_J, bool georaschep, int zon_stat)
	// Новая версия полёта нейтрального атома. Функция без лишнего кода и насыщена комментариями
	// Функция использует и геометрическое и физическое расщепление!!!
	// Описание переменных:
	// sens - датчик случайных чисел
	// x_0 ... -  координаты начала полёта нейтрального атома в данной ячейке
	// Vx ... -   скорость летящего нейтрального атома
	// now -      текущая ячейка
	// mu -       вес нейтрального атома
	// area -     какой атом летит
	// mu_start - начальный вес атома при инициализации (чтобы понимать, сколько от исходного осталось)
	// to_I - в какую область предназначался атом
{
	Cell* next_mb = nullptr;
	// Алгоритм вырубания лишних траектори (траектории вырубаются при пересечении новой зоны, далее 
	// вес лмбо увеличивается, либо частица уничтожается)
	// Не уверен, что to_I бывает больше, чем now->zona, кажется не должно такого быть

	if (false)//(area > 1 && to_I < now->zona)// (true && to_I < now->zona)
	{
		if (x_0 * Vx + y_0 * Vy + z_0 * Vz > 0.0 && to_I < now->zona) // Если летит вне солнца
		{
			double all = polar_angle(x_0, sqrt(kv(y_0) + kv(z_0)));
			to_I = now->zona;
			if (fabs(mu) < Mu[area][now->zona] * mu_start * 0.3 * sin(all))
			{
				if (fabs(mu) >= MyRandom(s1, s2, s3) * Mu[area][now->zona] * mu_start * 0.3 * sin(all))
				{
					mu = Mu[area][now->zona] * mu_start * sign(mu) * 0.3 * sin(all);
				}
				else
				{
					return;
				}
			}
		}
	}

	/*if (now->number == 34)
	{
		if (fabs(mu) > Mu[area][now->zona] * mu_start * 0.03)
		{
			this->mut_1.lock();
			this->k_1++;
			this->mut_1.unlock();
		}
	}*/

	double normV = sqrt(kvv(Vx, Vy, Vz));
	double dt = now->L / normV;  // Время на один мини-шаг в ячейке (планируется сделать много шагов в ячейке)
	double ddt = 0.001 * dt;
	// Нахождение времени до попадания в следующую ячейку
	//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "V = " << Vx << " " << Vy << " " << Vz << endl;
	//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "A" << endl;
	double a = 0.0, b = dt;
	double t_per = 100000000000.0;
	// Блок проверки перигелия к оси
	/*double t_per = (-Vy * y_0 - Vz * z_0) / (kv(Vy) + kv(Vz));
	if (t_per > 0.0000001)
	{
		if (now->belong(x_0 + t_per * Vx, sqrt(kv(y_0 + t_per * Vy) + kv(z_0 + t_per * Vz))) == true)
		{
			t_per = -1.0;
		}
	}
	else
	{
		t_per = -1.0;
	}*/

	//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "A" << endl;
	double tt1, tt2;
	double a1, a2, a3;

	// Пытаемся определить время выхода точно
	if (true)
	{
		double A, B, C;
		for (auto& i : now->Grans)
		{
			if (i->type == Axis)
			{
				continue;
			}
			A = i->aa;
			B = i->bb;
			C = i->cc;

			/*if (fabs(B) < 0.001)
			{
				tt1 = (-C / A - x_0) / Vx;
				if (std::fpclassify(tt1) == FP_NORMAL && tt1 > 0.0)
				{
					if (t_per > tt1)
					{
						if (now->belong(x_0 + (tt1 + min(ddt, 0.03 * tt1)) * Vx, sqrt(kv(y_0 + (tt1 + min(ddt, 0.03 * tt1)) * Vy) + kv(z_0 + (tt1 + min(ddt, 0.03 * tt1)) * Vz))) == false)
						{
							t_per = tt1;
							next_mb = i->Sosed;
						}
					}
				}
				continue;
			}*/
			/*cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "1 - " << A << " " << B << " " << C << endl;
			cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "2 - " << x_0 << " " << y_0 << " " << z_0 << endl;
			cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "3 - " << Vx << " " << Vy << " " << Vz << endl;*/
			a1 = 4.0 * kv(A * Vx * (C + A * x_0) - kv(B) * (Vy * y_0 + Vz * z_0)) - 4.0 * (kv(A * Vx) - //
				B * B * (kv(Vy) + kv(Vz))) * (kv(C + A * x_0) - B * B * (kv(y_0) + kv(z_0)));
			a2 = -A * Vx * (C + A * x_0) + B * B * (Vy * y_0 + Vz * z_0);
			a3 = A * A * Vx * Vx - B * B * (Vy * Vy + Vz * Vz);
			//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "4 - " << a1 << " " << a2 << " " << a3 << endl;
			if (fabs(a1) < 0.0000001)
			{
				a1 = 0.0;
			}
			if (a1 >= 0.0)
			{
				tt1 = (a2 + 0.5 * sqrt(a1)) / a3;
				tt2 = (a2 - 0.5 * sqrt(a1)) / a3;
				//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< A << " " << B << " " << C << endl;
				if (std::fpclassify(tt1) == FP_NORMAL && tt1 > 0.0)
				{
					if (t_per > tt1)
					{
						if (now->belong(x_0 + (tt1 + min(ddt, 0.005 * tt1)) * Vx, sqrt(kv(y_0 + (tt1 + min(ddt, 0.005 * tt1)) * Vy) + kv(z_0 + (tt1 + min(ddt, 0.01 * tt1)) * Vz))) == false)
						{
							t_per = tt1;
							next_mb = i->Sosed;
						}
					}
				}
				if (std::fpclassify(tt2) == FP_NORMAL && tt2 > 0.0)
				{
					if (t_per > tt2)
					{
						if (now->belong(x_0 + (tt2 + min(ddt, 0.005 * tt2)) * Vx, sqrt(kv(y_0 + (tt2 + min(ddt, 0.005 * tt2)) * Vy) + kv(z_0 + (tt2 + min(ddt, 0.01 * tt2)) * Vz))) == false)
						{
							t_per = tt2;
							next_mb = i->Sosed;
						}
					}
				}
			}
		}

		a = t_per - min(ddt, 0.005 * t_per);
		b = t_per + min(ddt, 0.005 * t_per);
	}

	bool hand = false;
	if (now->belong(x_0 + a * Vx, sqrt(kv(y_0 + a * Vy) + kv(z_0 + a * Vz))) == false)
	{
		hand = true;
	}
	else if (now->belong(x_0 + b * Vx, sqrt(kv(y_0 + b * Vy) + kv(z_0 + b * Vz))) == true)
	{
		hand = true;
	}
	else if (t_per > 100000000.0)
	{
		//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "ERROR  8568  " << t_per << endl;
		hand = true;
	}

	double k;

	// Находим время  в ячейке
	int lk = 0;
	if (hand)
	{
		a = 0.0;
		b = dt;
		t_per = 100000000000.0;
		for (auto& i : now->Grans)
		{
			a1 = -kv(i->a * Vx) * (kv(i->a * Vx) - kv(Vy) - kv(Vz)) * kv(Vz * y_0 - Vy * z_0);
			a2 = (kv(i->a * Vx) - kv(Vy) - kv(Vz)) * (Vz * z_0 + Vy * y_0);
			a3 = (kv(Vy) + kv(Vz)) * (-kv(i->a * Vx) + kv(Vy) + kv(Vz));
			if (a1 >= 0)
			{
				tt1 = (sqrt(a1) + a2) / a3;
				tt2 = (-sqrt(a1) + a2) / a3;
				if (std::fpclassify(tt1) == FP_NORMAL && tt1 > 0.0)
				{
					if (now->belong(x_0 + tt1 * Vx, sqrt(kv(y_0 + tt1 * Vy) + kv(z_0 + tt1 * Vz))) == false)
					{
						if (t_per > tt1)
						{
							t_per = tt1;
						}
					}
				}
				if (std::fpclassify(tt2) == FP_NORMAL && tt2 > 0.0)
				{
					if (now->belong(x_0 + tt2 * Vx, sqrt(kv(y_0 + tt2 * Vy) + kv(z_0 + tt2 * Vz))) == false)
					{
						if (t_per > tt2)
						{
							t_per = tt2;
						}
					}
				}
			}
		}

		while (now->belong(x_0 + b * Vx, sqrt(kv(y_0 + b * Vy) + kv(z_0 + b * Vz))) == true) // Слишком маленький шаг
		{
			lk++;
			if (lk > 1000)
			{
				cout << "Setka.cpp    " << "Fly_ex_Korol   " << "8055 dt = " << dt << " " << x_0 << " " << sqrt(kv(y_0) + kv(z_0)) << endl;
				cout << "Setka.cpp    " << "Fly_ex_Korol   " << normV << " " << Vx << " " << Vy << " " << Vz << endl;
				exit(-1);
			}
			b = b * 1.1;
		}

		if (t_per < b)
		{
			b = t_per;
		}

		k = (a + b) / 2.0;
		while (fabs(a - b) * normV > now->L / 300.0)
		{
			if (now->belong(x_0 + k * Vx, sqrt(kv(y_0 + k * Vy) + kv(z_0 + k * Vz))) == true)
			{
				a = k;
			}
			else
			{
				b = k;
			}
			k = (a + b) / 2.0;
		}
	}

	double time = b; // Время нахождения в ячейке


	bool slay = false;
	double t1, t2;
	//Нужно проверить пересечение со сферой
	double D = kv(2.0 * Vx * x_0 + 2.0 * Vy * y_0 + 2.0 * Vz * z_0) - 4.0 * kvv(Vx, Vy, Vz) * (kvv(x_0, y_0, z_0) - kv(Rmax_));
	if (D >= 0.0)
	{
		t1 = (-2.0 * (Vx * x_0 + Vy * y_0 + Vz * z_0) - sqrt(D)) / (2.0 * kv(normV));
		t2 = (-2.0 * (Vx * x_0 + Vy * y_0 + Vz * z_0) + sqrt(D)) / (2.0 * kv(normV));
		if (t1 > 0.0 && t1 < time)
		{
			if (x_0 + Vx * t1 >= 0.0)
			{
				time = t1;
				a = time * 0.999;
				b = time * 1.001;
				slay = true;
			}
		}
		if (t2 > 0.0 && t2 < time)
		{
			if (x_0 + Vx * t2 >= 0.0)  // Для правой полусферы
			{
				time = t2;
				a = time * 0.999;
				b = time * 1.001;
				slay = true;
			}
		}
	}
	// Теперь пересечение с правой границей
	if (fabs(Vx) > 0.000001)
	{
		t1 = -x_0 / Vx;
		if (t1 > 0.0 && t1 <= time && kvv(0.0, y_0 + t1 * Vy, z_0 + t1 * Vz) >= kv(Rmax_))
		{
			time = t1;
			a = time * 0.999;
			b = time * 1.001;
			slay = true;
		}
	}

	//Мини - проверка
	//for (double kk = 0.0; kk < a ; kk += time / 50.0)
	//{
	//	if (now->belong(x_0 + kk * Vx, sqrt(kv(y_0 + kk * Vy) + kv(z_0 + kk * Vz))) == false)
	//	{
	//		cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "ERROR  8280 iruhfuerfefr  " << x_0 + kk * Vx << " " << sqrt(kv(y_0 + kk * Vy) + kv(z_0 + kk * Vz)) << " " << now->contour[0]->x << " " << now->contour[0]->y << " " << kk << " " << time << endl;
	//		cout << "Setka.cpp    " << "Fly_ex_Korol   "<< a << " " << b << endl;
	//		cout << "Setka.cpp    " << "Fly_ex_Korol   "<< x_0 << " " << sqrt(kv(y_0) + kv(z_0)) << endl;
	//		cout << "Setka.cpp    " << "Fly_ex_Korol   "<< x_0 + 0.01 * Vx << " " << sqrt(kv(y_0 + 0.01 * Vy) + kv(z_0 + 0.01 * Vz)) << endl;
	//		double xx, yy;
	//		now->Get_Center(xx, yy);
	//		cout << "Setka.cpp    " << "Fly_ex_Korol   "<< xx << " " << yy << endl;
	//		cout << "Setka.cpp    " << "Fly_ex_Korol   "<< " FFF " << endl;
	//		for (auto& i : now->Grans)
	//		{
	//			if (i->type == Axis)
	//			{
	//				continue;
	//			}
	//			double A = i->aa;
	//			double B = i->bb;
	//			double C = i->cc;
	//			cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "1 - " << A << " " << B << " " << C << endl;
	//			a1 = 4.0 * kv(A * Vx * (C + A * x_0) - kv(B) * (Vy * y_0 + Vz * z_0)) - 4.0 * (kv(A * Vx) - //
	//				B * B * (kv(Vy) + kv(Vz))) * (kv(C + A * x_0) - B * B * (kv(y_0) + kv(z_0)));
	//			a2 = -A * Vx * (C + A * x_0) + B * B * (Vy * y_0 + Vz * z_0);
	//			a3 = A * A * Vx * Vx - B * B * (Vy * Vy + Vz * Vz);
	//			cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "4 - " << a1 << " " << a2 << " " << a3 << endl;
	//			if (fabs(a1) < 0.000000001)
	//			{
	//				a1 = 0.0;
	//			}
	//			if (a1 >= 0.0)
	//			{
	//				tt1 = (a2 + 0.5 * sqrt(a1)) / a3;
	//				tt2 = (a2 - 0.5 * sqrt(a1)) / a3;
	//				cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "tt1 = " << tt1 << endl;
	//				cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "tt2 = " << tt2 << endl;
	//			}
	//		}
	//		exit(-1);
	//	}
	//}

	// теперь время time, a, b уже определены
	// a - последний момент времени в ячейке, b - первый момент за пределами ячейки


	//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< a << " " << b << endl;
	double uz, uz_M, uz_E;								// средние скорости в интегралах
	double cp = sqrt(now->par[0].p / now->par[0].ro);
	double vx = now->par[0].u;							// Скорости плазмы в ячейке
	double vy = now->par[0].v;
	double ro = now->par[0].ro;


	double t_ex = 0.0;									// время до перезарядки
	// time - время нахождения атома в ячейке
	t2 = time;									// время после перезарядки (будет ниже)

	double x_ex = x_0, y_ex = y_0, z_ex = z_0;		    // координаты перезарядки
	double u1, u2, u3;

	double y_start = sqrt(kv(y_0) + kv(z_0));

	double l = sqrt(kvv(time * Vx, time * Vy, time * Vz));                    // Расстояние до выхода из ячейки

	double alpha = 0.0;								    // угол 
	double u = 0.0;                                     // модуль относительной скорости атома и плазмы
	double nu_ex = 0.0;								    // частота перезарядки
	double mu_ex = 0.0;									// вес перезаряженного атома
	double mu2 = mu;									// вес не-перезаряженного атома

	double drob = 5.0;     // Для того чтобы точнее учесть проворот скорости во время полёта

	if (true)//(ExCh == false) // Если перезарядки в этой ячейке ещё не было
	{
		double kappa = 0.0;
		//drob = min(fabs(polar_angle(y_0, z_0) - polar_angle(y_0 + time * Vy, z_0 + time * Vz)) / 0.017 + 5.0, 80.0);

		//if (y_start < 40.0/RR_)
		//{
		//	drob = 30.0;  // 30.0
		//}
		//else if (y_start < 60.0 / RR_)
		//{
		//	drob = 20.0;  // 20.0
		//}
		//else if (y_start < 80.0 / RR_)
		//{
		//	drob = 15.0;
		//}
		//else if (y_start < 120.0 / RR_)
		//{
		//	drob = 10.0;
		//}


		if (true)
		{

			for (double tt = 0.0; tt < time * 0.98; tt += time / drob)
			{
				double yy = sqrt(kv(y_0 + (tt + time / (2.0 * drob)) * Vy) + kv(z_0 + (tt + time / (2.0 * drob)) * Vz));

				alpha = polar_angle(y_0 + (tt + time / (2.0 * drob)) * Vy, z_0 + (tt + time / (2.0 * drob)) * Vz);
				u = sqrt(kvv(Vx - vx, Vy - vy * cos(alpha), Vz - vy * sin(alpha)));


				if (u / cp > 7.0)
				{
					uz = Velosity_1(u, cp);
					nu_ex = (ro * uz * sigma(uz)) / Kn_;
					//cout << x_0 << " " << y_start << " " << u / cp << " " << u << " " << cp <<  endl;
				}
				else
				{
					nu_ex = (ro * MK.int_1(u, cp)) / Kn_;        // Пробуем вычислять интеграллы численно

					//uz = Velosity_1(u, cp);
					//cout << (uz * sigma(uz)) << " " << MK.int_1(u, cp) << " " << u << " " << cp << endl;
				}
				kappa += (nu_ex * time / drob);
			}
		}
		else
		{
			alpha = polar_angle(y_0, z_0);
			u = sqrt(kvv(Vx - vx, Vy - vy * cos(alpha), Vz - vy * sin(alpha)));
			if (u / cp > 7.0)
			{
				uz = Velosity_1(u, cp);
				nu_ex = (ro * uz * sigma(uz)) / Kn_;
			}
			else
			{
				nu_ex = (ro * MK.int_1(u, cp)) / Kn_;        // Пробуем вычислять интеграллы численно
			}
			kappa = (nu_ex * time);
		}


		t_ex = -(time / kappa) * log(1.0 - MyRandom(s1, s2, s3) * (1.0 - exp(-kappa))); // Время до перезарядки
		t2 = time - t_ex;   // Время сколько лететь после того, как атом перезарядился
		mu_ex = mu * (1.0 - exp(-kappa)); // вес перезаряженного атома
		mu2 = mu - mu_ex;                 // вес оставшегося неперезаряженного атома

		x_ex = x_0 + t_ex * Vx;  // Координаты перезарядки
		y_ex = y_0 + t_ex * Vy;
		z_ex = z_0 + t_ex * Vz;

		alpha = polar_angle(y_ex, z_ex);

		// проверка, если перезарядка произошла за пределами ячейки (нужно немного сдвигать её в этом случае)
		if (now->belong(x_ex, sqrt(kv(y_ex) + kv(z_ex))) == false)
		{
			x_ex = x_0 + a * Vx;  // Координаты перезарядки
			y_ex = y_0 + a * Vy;
			z_ex = z_0 + a * Vz;
			t_ex = a;
			t2 = time - t_ex;
			alpha = polar_angle(y_ex, z_ex);
		}

		double r = sqrt(kvv(x_ex, y_ex, z_ex));
		double all = polar_angle(x_ex, sqrt(kv(y_ex) + kv(z_ex)));

		// для сбора статистики
		if (func_stat)
		{
			if (r < 100.0 / RR_ && r > 60.0 / RR_)
			{
				bool sec = false;
				bool sec2 = false;
				double tt;
				double a1 = Vx * x_0 + Vy * y_0 + Vz * z_0;
				double a2 = kvv(Vx, Vy, Vz);
				double a3 = kvv(x_0, y_0, z_0) - kv(R_stat / RR_);
				if (4.0 * a1 * a1 - 4.0 * a2 * a3 >= 0.0)
				{
					double t1 = (-2.0 * a1 + sqrt(4.0 * a1 * a1 - 4.0 * a2 * a3)) / (2.0 * a2);
					double t2 = (-2.0 * a1 - sqrt(4.0 * a1 * a1 - 4.0 * a2 * a3)) / (2.0 * a2);
					if (t1 > 0.000000001 && t1 < time)
					{
						sec = true;
						tt = t1;
						if (t1 < t_ex)
						{
							sec2 = true; // пролетел до перезярядки
						}
						else
						{
							sec2 = false;
						}
					}
					if (t2 > 0.000000001 && t2 < time)
					{
						sec = true;
						tt = t2;
						if (t2 < t_ex)
						{
							sec2 = true; // пролетел до перезярядки
						}
						else
						{
							sec2 = false;
						}
					}
				}

				if (sec == true)
				{
					double xx, yy, zz;
					xx = x_0 + tt * Vx;
					yy = y_0 + tt * Vy;
					zz = z_0 + tt * Vz;
					double Vr, Vphi, Vthe;
					double phi_ = polar_angle(yy, zz);  //polar_angle(xx, sqrt(kv(yy) + kv(zz)));
					spherical_skorost(yy, zz, xx, Vy, Vz, Vx, Vr, Vphi, Vthe);

					double Vn = fabs((Vx * xx + Vy * yy + Vz * zz) / (R_stat / RR_));
					//double Vn = sqrt(kvv(Vx, Vy, Vz));
					double Vrr = (Vy * cos(phi_) + Vz * sin(phi_));
					double mumu = 0.0;
					if (sec2)
					{
						mumu = mu / Vn;
					}
					else
					{
						mumu = mu2 / Vn;
					}

					mut_mom.lock();
					this->mu_mom[area][(int)(phi_ / (pi_ / Al_stat))] += mumu;
					this->Vx_mom[area][(int)(phi_ / (pi_ / Al_stat))] += mumu * Vx;
					this->Vy_mom[area][(int)(phi_ / (pi_ / Al_stat))] += mumu * Vrr;
					this->Vxx_mom[area][(int)(phi_ / (pi_ / Al_stat))] += mumu * Vx * Vx;
					this->Vyy_mom[area][(int)(phi_ / (pi_ / Al_stat))] += mumu * Vrr * Vrr;
					this->Vxy_mom[area][(int)(phi_ / (pi_ / Al_stat))] += mumu * Vx * Vrr;
					this->Vxxx_mom[area][(int)(phi_ / (pi_ / Al_stat))] += mumu * Vx * Vx * Vx;
					this->T_mom[area][(int)(phi_ / (pi_ / Al_stat))] += mumu * kvv(Vx, Vy, Vz);
					mut_mom.unlock();
				}
			}
		}

		double rr = sqrt(kvv(0.0, y_ex, z_ex));


		// Считаем источники для частицы до перезарядки

		u = sqrt(kvv(Vx - vx, Vy - vy * cos(alpha), Vz - vy * sin(alpha)));
		u1 = vx - Vx;
		u2 = vy * cos(alpha) - Vy;
		u3 = vy * sin(alpha) - Vz;
		double skalar = Vx * u1 + Vy * u2 + Vz * u3;

		// определим область перезарядки (по геометрическим зонам)
		int i_z = 0, i_alp = 0;
		i_z = geo_zones(r);
		i_alp = alpha_zones(x_ex, rr);

		int area2 = 0;
		// определим область в которой находится атом сейчас (это параметр самой ячейки)
		if (now->type == C_5)
		{
			area2 = 3;
		}
		else if (now->type == C_4)
		{
			area2 = 2;
		}
		else if (now->type == C_1 || now->type == C_centr)
		{
			area2 = 0;
		}
		else
		{
			area2 = 1;
		}


		now->mut.lock();

		// Блок расчёта ПОЛНОЙ функции распределения на сфере, если это нужно
		if (area == 3)
		{
			if (now->df_s4_bool == true)
			{
				now->df_s4->Add_point(Vx, Vy, Vz, y_ex, z_ex,  t_ex * mu);
			}
		}
		else if (area == 2)
		{
			if (now->df_s3_bool == true)
			{
				now->df_s3->Add_point(Vx, Vy, Vz, y_ex, z_ex, t_ex * mu);
			}
		}

		//now->par[0].w_m[area] += mu / max(sin(alpha), 0.3 * Sinus[i_alp]);
		now->par[0].F_n += t_ex * mu;
		now->par[0].F_u += t_ex * Vx * mu;
		now->par[0].F_T += t_ex * kvv(Vx, Vy, Vz) * mu;
		now->par[0].H_n[area] += t_ex * mu;
		now->par[0].H_u[area] += t_ex * Vx * mu;
		now->par[0].H_T[area] += t_ex * kvv(Vx, Vy, Vz) * mu;

		double Vrr = (Vy * cos(alpha) + Vz * sin(alpha));
		now->par[0].H_uu[area] += t_ex * Vx * Vx * mu;
		now->par[0].H_uv[area] += t_ex * Vrr * Vx * mu;
		now->par[0].H_vv[area] += t_ex * Vrr * Vrr * mu;
		now->par[0].H_uuu[area] += t_ex * Vx * Vx * Vx * mu;

		if (u / cp > 7.0)
		{
			uz = Velosity_1(u, cp);
			uz_M = Velosity_2(u, cp) / (uz * kv(cp) * cp * pi * sqrtpi_);
			uz_E = Velosity_3(u, cp);

			now->par[0].I_u += -mu_ex * uz_M * u1 / u;
			now->par[0].I_v += -mu_ex * uz_M * (u2 * cos(alpha) + u3 * sin(alpha)) / u;
			now->par[0].I_T += mu_ex * (-0.25 * (3.0 * kv(cp) + 2.0 * kv(u)) * (uz_E / uz) - //
				uz_M * skalar / u);
		}
		else
		{
			double k1 = MK.int_1(u, cp);
			double k2 = MK.int_2(u, cp);
			double k3 = MK.int_3(u, cp);
			now->par[0].I_u += mu_ex * (k2 / k1) * u1 / u;
			now->par[0].I_v += mu_ex * (k2 / k1) * (u2 * cos(alpha) + u3 * sin(alpha)) / u;
			now->par[0].I_T += mu_ex * (-0.5 * k3 / k1 + k2 / k1 * skalar / u);
		}

		//drob = min(fabs(polar_angle(y_0, z_0) - polar_angle(y_0 + t_ex * Vy, z_0 + t_ex * Vz)) / 0.017 + 5.0, 80.0);

		if (true)
		{
			for (double tt = 0.0; tt < t_ex * 0.99; tt += t_ex / drob)
			{
				alpha = polar_angle(y_0 + (tt + t_ex / (2.0 * drob)) * Vy, z_0 + (tt + t_ex / (2.0 * drob)) * Vz);

				now->par[0].F_v += (t_ex / drob) * (Vy * cos(alpha) + Vz * sin(alpha)) * mu;
				now->par[0].H_v[area] += (t_ex / drob) * (Vy * cos(alpha) + Vz * sin(alpha)) * mu;
			}
		}
		else
		{
			now->par[0].F_v += t_ex * (Vy * cos(alpha) + Vz * sin(alpha)) * mu;
			now->par[0].H_v[area] += t_ex * (Vy * cos(alpha) + Vz * sin(alpha)) * mu;
		}

		//now->par[0].I_u += -mu_ex * uz_M * (sigma2(uz_M, cp) / sigma2(uz, cp)) * u1 / u;
		//now->par[0].I_v += -mu_ex * uz_M * (sigma2(uz_M, cp) / sigma2(uz, cp)) * (u2 * cos(alpha) + u3 * sin(alpha)) / u;
		//now->par[0].I_T += mu_ex * (-0.25 * (3.0 * kv(cp) + 2.0 * kv(u)) * (uz_E / uz) * (sigma2(uz_E, cp) / sigma2(uz, cp)) - //
		//	uz_M * (sigma2(uz_M, cp) / sigma2(uz, cp)) * skalar / u);

		/*now->par[0].II_u += -mu_ex * u1;
		now->par[0].II_v += -mu_ex * (u2 * cos(alpha) + u3 * sin(alpha));
		now->par[0].II_T += mu_ex * 0.5 * (kvv(Vx, Vy, Vz) - kvv(vx, vy, 0.0));*/

		now->mut.unlock();


		double Ur, Uthe, Uphi, Vr, Vthe, Vphi;
		spherical_skorost(y_ex, z_ex, x_ex, vy * cos(alpha), vy * sin(alpha), vx, Ur, Uphi, Uthe);
		spherical_skorost(y_ex, z_ex, x_ex, Vy, Vz, Vx, Vr, Vphi, Vthe);
		double mu3;


		int stat_zone_ = now->zona;

		if (area != area2)
		{
			stat_zone_ = -1;
		}


		if (area2 == 0 || Ur / cp > 3.0)   // Без геометрического расщепления
		{
		aa:
			bool kj = true;
			mu3 = mu_ex;
			double alpha = polar_angle(y_ex, z_ex);
			double aa, bb, cc;  // Декартовы скорости после перезарядки
			vector <double> Wr(1);
			vector <double> Wthe(1);
			vector <double> Wphi(1);
			vector <double> mu_(1);
			MK.Change_Velosity4(s1, s2, s3, Ur / cp, Uthe / cp, Uphi / cp, Vr / cp, Vthe / cp, Vphi / cp, Wr, Wthe, Wphi, mu_, cp, r, 0);
			Wr[0] = Wr[0] * cp;
			Wthe[0] = Wthe[0] * cp;
			Wphi[0] = Wphi[0] * cp;
			dekard_skorost(y_ex, z_ex, x_ex, Wr[0], Wphi[0], Wthe[0], bb, cc, aa);
			double peregel, time_do_peregel;
			int ii_z;
			int ii_alp;

			if (Wr[0] < 0.0)
			{
				//wwt = sqrt(kv(Wphi[0]) + kv(Wthe[0]));
				//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "Letit vniz  " << x_ex * RR_ << " " << rr * RR_ << " " << mu3 << " " << area << " " << area2 << endl;
				//goto aa;
				time_do_peregel = (-aa * x_ex - bb * y_ex - cc * z_ex) / kvv(aa, bb, cc);
				peregel = sqrt(kvv(x_ex + aa * time_do_peregel, y_ex + bb * time_do_peregel, z_ex + cc * time_do_peregel));
				ii_z = this->geo_zones(peregel);                     // Номер зоны перегелия атома
				ii_alp = alpha_zones(x_ex + aa * time_do_peregel, sqrt(kvv(y_ex + bb * time_do_peregel, z_ex + cc * time_do_peregel, 0.0)));
				//kj = false;   // Убиваем траекторию, иначе может испортить статистику
			}
			else
			{
				ii_z = i_z;
				ii_alp = i_alp;
			}

			if (kj == true)
			{
				if (fabs(mu3) >= pred(area2, ii_z, all))
				{
					kj = true;
				}
				else
				{
					if (fabs(mu3) >= MyRandom(s1, s2, s3) * pred(area2, ii_z, all))
					{
						mu3 = pred(area2, ii_z, all) * sign(mu3);
						kj = true;
					}
					else
					{
						kj = false;
					}
				}
			}

			if (kj == true)
			{
				now->mut.lock();
				now->par[0].II_u += mu3 * (Vx - aa);
				now->par[0].II_v += mu3 * ((Vy - bb) * cos(alpha) + (Vz - cc) * sin(alpha));
				now->par[0].II_T += mu3 * 0.5 * (kvv(Vx, Vy, Vz) - kvv(aa, bb, cc));
				now->mut.unlock();
				Fly_exchenge_Imit_Korol(MK, s1, s2, s3, x_ex, y_ex, z_ex, aa, bb, cc, now, mu3, area2, false, mu_start, ii_z, ii_alp, true, stat_zone_);
			}
			//this->Fly_exchenge_Imit(MK, sens, x_ex, y_ex, z_ex, aa, bb, cc, now, mu3, -log(1.0 - sens->MakeRandom()),//
			//	0.0, area2, mu_start, to_I, 0);

		}
		else  // делаем геометрическое расщепление
		{
			int I = geo_zones(r, 1.2);  // Число дополнительных траекторий
			if (I > to_I)  // Для того чтобы не расщеплять атом, который летел в to_I область на атомы из более крупных областей
			{
				I = to_I;
			}

			vector <double> Wr(I + 1);
			vector <double> Wthe(I + 1);
			vector <double> Wphi(I + 1);
			vector <double> mu_(I + 1);

			//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< " C " << endl;
			bool bbb = MK.Change_Velosity4(s1, s2, s3, Ur / cp, Uthe / cp, Uphi / cp, Vr / cp, Vthe / cp, Vphi / cp, Wr, Wthe, Wphi, mu_, cp, r, I);
			//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< " D " << endl;

			for (int i = 0; i <= I; i++)
			{
				Wr[i] = Wr[i] * cp;
				Wthe[i] = Wthe[i] * cp;
				Wphi[i] = Wphi[i] * cp;
			}

			double aa, bb, cc, time_do_peregel, peregel;
			int ii_z, ii_alp;
			bool kj = true;
			//now->par[0].num_atoms++;
			for (int i = 0; i < I; i++)
			{
				dekard_skorost(y_ex, z_ex, x_ex, Wr[i], Wphi[i], Wthe[i], bb, cc, aa);
				mu3 = mu_[i] * mu_ex;
				time_do_peregel = (-aa * x_ex - bb * y_ex - cc * z_ex) / kvv(aa, bb, cc);
				peregel = sqrt(kvv(x_ex + aa * time_do_peregel, y_ex + bb * time_do_peregel, z_ex + cc * time_do_peregel));
				ii_z = this->geo_zones(peregel);                     // Номер зоны перегелия атома
				ii_alp = alpha_zones(x_ex + aa * time_do_peregel, sqrt(kvv(y_ex + bb * time_do_peregel, z_ex + cc * time_do_peregel, 0.0)));

				if (fabs(mu3) >= pred(area2, ii_z, all))
				{
					kj = true;
				}
				else
				{
					if (fabs(mu3) >= MyRandom(s1, s2, s3) * pred(area2, ii_z, all))
					{
						mu3 = pred(area2, ii_z, all);
						kj = true;
					}
					else
					{
						kj = false;
					}
				}
				if (kj == true)
				{
					now->mut.lock();
					now->par[0].II_u += mu3 * (Vx - aa);
					now->par[0].II_v += mu3 * ((Vy - bb) * cos(alpha) + (Vz - cc) * sin(alpha));
					now->par[0].II_T += mu3 * 0.5 * (kvv(Vx, Vy, Vz) - kvv(aa, bb, cc));
					now->mut.unlock();
					Fly_exchenge_Imit_Korol(MK, s1,s2,s3, x_ex, y_ex, z_ex, aa, bb, cc, now, mu3, area2, false, mu_start, ii_z, ii_alp, true, stat_zone_);
				}
			}

			if (bbb == true)  // Если вообще нужно запускать основной атом
			{
				mu3 = mu_[I] * mu_ex;

				// Нужно определить номер основного атома (от его перегелия)
				dekard_skorost(y_ex, z_ex, x_ex, Wr[I], Wphi[I], Wthe[I], bb, cc, aa);
				if (Wr[I] < 0.0)
				{
					time_do_peregel = (-aa * x_ex - bb * y_ex - cc * z_ex) / kvv(aa, bb, cc);
					peregel = sqrt(kvv(x_ex + aa * time_do_peregel, y_ex + bb * time_do_peregel, z_ex + cc * time_do_peregel));
					ii_z = this->geo_zones(peregel);                     // Номер зоны перегелия атома
					ii_alp = alpha_zones(x_ex + aa * time_do_peregel, sqrt(kvv(y_ex + bb * time_do_peregel, z_ex + cc * time_do_peregel, 0.0)));
				}
				else
				{
					ii_z = i_z;
					ii_alp = i_alp;
				}

				if (fabs(mu3) >= pred(area2, ii_z, all))
				{
					kj = true;
				}
				else
				{
					if (fabs(mu3) >= MyRandom(s1, s2, s3) * pred(area2, ii_z, all))
					{
						mu3 = pred(area2, ii_z, all) * sign(mu3);
						kj = true;
					}
					else
					{
						kj = false;
					}
				}

				if (kj == true)
				{
					now->mut.lock();
					now->par[0].II_u += mu3 * (Vx - aa);
					now->par[0].II_v += mu3 * ((Vy - bb) * cos(alpha) + (Vz - cc) * sin(alpha));
					now->par[0].II_T += mu3 * 0.5 * (kvv(Vx, Vy, Vz) - kvv(aa, bb, cc));
					now->mut.unlock();
					Fly_exchenge_Imit_Korol(MK, s1,s2,s3, x_ex, y_ex, z_ex, aa, bb, cc, now, mu3, area2, true, mu_start, ii_z, ii_alp, false, stat_zone_);
				}
			}

		}


		// Здесь можно добавить рулетку для оставшегося неперезаряженного атома
		if (fabs(mu2) < pred(area, to_I, all))
		{
			if (fabs(mu2) >= MyRandom(s1, s2, s3) * pred(area, to_I, all))
			{
				mu2 = pred(area, to_I, all) * sign(mu2);
			}
			else
			{
				return;
			}
		}

	}


	// этот блок кода будет работать и для частиц, которые не попали в предыдущий цикл
	// здесь t2, _ex будут иметь исходные значения, как при инициализации
	alpha = polar_angle(y_ex + 0.5 * t2 * Vy, z_ex + 0.5 * t2 * Vz);
	now->mut.lock();

	// Блок расчёта функции распределения, если это нужно
	if (area == 3)
	{
		if (now->df_s4_bool == true)
		{
			now->df_s4->Add_point(Vx, Vy, Vz, y_ex, z_ex, t2 * mu2);
		}
	}
	else if (area == 2)
	{
		if (now->df_s3_bool == true)
		{
			now->df_s3->Add_point(Vx, Vy, Vz, y_ex, z_ex, t2 * mu2);
		}
	}

	now->par[0].F_n += t2 * mu2;
	now->par[0].F_u += t2 * Vx * mu2;
	now->par[0].F_T += t2 * kvv(Vx, Vy, Vz) * mu2;
	now->par[0].H_n[area] += t2 * mu2;
	now->par[0].H_u[area] += t2 * Vx * mu2;
	now->par[0].H_T[area] += t2 * kvv(Vx, Vy, Vz) * mu2;

	double Vrr = (Vy * cos(alpha) + Vz * sin(alpha));
	now->par[0].H_uu[area] += t2 * Vx * Vx * mu2;
	now->par[0].H_uv[area] += t2 * Vrr * Vx * mu2;
	now->par[0].H_vv[area] += t2 * Vrr * Vrr * mu2;
	now->par[0].H_uuu[area] += t2 * Vx * Vx * Vx * mu2;

	//drob = min(fabs(polar_angle(y_ex, z_ex) - polar_angle(y_ex + t2 * Vy, z_ex + t2 * Vz)) / 0.017 + 5.0, 80.0);

	if (true)
	{
		for (double tt = 0.0; tt < t2 * 0.99; tt += t2 / drob)  // было 20 и норм работало
		{
			alpha = polar_angle(y_ex + (tt + t2 / (2.0 * drob)) * Vy, z_ex + (tt + t2 / (2.0 * drob)) * Vz);

			now->par[0].F_v += (t2 / drob) * (Vy * cos(alpha) + Vz * sin(alpha)) * mu2;
			now->par[0].H_v[area] += (t2 / drob) * (Vy * cos(alpha) + Vz * sin(alpha)) * mu2;
		}
	}
	else
	{
		now->par[0].F_v += t2 * (Vy * cos(alpha) + Vz * sin(alpha)) * mu2;
		now->par[0].H_v[area] += t2 * (Vy * cos(alpha) + Vz * sin(alpha)) * mu2;
	}
	now->mut.unlock();

	double mu3 = mu2;
	/*double X = x_ex + t2 * Vx;
	double Y = y_ex + t2 * Vy;
	double Z = z_ex + t2 * Vz;*/


	double X = x_0 + b * Vx;
	double Y = y_0 + b * Vy;
	double Z = z_0 + b * Vz;

	Cell* next = nullptr;             // В какую ячейку попадаем следующую
	double xk, yk;
	xk = x_0 + b * Vx;
	yk = sqrt(kv(y_0 + b * Vy) + kv(z_0 + b * Vz));
	double r_k = sqrt(kv(xk) + kv(yk));
	double r_0 = sqrt(kv(x_0) + kv(y_0) + kv(z_0));

	/*if (yk < 5.0)
	{
		cout << "Setka.cpp    " << "Fly_ex_Korol   "<< b << "  " << xk << endl;
	}*/

	if (slay == true)  //Если надо убить траекторию
	{
		return;
	}

	if (xk < Left_ + 2.0 / RR_)  // ЭТОГО НЕ БЫЛО
	{
		return;
	}

	if (yk > R5_ - 2.0 / RR_)  // ЭТОГО НЕ БЫЛО
	{
		return;
	}

	if (r_k >= Rmax_ - 3.0 / RR_ && xk >= -0.5 / RR_)
	{
		return;
		next = nullptr;
	}

	if (next_mb != nullptr)
	{
		if (next_mb->belong(xk, yk) == true)
		{
			next = next_mb;
		}
	}

	if (next == nullptr)
	{
		for (auto& i : now->Grans)
		{
			if (i->Sosed != nullptr)
			{
				if (i->Sosed->belong(xk, yk) == true)
				{
					next = i->Sosed;
				}
			}
		}
		if (next == nullptr)
		{
			for (auto& k : now->Grans)
			{
				if (k->Sosed != nullptr)
				{
					for (auto& i : k->Sosed->Grans)
					{
						if (i->Sosed != nullptr)
						{
							if (i->Sosed->belong(xk, yk) == true)
							{
								next = i->Sosed;
							}
						}
					}
				}
			}
		}
	}

	if (next == nullptr)
	{
		for (auto& i : this->All_Cells_zero)
		{
			if (i->belong(xk, yk) == true)
			{
				next = i;
			}
		}
	}

	if (next == nullptr)
	{
		this->Smut.lock();
		cout << "Setka.cpp    " << "Fly_ex_Korol   " << "Long   find  " << endl;
		cout << "Setka.cpp    " << "Fly_ex_Korol  end " << xk << " " << yk << endl;
		cout << "Setka.cpp    " << "Fly_ex_Korol  start " << x_0 << " " << y_start << endl;
		this->Smut.unlock();
		//exit(-10);
		//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< x_0 << " " << sqrt(kvv(0.0, y_0, z_0)) << endl;
		//t_per = 100000000.0;
		//if (true)
		//{
		//	double A, B, C;
		//	for (auto& i : now->Grans)
		//	{
		//		double xx, yy;
		//		i->Get_Center(xx, yy);
		//		cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "centr = " << xx << " " << yy << endl;
		//		cout << "Setka.cpp    " << "Fly_ex_Korol   "<< t_per << endl;
		//		if (i->type == Axis)
		//		{
		//			continue;
		//		}
		//		A = i->aa;
		//		B = i->bb;
		//		C = i->cc;
		//		if (fabs(B) < 0.00001)
		//		{
		//			tt1 = (-C / A - x_0) / Vx;
		//			cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "TT1 = " << tt1 << endl;
		//			if (std::fpclassify(tt1) == FP_NORMAL && tt1 > 0.0)
		//			{
		//				if (t_per > tt1)
		//				{
		//					if (now->belong(x_0 + 1.02 * tt1 * Vx, sqrt(kv(y_0 + 1.02 * tt1 * Vy) + kv(z_0 + 1.02 * tt1 * Vz))) == false)
		//					{
		//						t_per = tt1;
		//						next_mb = i->Sosed;
		//					}
		//				}
		//			}
		//			continue;
		//		}
		//		cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "1 - " << A << " " << B << " " << C << endl;
		//		//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "2 - " << x_0 << " " << y_0 << " " << z_0 << endl;
		//		//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "3 - " << Vx << " " << Vy << " " << Vz << endl;
		//		a1 = 4.0 * kv(A * Vx * (C + A * x_0) - kv(B) * (Vy * y_0 + Vz * z_0)) - 4.0 * (kv(A * Vx) - //
		//			B * B * (kv(Vy) + kv(Vz))) * (kv(C + A * x_0) - B * B * (kv(y_0) + kv(z_0)));
		//		a2 = -A * Vx * (C + A * x_0) + B * B * (Vy * y_0 + Vz * z_0);
		//		a3 = A * A * Vx * Vx - B * B * (Vy * Vy + Vz * Vz);
		//		cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "4 - " << a1 << " " << a2 << " " << a3 << endl;
		//		if (fabs(a1) < 0.000001)
		//		{
		//			a1 = 0.0;
		//		}
		//		if (a1 >= 0.0)
		//		{
		//			tt1 = (a2 + 0.5 * sqrt(a1)) / a3;
		//			tt2 = (a2 - 0.5 * sqrt(a1)) / a3;
		//			//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< A << " " << B << " " << C << endl;
		//			cout << "Setka.cpp    " << "Fly_ex_Korol   "<< tt1 << " " << tt2 << endl;
		//			if (std::fpclassify(tt1) == FP_NORMAL && tt1 > 0.0)
		//			{
		//				if (t_per > tt1)
		//				{
		//					if (now->belong(x_0 + 1.03 * tt1 * Vx, sqrt(kv(y_0 + 1.03 * tt1 * Vy) + kv(z_0 + 1.03 * tt1 * Vz))) == false)
		//					{
		//						cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "A " << endl;
		//						t_per = tt1;
		//					}
		//				}
		//			}
		//			if (std::fpclassify(tt2) == FP_NORMAL && tt2 > 0.0)
		//			{
		//				if (t_per > tt2)
		//				{
		//					if (now->belong(x_0 + 1.03 * tt2 * Vx, sqrt(kv(y_0 + 1.03 * tt2 * Vy) + kv(z_0 + 1.03 * tt2 * Vz))) == false)
		//					{
		//						cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "B " << endl;
		//						t_per = tt2;
		//					}
		//				}
		//			}
		//		}
		//	}
		//	cout << "Setka.cpp    " << "Fly_ex_Korol   "<< t_per << endl;
		//	xk = x_0 + t_per * Vx;
		//	yk = sqrt(kv(y_0 + t_per * Vy) + kv(z_0 + t_per * Vz));
		//	cout << "Setka.cpp    " << "Fly_ex_Korol   "<< xk << " " << yk << endl;

		//	a = t_per * 0.99;
		//	b = t_per * 1.01;
		//}




		//exit(-1);
		for (auto& i : this->All_Cells)
		{
			if (i->belong(xk, yk) == true)
			{
				next = i;
				break;
			}
		}
	}

	if (next != nullptr)
	{
		// Для подсчёта статистики весов
		if (mu_statistic)
		{
			if (now->zona != zon_stat)
			{
				Mu_stat[area][now->zona] += mu;
				I_stat[area][now->zona] ++;
				zon_stat = now->zona;
			}
		}

		Fly_exchenge_Imit_Korol(MK, s1, s2, s3, X, Y, Z, Vx, Vy, Vz, next, mu3, area, false, mu_start, to_I, to_J, georaschep, zon_stat);
	}
	else
	{

		if (xk > -1000.0 / RR_ && xk < 1000.0 / RR_ && yk < 1600.0 / RR_)
		{
			cout << "Setka.cpp    " << "Fly_ex_Korol   " << "ERROR  8670  poteryal atom" << endl;
			cout << "Setka.cpp    " << "Fly_ex_Korol   " << "SDFSDF    " << xk << " " << yk << endl;
		}
	}

	return;
}

void Setka::Fly_exchenge_Imit_Korol(MKmethod& MK, Sensor* sens,  bool** AZ, double x_0, double y_0, double z_0, double Vx, double Vy, double Vz, Cell* now, double mu, //
	int area, bool ExCh, const double& mu_start, int to_I, int to_J, bool georaschep, int zon_stat)
	// Новая версия полёта нейтрального атома. Функция без лишнего кода и насыщена комментариями
	// Функция использует и геометрическое и физическое расщепление!!!
	// Описание переменных:
	// sens - датчик случайных чисел
	// x_0 ... -  координаты начала полёта нейтрального атома в данной ячейке
	// Vx ... -   скорость летящего нейтрального атома
	// now -      текущая ячейка
	// mu -       вес нейтрального атома
	// area -     какой атом летит
	// mu_start - начальный вес атома при инициализации (чтобы понимать, сколько от исходного осталось)
	// to_I - в какую область предназначался атом
{
	//cout << "BBB" << endl;
	Cell* next_mb = nullptr;
	Sensor& adv = *sens;

	//cout << x_0 << " " << sqrt(kv(y_0) + kv(z_0)) << " " << mu << endl;

	// Алгоритм вырубания лишних траектори (траектории вырубаются при пересечении новой зоны, далее 
	// вес лмбо увеличивается, либо частица уничтожается)
	// Не уверен, что to_I бывает больше, чем now->zona, кажется не должно такого быть

	if (false)//(area > 1 && to_I < now->zona)// (true && to_I < now->zona)
	{
		if (x_0 * Vx + y_0 * Vy + z_0 * Vz > 0.0 && to_I < now->zona) // Если летит вне солнца
		{
			double all = polar_angle(x_0, sqrt(kv(y_0) + kv(z_0)));
			to_I = now->zona;
			if (fabs(mu) < Mu[area][now->zona] * mu_start * 0.3 * sin(all))
			{
				if (fabs(mu) >= adv.MakeRandom() * Mu[area][now->zona] * mu_start * 0.3 * sin(all))
				{
					mu = Mu[area][now->zona] * mu_start * sign(mu) * 0.3 * sin(all);
				}
				else
				{
					return;
				}
			}
		}
	}

	/*if (now->number == 34)
	{
		if (fabs(mu) > Mu[area][now->zona] * mu_start * 0.03)
		{
			this->mut_1.lock();
			this->k_1++;
			this->mut_1.unlock();
		}
	}*/

	double normV = sqrt(kvv(Vx, Vy, Vz));
	double dt = now->L / normV;  // Время на один мини-шаг в ячейке (планируется сделать много шагов в ячейке)
	double ddt = 0.001 * dt;
	// Нахождение времени до попадания в следующую ячейку
	//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "V = " << Vx << " " << Vy << " " << Vz << endl;
	//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "A" << endl;
	double a = 0.0, b = dt;
	double t_per = 100000000000.0;
	// Блок проверки перигелия к оси
	/*double t_per = (-Vy * y_0 - Vz * z_0) / (kv(Vy) + kv(Vz));
	if (t_per > 0.0000001)
	{
		if (now->belong(x_0 + t_per * Vx, sqrt(kv(y_0 + t_per * Vy) + kv(z_0 + t_per * Vz))) == true)
		{
			t_per = -1.0;
		}
	}
	else
	{
		t_per = -1.0;
	}*/

	//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "A" << endl;
	double tt1, tt2;
	double a1, a2, a3;

	// Пытаемся определить время выхода точно
	if (true)
	{
		double A, B, C;
		for (auto& i : now->Grans)
		{
			if (i->type == Axis)
			{
				continue;
			}
			A = i->aa;
			B = i->bb;
			C = i->cc;

			/*if (fabs(B) < 0.001)
			{
				tt1 = (-C / A - x_0) / Vx;
				if (std::fpclassify(tt1) == FP_NORMAL && tt1 > 0.0)
				{
					if (t_per > tt1)
					{
						if (now->belong(x_0 + (tt1 + min(ddt, 0.03 * tt1)) * Vx, sqrt(kv(y_0 + (tt1 + min(ddt, 0.03 * tt1)) * Vy) + kv(z_0 + (tt1 + min(ddt, 0.03 * tt1)) * Vz))) == false)
						{
							t_per = tt1;
							next_mb = i->Sosed;
						}
					}
				}
				continue;
			}*/
			/*cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "1 - " << A << " " << B << " " << C << endl;
			cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "2 - " << x_0 << " " << y_0 << " " << z_0 << endl;
			cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "3 - " << Vx << " " << Vy << " " << Vz << endl;*/
			a1 = 4.0 * kv(A * Vx * (C + A * x_0) - kv(B) * (Vy * y_0 + Vz * z_0)) - 4.0 * (kv(A * Vx) - //
				B * B * (kv(Vy) + kv(Vz))) * (kv(C + A * x_0) - B * B * (kv(y_0) + kv(z_0)));
			a2 = -A * Vx * (C + A * x_0) + B * B * (Vy * y_0 + Vz * z_0);
			a3 = A * A * Vx * Vx - B * B * (Vy * Vy + Vz * Vz);
			//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "4 - " << a1 << " " << a2 << " " << a3 << endl;
			if (fabs(a1) < 0.0000001)
			{
				a1 = 0.0;
			}
			if (a1 >= 0.0)
			{
				tt1 = (a2 + 0.5 * sqrt(a1)) / a3;
				tt2 = (a2 - 0.5 * sqrt(a1)) / a3;
				//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< A << " " << B << " " << C << endl;
				if (std::fpclassify(tt1) == FP_NORMAL && tt1 > 0.0)
				{
					if (t_per > tt1)
					{
						if (now->belong(x_0 + (tt1 + min(ddt, 0.005 * tt1)) * Vx, sqrt(kv(y_0 + (tt1 + min(ddt, 0.005 * tt1)) * Vy) + kv(z_0 + (tt1 + min(ddt, 0.01 * tt1)) * Vz))) == false)
						{
							t_per = tt1;
							next_mb = i->Sosed;
						}
					}
				}
				if (std::fpclassify(tt2) == FP_NORMAL && tt2 > 0.0)
				{
					if (t_per > tt2)
					{
						if (now->belong(x_0 + (tt2 + min(ddt, 0.005 * tt2)) * Vx, sqrt(kv(y_0 + (tt2 + min(ddt, 0.005 * tt2)) * Vy) + kv(z_0 + (tt2 + min(ddt, 0.01 * tt2)) * Vz))) == false)
						{
							t_per = tt2;
							next_mb = i->Sosed;
						}
					}
				}
			}
		}

		a = t_per - min(ddt, 0.005 * t_per);
		b = t_per + min(ddt, 0.005 * t_per);
	}

	bool hand = false;
	if (now->belong(x_0 + a * Vx, sqrt(kv(y_0 + a * Vy) + kv(z_0 + a * Vz))) == false)
	{
		hand = true;
	}
	else if (now->belong(x_0 + b * Vx, sqrt(kv(y_0 + b * Vy) + kv(z_0 + b * Vz))) == true)
	{
		hand = true;
	}
	else if (t_per > 100000000.0)
	{
		//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "ERROR  8568  " << t_per << endl;
		hand = true;
	}

	double k;

	// Находим время  в ячейке
	int lk = 0;
	if (hand)
	{
		a = 0.0;
		b = dt;
		t_per = 100000000000.0;
		for (auto& i : now->Grans)
		{
			a1 = -kv(i->a * Vx) * (kv(i->a * Vx) - kv(Vy) - kv(Vz)) * kv(Vz * y_0 - Vy * z_0);
			a2 = (kv(i->a * Vx) - kv(Vy) - kv(Vz)) * (Vz * z_0 + Vy * y_0);
			a3 = (kv(Vy) + kv(Vz)) * (-kv(i->a * Vx) + kv(Vy) + kv(Vz));
			if (a1 >= 0)
			{
				tt1 = (sqrt(a1) + a2) / a3;
				tt2 = (-sqrt(a1) + a2) / a3;
				if (std::fpclassify(tt1) == FP_NORMAL && tt1 > 0.0)
				{
					if (now->belong(x_0 + tt1 * Vx, sqrt(kv(y_0 + tt1 * Vy) + kv(z_0 + tt1 * Vz))) == false)
					{
						if (t_per > tt1)
						{
							t_per = tt1;
						}
					}
				}
				if (std::fpclassify(tt2) == FP_NORMAL && tt2 > 0.0)
				{
					if (now->belong(x_0 + tt2 * Vx, sqrt(kv(y_0 + tt2 * Vy) + kv(z_0 + tt2 * Vz))) == false)
					{
						if (t_per > tt2)
						{
							t_per = tt2;
						}
					}
				}
			}
		}

		while (now->belong(x_0 + b * Vx, sqrt(kv(y_0 + b * Vy) + kv(z_0 + b * Vz))) == true) // Слишком маленький шаг
		{
			lk++;
			if (lk > 1000)
			{
				cout << "Setka.cpp    " << "Fly_ex_Korol   " << "8055 dt = " << dt << " " << x_0 << " " << sqrt(kv(y_0) + kv(z_0)) << endl;
				cout << "Setka.cpp    " << "Fly_ex_Korol   " << normV << " " << Vx << " " << Vy << " " << Vz << endl;
				exit(-1);
			}
			b = b * 1.1;
		}

		if (t_per < b)
		{
			b = t_per;
		}

		k = (a + b) / 2.0;
		while (fabs(a - b) * normV > now->L / 300.0)
		{
			if (now->belong(x_0 + k * Vx, sqrt(kv(y_0 + k * Vy) + kv(z_0 + k * Vz))) == true)
			{
				a = k;
			}
			else
			{
				b = k;
			}
			k = (a + b) / 2.0;
		}
	}

	double time = b; // Время нахождения в ячейке


	bool slay = false;
	double t1, t2;
	//Нужно проверить пересечение со сферой
	double D = kv(2.0 * Vx * x_0 + 2.0 * Vy * y_0 + 2.0 * Vz * z_0) - 4.0 * kvv(Vx, Vy, Vz) * (kvv(x_0, y_0, z_0) - kv(Rmax_));
	if (D >= 0.0)
	{
		t1 = (-2.0 * (Vx * x_0 + Vy * y_0 + Vz * z_0) - sqrt(D)) / (2.0 * kv(normV));
		t2 = (-2.0 * (Vx * x_0 + Vy * y_0 + Vz * z_0) + sqrt(D)) / (2.0 * kv(normV));
		if (t1 > 0.0 && t1 < time)
		{
			if (x_0 + Vx * t1 >= 0.0)
			{
				time = t1;
				a = time * 0.999;
				b = time * 1.001;
				slay = true;
			}
		}
		if (t2 > 0.0 && t2 < time)
		{
			if (x_0 + Vx * t2 >= 0.0)  // Для правой полусферы
			{
				time = t2;
				a = time * 0.999;
				b = time * 1.001;
				slay = true;
			}
		}
	}
	// Теперь пересечение с правой границей
	if (fabs(Vx) > 0.000001)
	{
		t1 = -x_0 / Vx;
		if (t1 > 0.0 && t1 <= time && kvv(0.0, y_0 + t1 * Vy, z_0 + t1 * Vz) >= kv(Rmax_))
		{
			time = t1;
			a = time * 0.999;
			b = time * 1.001;
			slay = true;
		}
	}

	//Мини - проверка
	//for (double kk = 0.0; kk < a ; kk += time / 50.0)
	//{
	//	if (now->belong(x_0 + kk * Vx, sqrt(kv(y_0 + kk * Vy) + kv(z_0 + kk * Vz))) == false)
	//	{
	//		cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "ERROR  8280 iruhfuerfefr  " << x_0 + kk * Vx << " " << sqrt(kv(y_0 + kk * Vy) + kv(z_0 + kk * Vz)) << " " << now->contour[0]->x << " " << now->contour[0]->y << " " << kk << " " << time << endl;
	//		cout << "Setka.cpp    " << "Fly_ex_Korol   "<< a << " " << b << endl;
	//		cout << "Setka.cpp    " << "Fly_ex_Korol   "<< x_0 << " " << sqrt(kv(y_0) + kv(z_0)) << endl;
	//		cout << "Setka.cpp    " << "Fly_ex_Korol   "<< x_0 + 0.01 * Vx << " " << sqrt(kv(y_0 + 0.01 * Vy) + kv(z_0 + 0.01 * Vz)) << endl;
	//		double xx, yy;
	//		now->Get_Center(xx, yy);
	//		cout << "Setka.cpp    " << "Fly_ex_Korol   "<< xx << " " << yy << endl;
	//		cout << "Setka.cpp    " << "Fly_ex_Korol   "<< " FFF " << endl;
	//		for (auto& i : now->Grans)
	//		{
	//			if (i->type == Axis)
	//			{
	//				continue;
	//			}
	//			double A = i->aa;
	//			double B = i->bb;
	//			double C = i->cc;
	//			cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "1 - " << A << " " << B << " " << C << endl;
	//			a1 = 4.0 * kv(A * Vx * (C + A * x_0) - kv(B) * (Vy * y_0 + Vz * z_0)) - 4.0 * (kv(A * Vx) - //
	//				B * B * (kv(Vy) + kv(Vz))) * (kv(C + A * x_0) - B * B * (kv(y_0) + kv(z_0)));
	//			a2 = -A * Vx * (C + A * x_0) + B * B * (Vy * y_0 + Vz * z_0);
	//			a3 = A * A * Vx * Vx - B * B * (Vy * Vy + Vz * Vz);
	//			cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "4 - " << a1 << " " << a2 << " " << a3 << endl;
	//			if (fabs(a1) < 0.000000001)
	//			{
	//				a1 = 0.0;
	//			}
	//			if (a1 >= 0.0)
	//			{
	//				tt1 = (a2 + 0.5 * sqrt(a1)) / a3;
	//				tt2 = (a2 - 0.5 * sqrt(a1)) / a3;
	//				cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "tt1 = " << tt1 << endl;
	//				cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "tt2 = " << tt2 << endl;
	//			}
	//		}
	//		exit(-1);
	//	}
	//}

	// теперь время time, a, b уже определены
	// a - последний момент времени в ячейке, b - первый момент за пределами ячейки


	//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< a << " " << b << endl;
	double uz, uz_M, uz_E;								// средние скорости в интегралах
	double cp = sqrt(now->par[0].p / now->par[0].ro);
	double vx = now->par[0].u;							// Скорости плазмы в ячейке
	double vy = now->par[0].v;
	double ro = now->par[0].ro;


	double t_ex = 0.0;									// время до перезарядки
	// time - время нахождения атома в ячейке
	t2 = time;									// время после перезарядки (будет ниже)

	double x_ex = x_0, y_ex = y_0, z_ex = z_0;		    // координаты перезарядки
	double u1, u2, u3;

	double y_start = sqrt(kv(y_0) + kv(z_0));

	double l = sqrt(kvv(time * Vx, time * Vy, time * Vz));                    // Расстояние до выхода из ячейки

	double alpha = 0.0;								    // угол 
	double u = 0.0;                                     // модуль относительной скорости атома и плазмы
	double nu_ex = 0.0;								    // частота перезарядки
	double mu_ex = 0.0;									// вес перезаряженного атома
	double mu2 = mu;									// вес не-перезаряженного атома

	double drob = 5.0;     // Для того чтобы точнее учесть проворот скорости во время полёта

	//cout << "1" << endl;

	int i_z = 0, i_alp = 0;  // Области перезярядки по геометрическим зонам

	if (true)//(ExCh == false) // Если перезарядки в этой ячейке ещё не было
	{
		double kappa = 0.0;
		//drob = min(fabs(polar_angle(y_0, z_0) - polar_angle(y_0 + time * Vy, z_0 + time * Vz)) / 0.017 + 5.0, 80.0);

		//if (y_start < 40.0/RR_)
		//{
		//	drob = 30.0;  // 30.0
		//}
		//else if (y_start < 60.0 / RR_)
		//{
		//	drob = 20.0;  // 20.0
		//}
		//else if (y_start < 80.0 / RR_)
		//{
		//	drob = 15.0;
		//}
		//else if (y_start < 120.0 / RR_)
		//{
		//	drob = 10.0;
		//}


		if (true)
		{

			for (double tt = 0.0; tt < time * 0.98; tt += time / drob)
			{
				double yy = sqrt(kv(y_0 + (tt + time / (2.0 * drob)) * Vy) + kv(z_0 + (tt + time / (2.0 * drob)) * Vz));

				alpha = polar_angle(y_0 + (tt + time / (2.0 * drob)) * Vy, z_0 + (tt + time / (2.0 * drob)) * Vz);
				u = sqrt(kvv(Vx - vx, Vy - vy * cos(alpha), Vz - vy * sin(alpha)));


				if (u / cp > 7.0)
				{
					uz = Velosity_1(u, cp);
					nu_ex = (ro * uz * sigma(uz)) / Kn_;
					//cout << x_0 << " " << y_start << " " << u / cp << " " << u << " " << cp <<  endl;
				}
				else
				{
					nu_ex = (ro * MK.int_1(u, cp)) / Kn_;        // Пробуем вычислять интеграллы численно

					//uz = Velosity_1(u, cp);
					//cout << (uz * sigma(uz)) << " " << MK.int_1(u, cp) << " " << u << " " << cp << endl;
				}
				kappa += (nu_ex * time / drob);
			}
		}
		else
		{
			alpha = polar_angle(y_0, z_0);
			u = sqrt(kvv(Vx - vx, Vy - vy * cos(alpha), Vz - vy * sin(alpha)));
			if (u / cp > 7.0)
			{
				uz = Velosity_1(u, cp);
				nu_ex = (ro * uz * sigma(uz)) / Kn_;
			}
			else
			{
				nu_ex = (ro * MK.int_1(u, cp)) / Kn_;        // Пробуем вычислять интеграллы численно
			}
			kappa = (nu_ex * time);
		}


		t_ex = -(time / kappa) * log(1.0 - adv.MakeRandom() * (1.0 - exp(-kappa))); // Время до перезарядки
		t2 = time - t_ex;   // Время сколько лететь после того, как атом перезарядился
		mu_ex = mu * (1.0 - exp(-kappa)); // вес перезаряженного атома
		mu2 = mu - mu_ex;                 // вес оставшегося неперезаряженного атома

		x_ex = x_0 + t_ex * Vx;  // Координаты перезарядки
		y_ex = y_0 + t_ex * Vy;
		z_ex = z_0 + t_ex * Vz;

		alpha = polar_angle(y_ex, z_ex);

		// проверка, если перезарядка произошла за пределами ячейки (нужно немного сдвигать её в этом случае)
		if (now->belong(x_ex, sqrt(kv(y_ex) + kv(z_ex))) == false)
		{
			x_ex = x_0 + a * Vx;  // Координаты перезарядки
			y_ex = y_0 + a * Vy;
			z_ex = z_0 + a * Vz;
			t_ex = a;
			t2 = time - t_ex;
			alpha = polar_angle(y_ex, z_ex);
		}

		double r = sqrt(kvv(x_ex, y_ex, z_ex));
		double all = polar_angle(x_ex, sqrt(kv(y_ex) + kv(z_ex)));

		// для сбора статистики
		if (func_stat)
		{
			if (r < 100.0 / RR_ && r > 60.0 / RR_)
			{
				bool sec = false;
				bool sec2 = false;
				double tt;
				double a1 = Vx * x_0 + Vy * y_0 + Vz * z_0;
				double a2 = kvv(Vx, Vy, Vz);
				double a3 = kvv(x_0, y_0, z_0) - kv(R_stat / RR_);
				if (4.0 * a1 * a1 - 4.0 * a2 * a3 >= 0.0)
				{
					double t1 = (-2.0 * a1 + sqrt(4.0 * a1 * a1 - 4.0 * a2 * a3)) / (2.0 * a2);
					double t2 = (-2.0 * a1 - sqrt(4.0 * a1 * a1 - 4.0 * a2 * a3)) / (2.0 * a2);
					if (t1 > 0.000000001 && t1 < time)
					{
						sec = true;
						tt = t1;
						if (t1 < t_ex)
						{
							sec2 = true; // пролетел до перезярядки
						}
						else
						{
							sec2 = false;
						}
					}
					if (t2 > 0.000000001 && t2 < time)
					{
						sec = true;
						tt = t2;
						if (t2 < t_ex)
						{
							sec2 = true; // пролетел до перезярядки
						}
						else
						{
							sec2 = false;
						}
					}
				}

				if (sec == true)
				{
					double xx, yy, zz;
					xx = x_0 + tt * Vx;
					yy = y_0 + tt * Vy;
					zz = z_0 + tt * Vz;
					double Vr, Vphi, Vthe;
					double phi_ = polar_angle(yy, zz);  //polar_angle(xx, sqrt(kv(yy) + kv(zz)));
					spherical_skorost(yy, zz, xx, Vy, Vz, Vx, Vr, Vphi, Vthe);

					double Vn = fabs((Vx * xx + Vy * yy + Vz * zz) / (R_stat / RR_));
					//double Vn = sqrt(kvv(Vx, Vy, Vz));
					double Vrr = (Vy * cos(phi_) + Vz * sin(phi_));
					double mumu = 0.0;
					if (sec2)
					{
						mumu = mu / Vn;
					}
					else
					{
						mumu = mu2 / Vn;
					}

					mut_mom.lock();
					this->mu_mom[area][(int)(phi_ / (pi_ / Al_stat))] += mumu;
					this->Vx_mom[area][(int)(phi_ / (pi_ / Al_stat))] += mumu * Vx;
					this->Vy_mom[area][(int)(phi_ / (pi_ / Al_stat))] += mumu * Vrr;
					this->Vxx_mom[area][(int)(phi_ / (pi_ / Al_stat))] += mumu * Vx * Vx;
					this->Vyy_mom[area][(int)(phi_ / (pi_ / Al_stat))] += mumu * Vrr * Vrr;
					this->Vxy_mom[area][(int)(phi_ / (pi_ / Al_stat))] += mumu * Vx * Vrr;
					this->Vxxx_mom[area][(int)(phi_ / (pi_ / Al_stat))] += mumu * Vx * Vx * Vx;
					mut_mom.unlock();

					//if (fabs(mu) > pred(area, now->zona, all) * 0.3)
					//{
					//	if (phi_ < 0.17453293)  // Угол в 10 градусов
					//	{
					//		mut_stat.lock();
					//		this->V_r_stat[this->number_stat] = Vr;
					//		this->V_t_stat[this->number_stat] = Vthe;
					//		this->V_p_stat[this->number_stat] = Vphi;
					//		this->phi_stat[this->number_stat] = phi_;
					//		this->num_stat[this->number_stat] = area;
					//		if (sec2)
					//		{
					//			this->mu_stat[this->number_stat] = mu / Vn;
					//		}
					//		else
					//		{
					//			this->mu_stat[this->number_stat] = mu2 / Vn;
					//		}
					//		this->number_stat++;
					//		mut_stat.unlock();
					//	}
					//}
				}
			}
		}

		double rr = sqrt(kvv(0.0, y_ex, z_ex));

		//cout << "2" << endl;
		// Считаем источники для частицы до перезарядки

		u = sqrt(kvv(Vx - vx, Vy - vy * cos(alpha), Vz - vy * sin(alpha)));
		u1 = vx - Vx;
		u2 = vy * cos(alpha) - Vy;
		u3 = vy * sin(alpha) - Vz;
		double skalar = Vx * u1 + Vy * u2 + Vz * u3;

		// определим область перезарядки (по геометрическим зонам)
		i_z = geo_zones(r);
		i_alp = alpha_zones(x_ex, rr);

		int area2 = 0;
		// определим область в которой находится атом сейчас (это параметр самой ячейки)
		if (now->type == C_5)
		{
			area2 = 3;
		}
		else if (now->type == C_4)
		{
			area2 = 2;
		}
		else if (now->type == C_1 || now->type == C_centr)
		{
			area2 = 0;
		}
		else
		{
			area2 = 1;
		}


		now->mut.lock();

		// Блок расчёта функции распределения, если это нужно
		if (area == 3)
		{
			if (now->df_s4_bool == true)
			{
				now->df_s4->Add_point(Vx, Vy, Vz, y_ex, z_ex, t_ex * mu);
			}
		}
		else if (area == 2)
		{
			if (now->df_s3_bool == true)
			{
				now->df_s3->Add_point(Vx, Vy, Vz, y_ex, z_ex, t_ex * mu);
			}
		}
		else if (area == 1)
		{
			if (now->df_s2_bool == true)
			{
				now->df_s2->Add_point(Vx, Vy, Vz, y_ex, z_ex, t_ex * mu);
			}
		}
		else if (area == 0)
		{
			if (now->df_s1_bool == true)
			{
				now->df_s1->Add_point(Vx, Vy, Vz, y_ex, z_ex, t_ex * mu);
			}
		}

		//now->par[0].w_m[area] += mu / max(sin(alpha), 0.3 * Sinus[i_alp]);
		now->par[0].F_n += t_ex * mu;
		now->par[0].F_u += t_ex * Vx * mu;
		now->par[0].F_T += t_ex * kvv(Vx, Vy, Vz) * mu;
		now->par[0].H_n[area] += t_ex * mu;
		now->par[0].H_u[area] += t_ex * Vx * mu;
		now->par[0].H_T[area] += t_ex * kvv(Vx, Vy, Vz) * mu;

		double Vrr = (Vy * cos(alpha) + Vz * sin(alpha));
		now->par[0].H_uu[area] += t_ex * Vx * Vx * mu;
		now->par[0].H_uv[area] += t_ex * Vrr * Vx * mu;
		now->par[0].H_vv[area] += t_ex * Vrr * Vrr * mu;
		now->par[0].H_uuu[area] += t_ex * Vx * Vx * Vx * mu;

		if (u / cp > 7.0)
		{
			uz = Velosity_1(u, cp);
			uz_M = Velosity_2(u, cp) / (uz * kv(cp) * cp * pi * sqrtpi_);
			uz_E = Velosity_3(u, cp);

			now->par[0].I_u += -mu_ex * uz_M * u1 / u;
			now->par[0].I_v += -mu_ex * uz_M * (u2 * cos(alpha) + u3 * sin(alpha)) / u;
			now->par[0].I_T += mu_ex * (-0.25 * (3.0 * kv(cp) + 2.0 * kv(u)) * (uz_E / uz) - //
				uz_M * skalar / u);
		}
		else
		{
			double k1 = MK.int_1(u, cp);
			double k2 = MK.int_2(u, cp);
			double k3 = MK.int_3(u, cp);
			now->par[0].I_u += mu_ex * (k2 / k1) * u1 / u;
			now->par[0].I_v += mu_ex * (k2 / k1) * (u2 * cos(alpha) + u3 * sin(alpha)) / u;
			now->par[0].I_T += mu_ex * (-0.5 * k3 / k1 + k2 / k1 * skalar / u);
		}

		//drob = min(fabs(polar_angle(y_0, z_0) - polar_angle(y_0 + t_ex * Vy, z_0 + t_ex * Vz)) / 0.017 + 5.0, 80.0);

		if (true)
		{
			for (double tt = 0.0; tt < t_ex * 0.99; tt += t_ex / drob)
			{
				alpha = polar_angle(y_0 + (tt + t_ex / (2.0 * drob)) * Vy, z_0 + (tt + t_ex / (2.0 * drob)) * Vz);

				now->par[0].F_v += (t_ex / drob) * (Vy * cos(alpha) + Vz * sin(alpha)) * mu;
				now->par[0].H_v[area] += (t_ex / drob) * (Vy * cos(alpha) + Vz * sin(alpha)) * mu;
			}
		}
		else
		{
			now->par[0].F_v += t_ex * (Vy * cos(alpha) + Vz * sin(alpha)) * mu;
			now->par[0].H_v[area] += t_ex * (Vy * cos(alpha) + Vz * sin(alpha)) * mu;
		}

		//now->par[0].I_u += -mu_ex * uz_M * (sigma2(uz_M, cp) / sigma2(uz, cp)) * u1 / u;
		//now->par[0].I_v += -mu_ex * uz_M * (sigma2(uz_M, cp) / sigma2(uz, cp)) * (u2 * cos(alpha) + u3 * sin(alpha)) / u;
		//now->par[0].I_T += mu_ex * (-0.25 * (3.0 * kv(cp) + 2.0 * kv(u)) * (uz_E / uz) * (sigma2(uz_E, cp) / sigma2(uz, cp)) - //
		//	uz_M * (sigma2(uz_M, cp) / sigma2(uz, cp)) * skalar / u);

		/*now->par[0].II_u += -mu_ex * u1;
		now->par[0].II_v += -mu_ex * (u2 * cos(alpha) + u3 * sin(alpha));
		now->par[0].II_T += mu_ex * 0.5 * (kvv(Vx, Vy, Vz) - kvv(vx, vy, 0.0));*/

		now->mut.unlock();

		//cout << "3" << endl;
		double Ur, Uthe, Uphi, Vr, Vthe, Vphi;
		spherical_skorost(y_ex, z_ex, x_ex, vy * cos(alpha), vy * sin(alpha), vx, Ur, Uphi, Uthe);
		spherical_skorost(y_ex, z_ex, x_ex, Vy, Vz, Vx, Vr, Vphi, Vthe);
		double mu3;


		int stat_zone_ = now->zona;

		if (area != area2)
		{
			stat_zone_ = -1;
		}


		if (area2 == 0 || Ur / cp > 3.0)   // Без геометрического расщепления
		{
			//cout << "4" << endl;
		aa:
			bool kj = true;
			mu3 = mu_ex;
			double alpha = polar_angle(y_ex, z_ex);
			double aa, bb, cc;  // Декартовы скорости после перезарядки
			vector <double> Wr(1);
			vector <double> Wthe(1);
			vector <double> Wphi(1);
			vector <double> mu_(1);
			MK.Change_Velosity4(sens, Ur / cp, Uthe / cp, Uphi / cp, Vr / cp, Vthe / cp, Vphi / cp, Wr, Wthe, Wphi, mu_, cp, r, 0);
			Wr[0] = Wr[0] * cp;
			Wthe[0] = Wthe[0] * cp;
			Wphi[0] = Wphi[0] * cp;
			dekard_skorost(y_ex, z_ex, x_ex, Wr[0], Wphi[0], Wthe[0], bb, cc, aa);
			double peregel, time_do_peregel;
			int ii_z;     // В какую зону летит родившийся атом
			int ii_alp;   // В какую зону летит родившийся атом
			int to_ii;
			double r_per = r;  // радиус в перегелии

			if (Wr[0] < 0.0)
			{
				r_per = this->distination(x_ex, y_ex, z_ex, aa, bb, cc, area2, to_ii, ii_z, ii_alp);
			}
			else
			{
				to_ii = i_z;
				ii_z = i_z;
				ii_alp = i_alp;
			}

			double Muuu = max(min(Mu_[area2][ii_z][ii_alp] * 0.3 * SINKR[ii_alp] * kv(r_per / Rmax_), 
				Mu_[area2][i_z][i_alp] * 0.3 * SINKR[i_alp] * kv(r / Rmax_)), MinWes);
			//cout << "MUUU = " << Muuu << endl;

			if (kj == true)
			{
				if (fabs(mu3) >= Muuu)  // pred(area2, ii_z, all))
				{
					kj = true;
				}
				else
				{
					if (fabs(mu3) >= adv.MakeRandom() * Muuu)  // pred(area2, ii_z, all))
					{
						mu3 = Muuu * sign(mu3);  // pred(area2, ii_z, all) * sign(mu3);
						kj = true;
					}
					else
					{
						kj = false;
					}
				}
			}

			if (kj == true)
			{
				now->mut.lock();
				now->par[0].II_u += mu3 * (Vx - aa);
				now->par[0].II_v += mu3 * ((Vy - bb) * cos(alpha) + (Vz - cc) * sin(alpha));
				now->par[0].II_T += mu3 * 0.5 * (kvv(Vx, Vy, Vz) - kvv(aa, bb, cc));
				now->mut.unlock();

				if (area == area2)
				{
					bool** BZ = new bool* [I_];
					for (size_t i = 0; i < I_; i++)
					{
						BZ[i] = new bool[J_];
					}

					for (size_t i = 0; i < I_; i++)
					{
						for (size_t j = 0; j < J_; j++)
						{
							BZ[i][j] = AZ[i][j];
						}
					}
					BZ[i_z][i_alp] = true;
					//cout << "C" << endl;
					Fly_exchenge_Imit_Korol(MK, sens, BZ, x_ex, y_ex, z_ex, aa, bb, cc, now, mu3, area2, false, mu_start, to_ii, ii_alp, true, stat_zone_);
					for (int i = 0; i < I_; ++i) {
						delete[] BZ[i];
					}
					delete[] BZ;
				}
				else
				{
					bool** BZ = new bool* [I_];
					for (size_t i = 0; i < I_; i++)
					{
						BZ[i] = new bool[J_];
					}

					for (size_t i = 0; i < I_; i++)
					{
						for (size_t j = 0; j < J_; j++)
						{
							BZ[i][j] = false;
						}
					}

					//cout << "D" << endl;
					Fly_exchenge_Imit_Korol(MK, sens, BZ, x_ex, y_ex, z_ex, aa, bb, cc, now, mu3, area2, false, mu_start, ii_z, ii_alp, true, stat_zone_);
					for (int i = 0; i < I_; ++i) {
						delete[] BZ[i];
					}
					delete[] BZ;
				}
			}
			//this->Fly_exchenge_Imit(MK, sens, x_ex, y_ex, z_ex, aa, bb, cc, now, mu3, -log(1.0 - sens->MakeRandom()),//
			//	0.0, area2, mu_start, to_I, 0);

		}
		else  // делаем геометрическое расщепление
		{
		//cout << "5" << endl;
			int I = geo_zones(r, 1.2);  // Число дополнительных траекторий
			//if (I > to_I)  // Для того чтобы не расщеплять атом, который летел в to_I область на атомы из более крупных областей
			//{
			//	I = to_I;
			//}
			//cout << "5 1" << endl;
			vector <double> Wr(I + 1);
			vector <double> Wthe(I + 1);
			vector <double> Wphi(I + 1);
			vector <double> mu_(I + 1);
			//cout << "5 1 1   " << Ur << " " << cp << " " << Uphi << " " << Vr << " " << Vthe << " " << Vphi << " " << r << " " << I << endl;
			//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< " C " << endl;
			bool bbb = MK.Change_Velosity4(sens, Ur / cp, Uthe / cp, Uphi / cp, Vr / cp, Vthe / cp, Vphi / cp, Wr, Wthe, Wphi, mu_, cp, r, I);
			//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< " D " << endl;
			//cout << "5 2" << endl;
			for (int i = 0; i <= I; i++)
			{
				Wr[i] = Wr[i] * cp;
				Wthe[i] = Wthe[i] * cp;
				Wphi[i] = Wphi[i] * cp;
			}

			double aa, bb, cc, time_do_peregel, peregel;
			int ii_z, ii_alp;
			int to_ii;   // В какую зону по радиусу направлен
			bool kj = true;
			double r_per = r;
			//now->par[0].num_atoms++;
			//cout << "5 3" << endl;
			for (int i = 0; i < I; i++)
			{
				//cout << "A1" << endl;
				dekard_skorost(y_ex, z_ex, x_ex, Wr[i], Wphi[i], Wthe[i], bb, cc, aa);
				mu3 = mu_[i] * mu_ex;
				r_per = this->distination(x_ex, y_ex, z_ex, aa, bb, cc, area2, to_ii, ii_z, ii_alp);

				double Muuu = max(min(Mu_[area2][ii_z][ii_alp] * 0.3 * SINKR[ii_alp] * kv(r_per / Rmax_),
					Mu_[area2][i_z][i_alp] * 0.3 * SINKR[i_alp] * kv(r / Rmax_)), MinWes);

				if (fabs(mu3) >= Muuu)
				{
					kj = true;
				}
				else
				{
					if (fabs(mu3) >= adv.MakeRandom() * Muuu)
					{
						mu3 = Muuu;
						kj = true;
					}
					else
					{
						kj = false;
					}
				}
				//cout << "A2" << endl;
				if (kj == true)
				{
					now->mut.lock();
					now->par[0].II_u += mu3 * (Vx - aa);
					now->par[0].II_v += mu3 * ((Vy - bb) * cos(alpha) + (Vz - cc) * sin(alpha));
					now->par[0].II_T += mu3 * 0.5 * (kvv(Vx, Vy, Vz) - kvv(aa, bb, cc));
					now->mut.unlock();

					if (area == area2)
					{
						//cout << "S1" << endl;
						bool** BZ = new bool* [I_];
						for (size_t i = 0; i < I_; i++)
						{
							BZ[i] = new bool[J_];
						}

						for (size_t i = 0; i < I_; i++)
						{
							for (size_t j = 0; j < J_; j++)
							{
								BZ[i][j] = AZ[i][j];
							}
						}
						BZ[i_z][i_alp] = true;
						//cout << "S2" << endl;

						//cout << "E" << endl;
						Fly_exchenge_Imit_Korol(MK, sens, BZ, x_ex, y_ex, z_ex, aa, bb, cc, now, mu3, area2, false, mu_start, to_ii, ii_alp, true, stat_zone_);
						//cout << "E end" << endl;
						for (int i = 0; i < I_; ++i) {
							delete[] BZ[i];
						}
						delete[] BZ;
					}
					else
					{
						bool** BZ = new bool* [I_];
						for (size_t i = 0; i < I_; i++)
						{
							BZ[i] = new bool[J_];
						}

						for (size_t i = 0; i < I_; i++)
						{
							for (size_t j = 0; j < J_; j++)
							{
								BZ[i][j] = false;
							}
						}
						//cout << "F" << endl;
						Fly_exchenge_Imit_Korol(MK, sens, BZ, x_ex, y_ex, z_ex, aa, bb, cc, now, mu3, area2, false, mu_start, to_ii, ii_alp, true, stat_zone_);
						//cout << "F end" << endl;
						for (int i = 0; i < I_; ++i) {
							delete[] BZ[i];
						}
						delete[] BZ;
					}
				}
			}

			if (bbb == true)  // Если вообще нужно запускать основной атом
			{
				mu3 = mu_[I] * mu_ex;

				// Нужно определить номер основного атома (от его перегелия)
				dekard_skorost(y_ex, z_ex, x_ex, Wr[I], Wphi[I], Wthe[I], bb, cc, aa);
				double r_per = r;
				if (Wr[I] < 0.0)
				{
					r_per = this->distination(x_ex, y_ex, z_ex, aa, bb, cc, area2, to_ii, ii_z, ii_alp);
				}
				else
				{
					to_ii = i_z;
					ii_z = i_z;
					ii_alp = i_alp;
				}

				double Muuu = max(min(Mu_[area2][ii_z][ii_alp] * 0.3 * SINKR[ii_alp] * kv(r_per / Rmax_),
					Mu_[area2][i_z][i_alp] * 0.3 * SINKR[i_alp] * kv(r / Rmax_)), MinWes);
				//cout << "MUUU = " << Muuu << endl;

				if (fabs(mu3) >= Muuu)
				{
					kj = true;
				}
				else
				{
					if (fabs(mu3) >= adv.MakeRandom() * Muuu)
					{
						mu3 = Muuu * sign(mu3);
						kj = true;
					}
					else
					{
						kj = false;
					}
				}

				if (kj == true)
				{
					now->mut.lock();
					now->par[0].II_u += mu3 * (Vx - aa);
					now->par[0].II_v += mu3 * ((Vy - bb) * cos(alpha) + (Vz - cc) * sin(alpha));
					now->par[0].II_T += mu3 * 0.5 * (kvv(Vx, Vy, Vz) - kvv(aa, bb, cc));
					now->mut.unlock();
					//cout << "7" << endl;
					if (area == area2)
					{
						//cout << "do  G" << endl;
						bool** BZ = new bool* [I_];
						for (size_t i = 0; i < I_; i++)
						{
							BZ[i] = new bool[J_];
						}

						for (size_t i = 0; i < I_; i++)
						{
							for (size_t j = 0; j < J_; j++)
							{
								BZ[i][j] = AZ[i][j];
							}
						}
						BZ[i_z][i_alp] = true;

						//cout << "G" << endl;
						Fly_exchenge_Imit_Korol(MK, sens, BZ, x_ex, y_ex, z_ex, aa, bb, cc, now, mu3, area2,
							true, mu_start, ii_z, ii_alp, false, stat_zone_);
						//cout << "G end" << endl;
						for (int i = 0; i < I_; ++i) {
							delete[] BZ[i];
						}
						delete[] BZ;
					}
					else
					{
						//cout << "do H" << endl;
						bool** BZ = new bool* [I_];
						for (size_t i = 0; i < I_; i++)
						{
							BZ[i] = new bool[J_];
						}

						for (size_t i = 0; i < I_; i++)
						{
							for (size_t j = 0; j < J_; j++)
							{
								BZ[i][j] = false;
							}
						}
						
						//cout << "H" << endl;
						Fly_exchenge_Imit_Korol(MK, sens, BZ, x_ex, y_ex, z_ex, aa, bb, cc, now, mu3, area2,
							true, mu_start, ii_z, ii_alp, false, stat_zone_);
						//cout << "H end" << endl;

						for (int i = 0; i < I_; ++i) {
							delete[] BZ[i];
						}
						delete[] BZ;
					}
				}
			}

		}


		// Здесь можно добавить рулетку для оставшегося неперезаряженного атома
		double Muuu = max(min(Mu_[area2][to_I][to_J] * 0.3 * SINKR[to_J] * kv(Ri[to_I] / Rmax_),
			Mu_[area2][i_z][i_alp] * 0.3 * SINKR[i_alp] * kv(r / Rmax_)), MinWes);

		if (fabs(mu2) < Muuu)
		{
			if (fabs(mu2) >= adv.MakeRandom() * Muuu)
			{
				mu2 = Muuu * sign(mu2);
			}
			else
			{
				return;
			}
		}

	}

	


	// этот блок кода будет работать и для частиц, которые не попали в предыдущий цикл
	// здесь t2, _ex будут иметь исходные значения, как при инициализации
	alpha = polar_angle(y_ex + 0.5 * t2 * Vy, z_ex + 0.5 * t2 * Vz);
	now->mut.lock();

	// Блок расчёта функции распределения, если это нужно
	if (area == 3)
	{
		if (now->df_s4_bool == true)
		{
			now->df_s4->Add_point(Vx, Vy, Vz, y_ex, z_ex, t2 * mu2);
		}
	}
	else if (area == 2)
	{
		if (now->df_s3_bool == true)
		{
			now->df_s3->Add_point(Vx, Vy, Vz, y_ex, z_ex, t2 * mu2);
		}
	}
	else if (area == 1)
	{
		if (now->df_s2_bool == true)
		{
			now->df_s2->Add_point(Vx, Vy, Vz, y_ex, z_ex, t2 * mu2);
		}
	}
	else if (area == 0)
	{
		if (now->df_s1_bool == true)
		{
			now->df_s1->Add_point(Vx, Vy, Vz, y_ex, z_ex, t2 * mu2);
		}
	}

	now->par[0].F_n += t2 * mu2;
	now->par[0].F_u += t2 * Vx * mu2;
	now->par[0].F_T += t2 * kvv(Vx, Vy, Vz) * mu2;
	now->par[0].H_n[area] += t2 * mu2;
	now->par[0].H_u[area] += t2 * Vx * mu2;
	now->par[0].H_T[area] += t2 * kvv(Vx, Vy, Vz) * mu2;

	double Vrr = (Vy * cos(alpha) + Vz * sin(alpha));
	now->par[0].H_uu[area] += t2 * Vx * Vx * mu2;
	now->par[0].H_uv[area] += t2 * Vrr * Vx * mu2;
	now->par[0].H_vv[area] += t2 * Vrr * Vrr * mu2;
	now->par[0].H_uuu[area] += t2 * Vx * Vx * Vx * mu2;

	//drob = min(fabs(polar_angle(y_ex, z_ex) - polar_angle(y_ex + t2 * Vy, z_ex + t2 * Vz)) / 0.017 + 5.0, 80.0);

	if (true)
	{
		for (double tt = 0.0; tt < t2 * 0.99; tt += t2 / drob)  // было 20 и норм работало
		{
			alpha = polar_angle(y_ex + (tt + t2 / (2.0 * drob)) * Vy, z_ex + (tt + t2 / (2.0 * drob)) * Vz);

			now->par[0].F_v += (t2 / drob) * (Vy * cos(alpha) + Vz * sin(alpha)) * mu2;
			now->par[0].H_v[area] += (t2 / drob) * (Vy * cos(alpha) + Vz * sin(alpha)) * mu2;
		}
	}
	else
	{
		now->par[0].F_v += t2 * (Vy * cos(alpha) + Vz * sin(alpha)) * mu2;
		now->par[0].H_v[area] += t2 * (Vy * cos(alpha) + Vz * sin(alpha)) * mu2;
	}
	now->mut.unlock();

	double mu3 = mu2;
	/*double X = x_ex + t2 * Vx;
	double Y = y_ex + t2 * Vy;
	double Z = z_ex + t2 * Vz;*/
	double X = x_0 + b * Vx;
	double Y = y_0 + b * Vy;
	double Z = z_0 + b * Vz;

	Cell* next = nullptr;             // В какую ячейку попадаем следующую
	double xk, yk;
	xk = x_0 + b * Vx;
	yk = sqrt(kv(y_0 + b * Vy) + kv(z_0 + b * Vz));
	double r_k = sqrt(kv(xk) + kv(yk));
	double r_0 = sqrt(kv(x_0) + kv(y_0) + kv(z_0));

	/*if (yk < 5.0)
	{
		cout << "Setka.cpp    " << "Fly_ex_Korol   "<< b << "  " << xk << endl;
	}*/

	if (slay == true)  //Если надо убить траекторию
	{
		return;
	}

	if (xk < Left_ + 2.0 / RR_)  // ЭТОГО НЕ БЫЛО
	{
		return;
	}

	if (yk > R5_ - 2.0 / RR_)  // ЭТОГО НЕ БЫЛО
	{
		return;
	}

	if (r_k >= Rmax_ - 3.0 / RR_ && xk >= -0.5 / RR_)
	{
		return;
		next = nullptr;
	}

	if (next_mb != nullptr)
	{
		if (next_mb->belong(xk, yk) == true)
		{
			next = next_mb;
		}
	}

	if (next == nullptr)
	{
		for (auto& i : now->Grans)
		{
			if (i->Sosed != nullptr)
			{
				if (i->Sosed->belong(xk, yk) == true)
				{
					next = i->Sosed;
				}
			}
		}
		if (next == nullptr)
		{
			for (auto& k : now->Grans)
			{
				if (k->Sosed != nullptr)
				{
					for (auto& i : k->Sosed->Grans)
					{
						if (i->Sosed != nullptr)
						{
							if (i->Sosed->belong(xk, yk) == true)
							{
								next = i->Sosed;
							}
						}
					}
				}
			}
		}
	}

	if (next == nullptr)
	{
		for (auto& i : this->All_Cells_zero)
		{
			if (i->belong(xk, yk) == true)
			{
				next = i;
			}
		}
	}

	if (next == nullptr)
	{
		this->Smut.lock();
		cout << "Setka.cpp    " << "Fly_ex_Korol   " << "Long   find  " << endl;
		cout << "Setka.cpp    " << "Fly_ex_Korol  end " << xk << " " << yk << endl;
		cout << "Setka.cpp    " << "Fly_ex_Korol  start " << x_0 << " " << y_start << endl;
		this->Smut.unlock();
		//exit(-10);
		//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< x_0 << " " << sqrt(kvv(0.0, y_0, z_0)) << endl;
		//t_per = 100000000.0;
		//if (true)
		//{
		//	double A, B, C;
		//	for (auto& i : now->Grans)
		//	{
		//		double xx, yy;
		//		i->Get_Center(xx, yy);
		//		cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "centr = " << xx << " " << yy << endl;
		//		cout << "Setka.cpp    " << "Fly_ex_Korol   "<< t_per << endl;
		//		if (i->type == Axis)
		//		{
		//			continue;
		//		}
		//		A = i->aa;
		//		B = i->bb;
		//		C = i->cc;
		//		if (fabs(B) < 0.00001)
		//		{
		//			tt1 = (-C / A - x_0) / Vx;
		//			cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "TT1 = " << tt1 << endl;
		//			if (std::fpclassify(tt1) == FP_NORMAL && tt1 > 0.0)
		//			{
		//				if (t_per > tt1)
		//				{
		//					if (now->belong(x_0 + 1.02 * tt1 * Vx, sqrt(kv(y_0 + 1.02 * tt1 * Vy) + kv(z_0 + 1.02 * tt1 * Vz))) == false)
		//					{
		//						t_per = tt1;
		//						next_mb = i->Sosed;
		//					}
		//				}
		//			}
		//			continue;
		//		}
		//		cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "1 - " << A << " " << B << " " << C << endl;
		//		//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "2 - " << x_0 << " " << y_0 << " " << z_0 << endl;
		//		//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "3 - " << Vx << " " << Vy << " " << Vz << endl;
		//		a1 = 4.0 * kv(A * Vx * (C + A * x_0) - kv(B) * (Vy * y_0 + Vz * z_0)) - 4.0 * (kv(A * Vx) - //
		//			B * B * (kv(Vy) + kv(Vz))) * (kv(C + A * x_0) - B * B * (kv(y_0) + kv(z_0)));
		//		a2 = -A * Vx * (C + A * x_0) + B * B * (Vy * y_0 + Vz * z_0);
		//		a3 = A * A * Vx * Vx - B * B * (Vy * Vy + Vz * Vz);
		//		cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "4 - " << a1 << " " << a2 << " " << a3 << endl;
		//		if (fabs(a1) < 0.000001)
		//		{
		//			a1 = 0.0;
		//		}
		//		if (a1 >= 0.0)
		//		{
		//			tt1 = (a2 + 0.5 * sqrt(a1)) / a3;
		//			tt2 = (a2 - 0.5 * sqrt(a1)) / a3;
		//			//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< A << " " << B << " " << C << endl;
		//			cout << "Setka.cpp    " << "Fly_ex_Korol   "<< tt1 << " " << tt2 << endl;
		//			if (std::fpclassify(tt1) == FP_NORMAL && tt1 > 0.0)
		//			{
		//				if (t_per > tt1)
		//				{
		//					if (now->belong(x_0 + 1.03 * tt1 * Vx, sqrt(kv(y_0 + 1.03 * tt1 * Vy) + kv(z_0 + 1.03 * tt1 * Vz))) == false)
		//					{
		//						cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "A " << endl;
		//						t_per = tt1;
		//					}
		//				}
		//			}
		//			if (std::fpclassify(tt2) == FP_NORMAL && tt2 > 0.0)
		//			{
		//				if (t_per > tt2)
		//				{
		//					if (now->belong(x_0 + 1.03 * tt2 * Vx, sqrt(kv(y_0 + 1.03 * tt2 * Vy) + kv(z_0 + 1.03 * tt2 * Vz))) == false)
		//					{
		//						cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "B " << endl;
		//						t_per = tt2;
		//					}
		//				}
		//			}
		//		}
		//	}
		//	cout << "Setka.cpp    " << "Fly_ex_Korol   "<< t_per << endl;
		//	xk = x_0 + t_per * Vx;
		//	yk = sqrt(kv(y_0 + t_per * Vy) + kv(z_0 + t_per * Vz));
		//	cout << "Setka.cpp    " << "Fly_ex_Korol   "<< xk << " " << yk << endl;

		//	a = t_per * 0.99;
		//	b = t_per * 1.01;
		//}




		//exit(-1);
		for (auto& i : this->All_Cells)
		{
			if (i->belong(xk, yk) == true)
			{
				next = i;
				break;
			}
		}
	}

	if (next != nullptr)
	{
		if (AZ[i_z][i_alp] == false)
		{
			AZ[i_z][i_alp] = true;
			m_m.lock();
			Mu_statistic[area][i_z][i_alp] += mu/max(0.3 * SINKR[i_alp], sin(alpha));
			m_m.unlock();
		}

		if (next->zona > now->zona)
		{
			double Muuu = Mu_[area][next->zona][i_alp] * 0.3 * SINKR[i_alp] * kv(Ri[next->zona] / Rmax_);
			if (fabs(mu3) < Muuu)
			{
				if (fabs(mu3) >= adv.MakeRandom() * Muuu)
				{
					mu3 = Muuu * sign(mu3);
				}
				else
				{
					return;
				}
			}
		}

		//cout << "continue " << endl;
		Fly_exchenge_Imit_Korol(MK, sens, AZ, X, Y, Z, Vx, Vy, Vz, next, mu3, area, false, mu_start, to_I, to_J, georaschep, zon_stat);
		//cout << "continue END" << endl;
	}
	else
	{

		if (xk > -1000.0 / RR_ && xk < 1000.0 / RR_ && yk < 1600.0 / RR_)
		{
			cout << "Setka.cpp    " << "Fly_ex_Korol   " << "ERROR  8670  poteryal atom" << endl;
			cout << "Setka.cpp    " << "Fly_ex_Korol   " << "SDFSDF    " << xk << " " << yk << endl;
		}
	}

	return;
}

void Setka::Fly_exchenge_Imit_Korol_auto_weight(MKmethod& MK, int& s1, int& s2, int& s3, double x_0, double y_0, double z_0, double Vx, double Vy, double Vz, Cell* now, double mu, //
	int area, bool ExCh, const double& mu_start, int to_I, int to_J, bool georaschep, int zon_stat)
{
	// Новая версия полёта нейтрального атома. Функция без лишнего кода и насыщена комментариями
	// Функция использует и геометрическое и физическое расщепление!!!
	// Описание переменных:
	// sens - датчик случайных чисел
	// x_0 ... -  координаты начала полёта нейтрального атома в данной ячейке
	// Vx ... -   скорость летящего нейтрального атома
	// now -      текущая ячейка
	// mu -       вес нейтрального атома
	// area -     какой атом летит
	// mu_start - начальный вес атома при инициализации (чтобы понимать, сколько от исходного осталось)
	// to_I - в какую область предназначался атом
	//cout << "B" << endl;
	Cell* next_mb = nullptr;


	double normV = sqrt(kvv(Vx, Vy, Vz));
	double dt = now->L / normV;  // Время на один мини-шаг в ячейке (планируется сделать много шагов в ячейке)
	double ddt = 0.001 * dt;
	// Нахождение времени до попадания в следующую ячейку
	//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "V = " << Vx << " " << Vy << " " << Vz << endl;
	//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "A" << endl;
	double a = 0.0, b = dt;
	double t_per = 100000000000.0;
	// Блок проверки перигелия к оси
	/*double t_per = (-Vy * y_0 - Vz * z_0) / (kv(Vy) + kv(Vz));
	if (t_per > 0.0000001)
	{
		if (now->belong(x_0 + t_per * Vx, sqrt(kv(y_0 + t_per * Vy) + kv(z_0 + t_per * Vz))) == true)
		{
			t_per = -1.0;
		}
	}
	else
	{
		t_per = -1.0;
	}*/

	//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "A" << endl;
	double tt1, tt2;
	double a1, a2, a3;

	// Пытаемся определить время выхода точно
	if (true)
	{
		double A, B, C;
		for (auto& i : now->Grans)
		{
			if (i->type == Axis)
			{
				continue;
			}
			A = i->aa;
			B = i->bb;
			C = i->cc;

			/*if (fabs(B) < 0.001)
			{
				tt1 = (-C / A - x_0) / Vx;
				if (std::fpclassify(tt1) == FP_NORMAL && tt1 > 0.0)
				{
					if (t_per > tt1)
					{
						if (now->belong(x_0 + (tt1 + min(ddt, 0.03 * tt1)) * Vx, sqrt(kv(y_0 + (tt1 + min(ddt, 0.03 * tt1)) * Vy) + kv(z_0 + (tt1 + min(ddt, 0.03 * tt1)) * Vz))) == false)
						{
							t_per = tt1;
							next_mb = i->Sosed;
						}
					}
				}
				continue;
			}*/
			/*cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "1 - " << A << " " << B << " " << C << endl;
			cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "2 - " << x_0 << " " << y_0 << " " << z_0 << endl;
			cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "3 - " << Vx << " " << Vy << " " << Vz << endl;*/
			a1 = 4.0 * kv(A * Vx * (C + A * x_0) - kv(B) * (Vy * y_0 + Vz * z_0)) - 4.0 * (kv(A * Vx) - //
				B * B * (kv(Vy) + kv(Vz))) * (kv(C + A * x_0) - B * B * (kv(y_0) + kv(z_0)));
			a2 = -A * Vx * (C + A * x_0) + B * B * (Vy * y_0 + Vz * z_0);
			a3 = A * A * Vx * Vx - B * B * (Vy * Vy + Vz * Vz);
			//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "4 - " << a1 << " " << a2 << " " << a3 << endl;
			if (fabs(a1) < 0.0000001)
			{
				a1 = 0.0;
			}
			if (a1 >= 0.0)
			{
				tt1 = (a2 + 0.5 * sqrt(a1)) / a3;
				tt2 = (a2 - 0.5 * sqrt(a1)) / a3;
				//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< A << " " << B << " " << C << endl;
				if (std::fpclassify(tt1) == FP_NORMAL && tt1 > 0.0)
				{
					if (t_per > tt1)
					{
						if (now->belong(x_0 + (tt1 + min(ddt, 0.005 * tt1)) * Vx, sqrt(kv(y_0 + (tt1 + min(ddt, 0.005 * tt1)) * Vy) + kv(z_0 + (tt1 + min(ddt, 0.01 * tt1)) * Vz))) == false)
						{
							t_per = tt1;
							next_mb = i->Sosed;
						}
					}
				}
				if (std::fpclassify(tt2) == FP_NORMAL && tt2 > 0.0)
				{
					if (t_per > tt2)
					{
						if (now->belong(x_0 + (tt2 + min(ddt, 0.005 * tt2)) * Vx, sqrt(kv(y_0 + (tt2 + min(ddt, 0.005 * tt2)) * Vy) + kv(z_0 + (tt2 + min(ddt, 0.01 * tt2)) * Vz))) == false)
						{
							t_per = tt2;
							next_mb = i->Sosed;
						}
					}
				}
			}
		}

		a = t_per - min(ddt, 0.005 * t_per);
		b = t_per + min(ddt, 0.005 * t_per);
	}

	bool hand = false;
	if (now->belong(x_0 + a * Vx, sqrt(kv(y_0 + a * Vy) + kv(z_0 + a * Vz))) == false)
	{
		hand = true;
	}
	else if (now->belong(x_0 + b * Vx, sqrt(kv(y_0 + b * Vy) + kv(z_0 + b * Vz))) == true)
	{
		hand = true;
	}
	else if (t_per > 100000000.0)
	{
		//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "ERROR  8568  " << t_per << endl;
		hand = true;
	}

	double k;

	// Находим время  в ячейке
	int lk = 0;
	if (hand)
	{
		a = 0.0;
		b = dt;
		t_per = 100000000000.0;
		for (auto& i : now->Grans)
		{
			a1 = -kv(i->a * Vx) * (kv(i->a * Vx) - kv(Vy) - kv(Vz)) * kv(Vz * y_0 - Vy * z_0);
			a2 = (kv(i->a * Vx) - kv(Vy) - kv(Vz)) * (Vz * z_0 + Vy * y_0);
			a3 = (kv(Vy) + kv(Vz)) * (-kv(i->a * Vx) + kv(Vy) + kv(Vz));
			if (a1 >= 0)
			{
				tt1 = (sqrt(a1) + a2) / a3;
				tt2 = (-sqrt(a1) + a2) / a3;
				if (std::fpclassify(tt1) == FP_NORMAL && tt1 > 0.0)
				{
					if (now->belong(x_0 + tt1 * Vx, sqrt(kv(y_0 + tt1 * Vy) + kv(z_0 + tt1 * Vz))) == false)
					{
						if (t_per > tt1)
						{
							t_per = tt1;
						}
					}
				}
				if (std::fpclassify(tt2) == FP_NORMAL && tt2 > 0.0)
				{
					if (now->belong(x_0 + tt2 * Vx, sqrt(kv(y_0 + tt2 * Vy) + kv(z_0 + tt2 * Vz))) == false)
					{
						if (t_per > tt2)
						{
							t_per = tt2;
						}
					}
				}
			}
		}

		while (now->belong(x_0 + b * Vx, sqrt(kv(y_0 + b * Vy) + kv(z_0 + b * Vz))) == true) // Слишком маленький шаг
		{
			lk++;
			if (lk > 1000)
			{
				cout << "Setka.cpp    " << "Fly_ex_Korol   " << "8055 dt = " << dt << " " << x_0 << " " << sqrt(kv(y_0) + kv(z_0)) << endl;
				cout << "Setka.cpp    " << "Fly_ex_Korol   " << normV << " " << Vx << " " << Vy << " " << Vz << endl;
				exit(-1);
			}
			b = b * 1.1;
		}

		if (t_per < b)
		{
			b = t_per;
		}

		k = (a + b) / 2.0;
		while (fabs(a - b) * normV > now->L / 300.0)
		{
			if (now->belong(x_0 + k * Vx, sqrt(kv(y_0 + k * Vy) + kv(z_0 + k * Vz))) == true)
			{
				a = k;
			}
			else
			{
				b = k;
			}
			k = (a + b) / 2.0;
		}
	}

	double time = b; // Время нахождения в ячейке


	bool slay = false;
	double t1, t2;
	//Нужно проверить пересечение со сферой
	double D = kv(2.0 * Vx * x_0 + 2.0 * Vy * y_0 + 2.0 * Vz * z_0) - 4.0 * kvv(Vx, Vy, Vz) * (kvv(x_0, y_0, z_0) - kv(Rmax_));
	if (D >= 0.0)
	{
		t1 = (-2.0 * (Vx * x_0 + Vy * y_0 + Vz * z_0) - sqrt(D)) / (2.0 * kv(normV));
		t2 = (-2.0 * (Vx * x_0 + Vy * y_0 + Vz * z_0) + sqrt(D)) / (2.0 * kv(normV));
		if (t1 > 0.0 && t1 < time)
		{
			if (x_0 + Vx * t1 >= 0.0)
			{
				time = t1;
				a = time * 0.999;
				b = time * 1.001;
				slay = true;
			}
		}
		if (t2 > 0.0 && t2 < time)
		{
			if (x_0 + Vx * t2 >= 0.0)  // Для правой полусферы
			{
				time = t2;
				a = time * 0.999;
				b = time * 1.001;
				slay = true;
			}
		}
	}
	// Теперь пересечение с правой границей
	if (fabs(Vx) > 0.000001)
	{
		t1 = -x_0 / Vx;
		if (t1 > 0.0 && t1 <= time && kvv(0.0, y_0 + t1 * Vy, z_0 + t1 * Vz) >= kv(Rmax_))
		{
			time = t1;
			a = time * 0.999;
			b = time * 1.001;
			slay = true;
		}
	}


	// теперь время time, a, b уже определены
	// a - последний момент времени в ячейке, b - первый момент за пределами ячейки


	//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< a << " " << b << endl;
	double uz, uz_M, uz_E;								// средние скорости в интегралах
	double cp = sqrt(now->par[0].p / now->par[0].ro);
	double vx = now->par[0].u;							// Скорости плазмы в ячейке
	double vy = now->par[0].v;
	double ro = now->par[0].ro;


	double t_ex = 0.0;									// время до перезарядки
	// time - время нахождения атома в ячейке
	t2 = time;									// время после перезарядки (будет ниже)

	double x_ex = x_0, y_ex = y_0, z_ex = z_0;		    // координаты перезарядки
	double u1, u2, u3;

	double y_start = sqrt(kv(y_0) + kv(z_0));

	double l = sqrt(kvv(time * Vx, time * Vy, time * Vz));                    // Расстояние до выхода из ячейки

	double alpha = 0.0;								    // угол 
	double u = 0.0;                                     // модуль относительной скорости атома и плазмы
	double nu_ex = 0.0;								    // частота перезарядки
	double mu_ex = 0.0;									// вес перезаряженного атома
	double mu2 = mu;									// вес не-перезаряженного атома

	double drob = 5.0;     // Для того чтобы точнее учесть проворот скорости во время полёта

	int i_z = 0, i_alp = 0;  // Области перезярядки по геометрическим зонам

	if (true)//(ExCh == false) // Если перезарядки в этой ячейке ещё не было
	{
		double kappa = 0.0;
		//drob = min(fabs(polar_angle(y_0, z_0) - polar_angle(y_0 + time * Vy, z_0 + time * Vz)) / 0.017 + 5.0, 80.0);

		if (true)
		{

			for (double tt = 0.0; tt < time * 0.98; tt += time / drob)
			{
				double yy = sqrt(kv(y_0 + (tt + time / (2.0 * drob)) * Vy) + kv(z_0 + (tt + time / (2.0 * drob)) * Vz));

				alpha = polar_angle(y_0 + (tt + time / (2.0 * drob)) * Vy, z_0 + (tt + time / (2.0 * drob)) * Vz);
				u = sqrt(kvv(Vx - vx, Vy - vy * cos(alpha), Vz - vy * sin(alpha)));


				if (u / cp > 7.0)
				{
					uz = Velosity_1(u, cp);
					nu_ex = (ro * uz * sigma(uz)) / Kn_;
					//cout << x_0 << " " << y_start << " " << u / cp << " " << u << " " << cp <<  endl;
				}
				else
				{
					nu_ex = (ro * MK.int_1(u, cp)) / Kn_;        // Пробуем вычислять интеграллы численно

					//uz = Velosity_1(u, cp);
					//cout << (uz * sigma(uz)) << " " << MK.int_1(u, cp) << " " << u << " " << cp << endl;
				}
				kappa += (nu_ex * time / drob);
			}
		}
		else
		{
			alpha = polar_angle(y_0, z_0);
			u = sqrt(kvv(Vx - vx, Vy - vy * cos(alpha), Vz - vy * sin(alpha)));
			if (u / cp > 7.0)
			{
				uz = Velosity_1(u, cp);
				nu_ex = (ro * uz * sigma(uz)) / Kn_;
			}
			else
			{
				nu_ex = (ro * MK.int_1(u, cp)) / Kn_;        // Пробуем вычислять интеграллы численно
			}
			kappa = (nu_ex * time);
		}


		t_ex = -(time / kappa) * log(1.0 - MyRandom(s1, s2, s3) * (1.0 - exp(-kappa))); // Время до перезарядки
		t2 = time - t_ex;   // Время сколько лететь после того, как атом перезарядился
		mu_ex = mu * (1.0 - exp(-kappa)); // вес перезаряженного атома
		mu2 = mu - mu_ex;                 // вес оставшегося неперезаряженного атома

		x_ex = x_0 + t_ex * Vx;  // Координаты перезарядки
		y_ex = y_0 + t_ex * Vy;
		z_ex = z_0 + t_ex * Vz;

		alpha = polar_angle(y_ex, z_ex);

		// проверка, если перезарядка произошла за пределами ячейки (нужно немного сдвигать её в этом случае)
		if (now->belong(x_ex, sqrt(kv(y_ex) + kv(z_ex))) == false)
		{
			x_ex = x_0 + a * Vx;  // Координаты перезарядки
			y_ex = y_0 + a * Vy;
			z_ex = z_0 + a * Vz;
			t_ex = a;
			t2 = time - t_ex;
			alpha = polar_angle(y_ex, z_ex);
		}

		double r = sqrt(kvv(x_ex, y_ex, z_ex));
		double all = polar_angle(x_ex, sqrt(kv(y_ex) + kv(z_ex)));

		// Находим пересечение с лучами зрения
		if (false)
		{
			double xk = now->x_center;
			double yk = now->y_center;
			double alf = now->alf_center;
			double aa = tan(alf);


			double nn = sqrt(kvv(x_ex, y_ex, z_ex));
			double e1 = x_ex / nn;
			double e2 = y_ex / nn;
			double e3 = z_ex / nn;
			
			double Vu = Vx * e1 + Vy * e2 + Vz * e3;

			now->pogloshenie[area][min(pogl_rad_ - 1, max(0 , (int)( (Vu - pogVmin) / ( (pogVmax - pogVmin) / pogl_rad_) )))] += t_ex * mu_ex + mu2 * time;
		}

		// для сбора статистики
		if (func_stat)
		{
			if (r < 100.0 / RR_ && r > 60.0 / RR_)
			{
				bool sec = false;
				bool sec2 = false;
				double tt;
				double a1 = Vx * x_0 + Vy * y_0 + Vz * z_0;
				double a2 = kvv(Vx, Vy, Vz);
				double a3 = kvv(x_0, y_0, z_0) - kv(R_stat / RR_);
				if (4.0 * a1 * a1 - 4.0 * a2 * a3 >= 0.0)
				{
					double t1 = (-2.0 * a1 + sqrt(4.0 * a1 * a1 - 4.0 * a2 * a3)) / (2.0 * a2);
					double t2 = (-2.0 * a1 - sqrt(4.0 * a1 * a1 - 4.0 * a2 * a3)) / (2.0 * a2);
					if (t1 > 0.000000001 && t1 < time)
					{
						sec = true;
						tt = t1;
						if (t1 < t_ex)
						{
							sec2 = true; // пролетел до перезярядки
						}
						else
						{
							sec2 = false;
						}
					}
					if (t2 > 0.000000001 && t2 < time)
					{
						sec = true;
						tt = t2;
						if (t2 < t_ex)
						{
							sec2 = true; // пролетел до перезярядки
						}
						else
						{
							sec2 = false;
						}
					}
				}

				if (sec == true)
				{
					double xx, yy, zz;
					xx = x_0 + tt * Vx;
					yy = y_0 + tt * Vy;
					zz = z_0 + tt * Vz;
					double Vr, Vphi, Vthe;
					double phi_ = polar_angle(yy, zz);  //polar_angle(xx, sqrt(kv(yy) + kv(zz)));
					double the_ = polar_angle(xx, sqrt(kv(yy) + kv(zz)));
					spherical_skorost(yy, zz, xx, Vy, Vz, Vx, Vr, Vphi, Vthe);

					double Vn = fabs((Vx * xx + Vy * yy + Vz * zz) / (R_stat / RR_));
					//double Vn = sqrt(kvv(Vx, Vy, Vz));
					double Vrr = (Vy * cos(phi_) + Vz * sin(phi_));
					double mumu = 0.0;
					if (sec2)
					{
						mumu = mu / Vn;
					}
					else
					{
						mumu = mu2 / Vn;
					}

					mut_mom.lock();
					this->mu_mom[area][(int)(the_ / (pi_ / Al_stat))] += mumu;
					this->Vx_mom[area][(int)(the_ / (pi_ / Al_stat))] += mumu * Vx;
					this->Vy_mom[area][(int)(the_ / (pi_ / Al_stat))] += mumu * Vrr;
					this->Vxx_mom[area][(int)(the_ / (pi_ / Al_stat))] += mumu * Vx * Vx;
					this->Vyy_mom[area][(int)(the_ / (pi_ / Al_stat))] += mumu * Vrr * Vrr;
					this->Vxy_mom[area][(int)(the_ / (pi_ / Al_stat))] += mumu * Vx * Vrr;
					this->Vxxx_mom[area][(int)(the_ / (pi_ / Al_stat))] += mumu * Vx * Vx * Vx;
					this->T_mom[area][(int)(the_ / (pi_ / Al_stat))] += mumu * kvv(Vx, Vy, Vz);
					mut_mom.unlock();

					
				}
			}
		}

		double rr = sqrt(kvv(0.0, y_ex, z_ex));


		// Считаем источники для частицы до перезарядки

		u = sqrt(kvv(Vx - vx, Vy - vy * cos(alpha), Vz - vy * sin(alpha)));
		u1 = vx - Vx;
		u2 = vy * cos(alpha) - Vy;
		u3 = vy * sin(alpha) - Vz;
		double skalar = Vx * u1 + Vy * u2 + Vz * u3;

		// определим область перезарядки (по геометрическим зонам)
		i_z = geo_zones(r);
		i_alp = alpha_zones(x_ex, rr);

		int area2 = 0;
		// определим область в которой находится атом сейчас (это параметр самой ячейки)
		if (now->type == C_5)
		{
			area2 = 3;
		}
		else if (now->type == C_4)
		{
			area2 = 2;
		}
		else if (now->type == C_1 || now->type == C_centr)
		{
			area2 = 0;
		}
		else
		{
			area2 = 1;
		}


		now->mut.lock();

		// Блок расчёта функции распределения, если это нужно
		if (area == 3)
		{
			if (now->df_s4_bool == true)
			{
				now->df_s4->Add_point(Vx, Vy, Vz, y_ex, z_ex, t_ex * mu);
			}
		}
		else if (area == 2)
		{
			if (now->df_s3_bool == true)
			{
				now->df_s3->Add_point(Vx, Vy, Vz, y_ex, z_ex, t_ex * mu);
			}
		}
		else if (area == 1)
		{
			if (now->df_s2_bool == true)
			{
				now->df_s2->Add_point(Vx, Vy, Vz, y_ex, z_ex, t_ex * mu);
			}
		}
		else if (area == 0)
		{
			if (now->df_s1_bool == true)
			{
				now->df_s1->Add_point(Vx, Vy, Vz, y_ex, z_ex, t_ex * mu);
			}
		}

		//now->par[0].w_m[area] += mu / max(sin(alpha), 0.3 * Sinus[i_alp]);
		now->par[0].F_n += t_ex * mu;
		now->par[0].F_u += t_ex * Vx * mu;
		now->par[0].F_T += t_ex * kvv(Vx, Vy, Vz) * mu;
		now->par[0].H_n[area] += t_ex * mu;
		now->par[0].H_u[area] += t_ex * Vx * mu;
		now->par[0].H_T[area] += t_ex * kvv(Vx, Vy, Vz) * mu;

		double Vrr = (Vy * cos(alpha) + Vz * sin(alpha));
		now->par[0].H_uu[area] += t_ex * Vx * Vx * mu;
		now->par[0].H_uv[area] += t_ex * Vrr * Vx * mu;
		now->par[0].H_vv[area] += t_ex * Vrr * Vrr * mu;
		now->par[0].H_uuu[area] += t_ex * Vx * Vx * Vx * mu;

		if (u / cp > 7.0)
		{
			uz = Velosity_1(u, cp);
			uz_M = Velosity_2(u, cp) / (uz * kv(cp) * cp * pi * sqrtpi_);
			uz_E = Velosity_3(u, cp);

			now->par[0].I_u += -mu_ex * uz_M * u1 / u;
			now->par[0].I_v += -mu_ex * uz_M * (u2 * cos(alpha) + u3 * sin(alpha)) / u;
			now->par[0].I_T += mu_ex * (-0.25 * (3.0 * kv(cp) + 2.0 * kv(u)) * (uz_E / uz) - //
				uz_M * skalar / u);
		}
		else
		{
			double k1 = MK.int_1(u, cp);
			double k2 = MK.int_2(u, cp);
			double k3 = MK.int_3(u, cp);
			now->par[0].I_u += mu_ex * (k2 / k1) * u1 / u;
			now->par[0].I_v += mu_ex * (k2 / k1) * (u2 * cos(alpha) + u3 * sin(alpha)) / u;
			now->par[0].I_T += mu_ex * (-0.5 * k3 / k1 + k2 / k1 * skalar / u);
		}

		//drob = min(fabs(polar_angle(y_0, z_0) - polar_angle(y_0 + t_ex * Vy, z_0 + t_ex * Vz)) / 0.017 + 5.0, 80.0);

		if (true)
		{
			for (double tt = 0.0; tt < t_ex * 0.99; tt += t_ex / drob)
			{
				alpha = polar_angle(y_0 + (tt + t_ex / (2.0 * drob)) * Vy, z_0 + (tt + t_ex / (2.0 * drob)) * Vz);

				now->par[0].F_v += (t_ex / drob) * (Vy * cos(alpha) + Vz * sin(alpha)) * mu;
				now->par[0].H_v[area] += (t_ex / drob) * (Vy * cos(alpha) + Vz * sin(alpha)) * mu;
			}
		}
		else
		{
			now->par[0].F_v += t_ex * (Vy * cos(alpha) + Vz * sin(alpha)) * mu;
			now->par[0].H_v[area] += t_ex * (Vy * cos(alpha) + Vz * sin(alpha)) * mu;
		}

		//now->par[0].I_u += -mu_ex * uz_M * (sigma2(uz_M, cp) / sigma2(uz, cp)) * u1 / u;
		//now->par[0].I_v += -mu_ex * uz_M * (sigma2(uz_M, cp) / sigma2(uz, cp)) * (u2 * cos(alpha) + u3 * sin(alpha)) / u;
		//now->par[0].I_T += mu_ex * (-0.25 * (3.0 * kv(cp) + 2.0 * kv(u)) * (uz_E / uz) * (sigma2(uz_E, cp) / sigma2(uz, cp)) - //
		//	uz_M * (sigma2(uz_M, cp) / sigma2(uz, cp)) * skalar / u);

		/*now->par[0].II_u += -mu_ex * u1;
		now->par[0].II_v += -mu_ex * (u2 * cos(alpha) + u3 * sin(alpha));
		now->par[0].II_T += mu_ex * 0.5 * (kvv(Vx, Vy, Vz) - kvv(vx, vy, 0.0));*/

		now->mut.unlock();


		double Ur, Uthe, Uphi, Vr, Vthe, Vphi;
		spherical_skorost(y_ex, z_ex, x_ex, vy * cos(alpha), vy * sin(alpha), vx, Ur, Uphi, Uthe);
		spherical_skorost(y_ex, z_ex, x_ex, Vy, Vz, Vx, Vr, Vphi, Vthe);
		double mu3;


		int stat_zone_ = now->zona;

		if (area != area2)
		{
			stat_zone_ = -1;
		}


		if (area2 == 0 || Ur / cp > 3.0)   // Без геометрического расщепления
		{
		aa:
			bool kj = true;
			mu3 = mu_ex;
			double alpha = polar_angle(y_ex, z_ex);
			double aa, bb, cc;  // Декартовы скорости после перезарядки
			vector <double> Wr(1);
			vector <double> Wthe(1);
			vector <double> Wphi(1);
			vector <double> mu_(1);
			MK.Change_Velosity4(s1, s2, s3, Ur / cp, Uthe / cp, Uphi / cp, Vr / cp, Vthe / cp, Vphi / cp, Wr, Wthe, Wphi, mu_, cp, r, 0);
			Wr[0] = Wr[0] * cp;
			Wthe[0] = Wthe[0] * cp;
			Wphi[0] = Wphi[0] * cp;
			dekard_skorost(y_ex, z_ex, x_ex, Wr[0], Wphi[0], Wthe[0], bb, cc, aa);
			double peregel, time_do_peregel;
			int ii_z;     // В какую зону летит родившийся атом
			int ii_alp;   // В какую зону летит родившийся атом
			int to_ii;
			double r_per = r;  // радиус в перегелии

			if (Wr[0] < 0.0)
			{
				r_per = this->distination(x_ex, y_ex, z_ex, aa, bb, cc, area2, to_ii, ii_z, ii_alp);
			}
			else
			{
				to_ii = i_z;
				ii_z = i_z;
				ii_alp = i_alp;
			}

			double Muuu = max(min(Mu_[area2][ii_z][ii_alp] * 0.3 * SINKR[ii_alp] * kv(r_per / Rmax_),
				Mu_[area2][i_z][i_alp] * 0.3 * SINKR[i_alp] * kv(r / Rmax_)), MinWes);

			if (kj == true)
			{
				if (fabs(mu3) >= Muuu)  // pred(area2, ii_z, all))
				{
					kj = true;
				}
				else
				{
					if (fabs(mu3) >= MyRandom(s1, s2, s3) * Muuu)  // pred(area2, ii_z, all))
					{
						mu3 = Muuu * sign(mu3);  // pred(area2, ii_z, all) * sign(mu3);
						kj = true;
					}
					else
					{
						kj = false;
					}
				}
			}

			if (kj == true)
			{
				now->mut.lock();
				now->par[0].II_u += mu3 * (Vx - aa);
				now->par[0].II_v += mu3 * ((Vy - bb) * cos(alpha) + (Vz - cc) * sin(alpha));
				now->par[0].II_T += mu3 * 0.5 * (kvv(Vx, Vy, Vz) - kvv(aa, bb, cc));
				now->mut.unlock();

				Fly_exchenge_Imit_Korol_auto_weight(MK, s1,s2,s3, x_ex, y_ex, z_ex, aa, bb, cc, now, mu3, area2, false, mu_start, to_ii, ii_alp, true, stat_zone_);
			}
			//this->Fly_exchenge_Imit(MK, sens, x_ex, y_ex, z_ex, aa, bb, cc, now, mu3, -log(1.0 - sens->MakeRandom()),//
			//	0.0, area2, mu_start, to_I, 0);

		}
		else  // делаем геометрическое расщепление
		{
			int I = geo_zones(r, 1.2);  // Число дополнительных траекторий
			//if (I > to_I)  // Для того чтобы не расщеплять атом, который летел в to_I область на атомы из более крупных областей
			//{
			//	I = to_I;
			//}

			vector <double> Wr(I + 1);
			vector <double> Wthe(I + 1);
			vector <double> Wphi(I + 1);
			vector <double> mu_(I + 1);

			//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< " C " << endl;
			bool bbb = MK.Change_Velosity4(s1, s2, s3, Ur / cp, Uthe / cp, Uphi / cp, Vr / cp, Vthe / cp, Vphi / cp, Wr, Wthe, Wphi, mu_, cp, r, I);
			//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< " D " << endl;

			for (int i = 0; i <= I; i++)
			{
				Wr[i] = Wr[i] * cp;
				Wthe[i] = Wthe[i] * cp;
				Wphi[i] = Wphi[i] * cp;
			}

			double aa, bb, cc, time_do_peregel, peregel;
			int ii_z, ii_alp;
			int to_ii;   // В какую зону по радиусу направлен
			bool kj = true;
			double r_per = r;
			//now->par[0].num_atoms++;
			for (int i = 0; i < I; i++)
			{
				dekard_skorost(y_ex, z_ex, x_ex, Wr[i], Wphi[i], Wthe[i], bb, cc, aa);
				mu3 = mu_[i] * mu_ex;
				r_per = this->distination(x_ex, y_ex, z_ex, aa, bb, cc, area2, to_ii, ii_z, ii_alp);

				double Muuu = max(min(Mu_[area2][ii_z][ii_alp] * 0.3 * SINKR[ii_alp] * kv(r_per / Rmax_),
					Mu_[area2][i_z][i_alp] * 0.3 * SINKR[i_alp] * kv(r / Rmax_)), MinWes);

				if (fabs(mu3) >= Muuu)
				{
					kj = true;
				}
				else
				{
					if (fabs(mu3) >= MyRandom(s1, s2, s3) * Muuu)
					{
						mu3 = Muuu;
						kj = true;
					}
					else
					{
						kj = false;
					}
				}
				if (kj == true)
				{
					now->mut.lock();
					now->par[0].II_u += mu3 * (Vx - aa);
					now->par[0].II_v += mu3 * ((Vy - bb) * cos(alpha) + (Vz - cc) * sin(alpha));
					now->par[0].II_T += mu3 * 0.5 * (kvv(Vx, Vy, Vz) - kvv(aa, bb, cc));
					now->mut.unlock();

					Fly_exchenge_Imit_Korol_auto_weight(MK, s1, s2, s3, x_ex, y_ex, z_ex, aa, bb, cc, now, mu3, area2, false, mu_start, to_ii, ii_alp, true, stat_zone_);
				}
			}

			if (bbb == true)  // Если вообще нужно запускать основной атом
			{
				mu3 = mu_[I] * mu_ex;

				// Нужно определить номер основного атома (от его перегелия)
				dekard_skorost(y_ex, z_ex, x_ex, Wr[I], Wphi[I], Wthe[I], bb, cc, aa);
				double r_per = r;
				if (Wr[I] < 0.0)
				{
					r_per = this->distination(x_ex, y_ex, z_ex, aa, bb, cc, area2, to_ii, ii_z, ii_alp);
				}
				else
				{
					to_ii = i_z;
					ii_z = i_z;
					ii_alp = i_alp;
				}

				double Muuu = max(min(Mu_[area2][ii_z][ii_alp] * 0.3 * SINKR[ii_alp] * kv(r_per / Rmax_),
					Mu_[area2][i_z][i_alp] * 0.3 * SINKR[i_alp] * kv(r / Rmax_)), MinWes);

				if (fabs(mu3) >= Muuu)
				{
					kj = true;
				}
				else
				{
					if (fabs(mu3) >= MyRandom(s1, s2, s3) * Muuu)
					{
						mu3 = Muuu * sign(mu3);
						kj = true;
					}
					else
					{
						kj = false;
					}
				}

				if (kj == true)
				{
					now->mut.lock();
					now->par[0].II_u += mu3 * (Vx - aa);
					now->par[0].II_v += mu3 * ((Vy - bb) * cos(alpha) + (Vz - cc) * sin(alpha));
					now->par[0].II_T += mu3 * 0.5 * (kvv(Vx, Vy, Vz) - kvv(aa, bb, cc));
					now->mut.unlock();

					Fly_exchenge_Imit_Korol_auto_weight(MK, s1, s2, s3, x_ex, y_ex, z_ex, aa, bb, cc, now, mu3, area2,
						true, mu_start, ii_z, ii_alp, false, stat_zone_);
				}
			}

		}


		// Здесь можно добавить рулетку для оставшегося неперезаряженного атома
		double Muuu = max(min(Mu_[area2][to_I][to_J] * 0.3 * SINKR[to_J] * kv(Ri[to_I] / Rmax_),
			Mu_[area2][i_z][i_alp] * 0.3 * SINKR[i_alp] * kv(r / Rmax_)), MinWes);

		if (fabs(mu2) < Muuu)
		{
			if (fabs(mu2) >= MyRandom(s1, s2, s3) * Muuu)
			{
				mu2 = Muuu * sign(mu2);
			}
			else
			{
				return;
			}
		}

	}




	// этот блок кода будет работать и для частиц, которые не попали в предыдущий цикл
	// здесь t2, _ex будут иметь исходные значения, как при инициализации
	alpha = polar_angle(y_ex + 0.5 * t2 * Vy, z_ex + 0.5 * t2 * Vz);
	now->mut.lock();

	// Блок расчёта функции распределения, если это нужно
	if (area == 3)
	{
		if (now->df_s4_bool == true)
		{
			now->df_s4->Add_point(Vx, Vy, Vz, y_ex, z_ex, t2 * mu2);
		}
	}
	else if (area == 2)
	{
		if (now->df_s3_bool == true)
		{
			now->df_s3->Add_point(Vx, Vy, Vz, y_ex, z_ex, t2 * mu2);
		}
	}
	else if (area == 1)
	{
		if (now->df_s2_bool == true)
		{
			now->df_s2->Add_point(Vx, Vy, Vz, y_ex, z_ex, t2 * mu2);
		}
	}
	else if (area == 0)
	{
		if (now->df_s1_bool == true)
		{
			now->df_s1->Add_point(Vx, Vy, Vz, y_ex, z_ex, t2 * mu2);
		}
	}

	now->par[0].F_n += t2 * mu2;
	now->par[0].F_u += t2 * Vx * mu2;
	now->par[0].F_T += t2 * kvv(Vx, Vy, Vz) * mu2;
	now->par[0].H_n[area] += t2 * mu2;
	now->par[0].H_u[area] += t2 * Vx * mu2;
	now->par[0].H_T[area] += t2 * kvv(Vx, Vy, Vz) * mu2;

	double Vrr = (Vy * cos(alpha) + Vz * sin(alpha));
	now->par[0].H_uu[area] += t2 * Vx * Vx * mu2;
	now->par[0].H_uv[area] += t2 * Vrr * Vx * mu2;
	now->par[0].H_vv[area] += t2 * Vrr * Vrr * mu2;
	now->par[0].H_uuu[area] += t2 * Vx * Vx * Vx * mu2;

	//drob = min(fabs(polar_angle(y_ex, z_ex) - polar_angle(y_ex + t2 * Vy, z_ex + t2 * Vz)) / 0.017 + 5.0, 80.0);

	if (true)
	{
		for (double tt = 0.0; tt < t2 * 0.99; tt += t2 / drob)  // было 20 и норм работало
		{
			alpha = polar_angle(y_ex + (tt + t2 / (2.0 * drob)) * Vy, z_ex + (tt + t2 / (2.0 * drob)) * Vz);

			now->par[0].F_v += (t2 / drob) * (Vy * cos(alpha) + Vz * sin(alpha)) * mu2;
			now->par[0].H_v[area] += (t2 / drob) * (Vy * cos(alpha) + Vz * sin(alpha)) * mu2;
		}
	}
	else
	{
		now->par[0].F_v += t2 * (Vy * cos(alpha) + Vz * sin(alpha)) * mu2;
		now->par[0].H_v[area] += t2 * (Vy * cos(alpha) + Vz * sin(alpha)) * mu2;
	}
	now->mut.unlock();

	double mu3 = mu2;
	/*double X = x_ex + t2 * Vx;
	double Y = y_ex + t2 * Vy;
	double Z = z_ex + t2 * Vz;*/

	bool goto1 = false;
aa1:

	double X = x_0 + b * Vx;
	double Y = y_0 + b * Vy;
	double Z = z_0 + b * Vz;

	if (goto1 == true && X < 2.0 && X > -3.0)
	{
		if (Vx > 0)
		{
			X = X + 0.00001;
		}
		else
		{
			X = X - 0.00001;
		}
	}

	Cell* next = nullptr;             // В какую ячейку попадаем следующую
	double xk, yk;
	xk = x_0 + b * Vx;
	yk = sqrt(kv(y_0 + b * Vy) + kv(z_0 + b * Vz));
	double r_k = sqrt(kv(xk) + kv(yk));
	double r_0 = sqrt(kv(x_0) + kv(y_0) + kv(z_0));

	/*if (yk < 5.0)
	{
		cout << "Setka.cpp    " << "Fly_ex_Korol   "<< b << "  " << xk << endl;
	}*/

	if (slay == true)  //Если надо убить траекторию
	{
		return;
	}

	if (xk < Left_ + 2.0 / RR_)  // ЭТОГО НЕ БЫЛО
	{
		return;
	}

	if (yk > R5_ - 2.0 / RR_)  // ЭТОГО НЕ БЫЛО
	{
		return;
	}

	if (r_k >= Rmax_ - 3.0 / RR_ && xk >= -0.5 / RR_)
	{
		return;
		next = nullptr;
	}

	if (next_mb != nullptr)
	{
		if (next_mb->belong(xk, yk) == true)
		{
			next = next_mb;
		}
	}

	if (next == nullptr)
	{
		for (auto& i : now->Grans)
		{
			if (i->Sosed != nullptr)
			{
				if (i->Sosed->belong(xk, yk) == true)
				{
					next = i->Sosed;
				}
			}
		}
		if (next == nullptr)
		{
			for (auto& k : now->Grans)
			{
				if (k->Sosed != nullptr)
				{
					for (auto& i : k->Sosed->Grans)
					{
						if (i->Sosed != nullptr)
						{
							if (i->Sosed->belong(xk, yk) == true)
							{
								next = i->Sosed;
							}
						}
					}
				}
			}
		}
	}

	if (next == nullptr)
	{
		for (auto& i : this->All_Cells_zero)
		{
			if (i->belong(xk, yk) == true)
			{
				next = i;
			}
		}
	}

	if (next == nullptr)
	{
		if (goto1 == false)
		{
			goto1 = true;
			b = b * 1.01;
			b = b + 0.00001;
			goto aa1;
		}
		this->Smut.lock();
		cout << "Setka.cpp    " << "1ds Fly_ex_Korol   " << "Long   find  " << endl;
		cout << "Setka.cpp    " << "1ds Fly_ex_Korol  end " << xk << " " << yk << endl;
		cout << "Setka.cpp    " << "1ds Fly_ex_Korol  start " << x_0 << " " << y_start << endl;
		this->Smut.unlock();
		//exit(-10);
		//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< x_0 << " " << sqrt(kvv(0.0, y_0, z_0)) << endl;
		//t_per = 100000000.0;
		//if (true)
		//{
		//	double A, B, C;
		//	for (auto& i : now->Grans)
		//	{
		//		double xx, yy;
		//		i->Get_Center(xx, yy);
		//		cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "centr = " << xx << " " << yy << endl;
		//		cout << "Setka.cpp    " << "Fly_ex_Korol   "<< t_per << endl;
		//		if (i->type == Axis)
		//		{
		//			continue;
		//		}
		//		A = i->aa;
		//		B = i->bb;
		//		C = i->cc;
		//		if (fabs(B) < 0.00001)
		//		{
		//			tt1 = (-C / A - x_0) / Vx;
		//			cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "TT1 = " << tt1 << endl;
		//			if (std::fpclassify(tt1) == FP_NORMAL && tt1 > 0.0)
		//			{
		//				if (t_per > tt1)
		//				{
		//					if (now->belong(x_0 + 1.02 * tt1 * Vx, sqrt(kv(y_0 + 1.02 * tt1 * Vy) + kv(z_0 + 1.02 * tt1 * Vz))) == false)
		//					{
		//						t_per = tt1;
		//						next_mb = i->Sosed;
		//					}
		//				}
		//			}
		//			continue;
		//		}
		//		cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "1 - " << A << " " << B << " " << C << endl;
		//		//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "2 - " << x_0 << " " << y_0 << " " << z_0 << endl;
		//		//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "3 - " << Vx << " " << Vy << " " << Vz << endl;
		//		a1 = 4.0 * kv(A * Vx * (C + A * x_0) - kv(B) * (Vy * y_0 + Vz * z_0)) - 4.0 * (kv(A * Vx) - //
		//			B * B * (kv(Vy) + kv(Vz))) * (kv(C + A * x_0) - B * B * (kv(y_0) + kv(z_0)));
		//		a2 = -A * Vx * (C + A * x_0) + B * B * (Vy * y_0 + Vz * z_0);
		//		a3 = A * A * Vx * Vx - B * B * (Vy * Vy + Vz * Vz);
		//		cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "4 - " << a1 << " " << a2 << " " << a3 << endl;
		//		if (fabs(a1) < 0.000001)
		//		{
		//			a1 = 0.0;
		//		}
		//		if (a1 >= 0.0)
		//		{
		//			tt1 = (a2 + 0.5 * sqrt(a1)) / a3;
		//			tt2 = (a2 - 0.5 * sqrt(a1)) / a3;
		//			//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< A << " " << B << " " << C << endl;
		//			cout << "Setka.cpp    " << "Fly_ex_Korol   "<< tt1 << " " << tt2 << endl;
		//			if (std::fpclassify(tt1) == FP_NORMAL && tt1 > 0.0)
		//			{
		//				if (t_per > tt1)
		//				{
		//					if (now->belong(x_0 + 1.03 * tt1 * Vx, sqrt(kv(y_0 + 1.03 * tt1 * Vy) + kv(z_0 + 1.03 * tt1 * Vz))) == false)
		//					{
		//						cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "A " << endl;
		//						t_per = tt1;
		//					}
		//				}
		//			}
		//			if (std::fpclassify(tt2) == FP_NORMAL && tt2 > 0.0)
		//			{
		//				if (t_per > tt2)
		//				{
		//					if (now->belong(x_0 + 1.03 * tt2 * Vx, sqrt(kv(y_0 + 1.03 * tt2 * Vy) + kv(z_0 + 1.03 * tt2 * Vz))) == false)
		//					{
		//						cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "B " << endl;
		//						t_per = tt2;
		//					}
		//				}
		//			}
		//		}
		//	}
		//	cout << "Setka.cpp    " << "Fly_ex_Korol   "<< t_per << endl;
		//	xk = x_0 + t_per * Vx;
		//	yk = sqrt(kv(y_0 + t_per * Vy) + kv(z_0 + t_per * Vz));
		//	cout << "Setka.cpp    " << "Fly_ex_Korol   "<< xk << " " << yk << endl;

		//	a = t_per * 0.99;
		//	b = t_per * 1.01;
		//}




		//exit(-1);
		for (auto& i : this->All_Cells)
		{
			if (i->belong(xk, yk) == true)
			{
				next = i;
				break;
			}
		}
	}

	if (next != nullptr)
	{
		if (next->zona > now->zona)
		{
			double Muuu = Mu_[area][next->zona][i_alp] * 0.3 * SINKR[i_alp] * kv(Ri[next->zona] / Rmax_);
			if (fabs(mu3) < Muuu)
			{
				if (fabs(mu3) >= MyRandom(s1, s2, s3) * Muuu)
				{
					mu3 = Muuu * sign(mu3);
				}
				else
				{
					return;
				}
			}
		}


		if (mu_statistic)
		{
			if (now->zona != zon_stat)
			{
				Mu_stat[area][now->zona] += mu;
				I_stat[area][now->zona]++;
				zon_stat = now->zona;
			}
		}

		Fly_exchenge_Imit_Korol_auto_weight(MK, s1, s2, s3, X, Y, Z, Vx, Vy, Vz, next, mu3, area, false, mu_start, to_I, to_J, georaschep, zon_stat);
	}
	else
	{

		if (xk > -1000.0 / RR_ && xk < 1000.0 / RR_ && yk < 1600.0 / RR_)
		{
			cout << "Setka.cpp    " << "Fly_ex_Korol   " << "ERROR  8670  poteryal atom" << endl;
			cout << "Setka.cpp    " << "Fly_ex_Korol   " << "SDFSDF    " << xk << " " << yk << endl;
		}
	}

	return;
}


void Setka::Fly_exchenge_Imit_Korol_2(MKmethod& MK, Sensor* sens, double x_0, double y_0, double z_0, double Vx, double Vy, double Vz, Cell* now, double mu, double KSI, //
	double I_do, int area, const double& mu_start)
{
	// Новая версия полёта нейтрального атома. Функция без лишнего кода и насыщена комментариями
	// Без геометрического и физического расщепления, полностью имитационный метод

	Cell* next_mb;
	double normV, dt, ddt, a, b, t_per, tt1, tt2, a1, a2, a3, A, B, C;
	double k, time, t1, t2, D;
	bool hand, slay, kj;
	double uz, uz_M, uz_E, cp, vx, vy, ro, t_ex;
	int lk, area2, num;

	double x_ex, y_ex, z_ex, u1, u2, u3, y_start, l, u, nu_ex, mu_ex, mu2, mu3, I, sig, ksi, k1, k2, k3, Vrr;
	double r, rr, alpha, aa, bb, cc, X, Y, Z, xk, yk, r_k, r_0;
	double Ur, Uthe, Uphi, Vr, Vthe, Vphi;
	Cell* next;

	vector <double> Wr(1);
	vector <double> Wthe(1);
	vector <double> Wphi(1);
	vector <double> mu_(1);


my_start:

	if (true)
	{
		// Блок для того, чтобы все переменные удалялись при выходе из блока
		next_mb = nullptr;

		normV = sqrt(kvv(Vx, Vy, Vz));
		dt = now->L / normV;  // Время на один мини-шаг в ячейке (планируется сделать много шагов в ячейке)
		ddt = 0.001 * dt;
		// Нахождение времени до попадания в следующую ячейку
		//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "V = " << Vx << " " << Vy << " " << Vz << endl;
		//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "A" << endl;
		a = 0.0;
		b = dt;
		t_per = 100000000000.0;
		// Блок проверки перигелия к оси
		/*double t_per = (-Vy * y_0 - Vz * z_0) / (kv(Vy) + kv(Vz));
		if (t_per > 0.0000001)
		{
			if (now->belong(x_0 + t_per * Vx, sqrt(kv(y_0 + t_per * Vy) + kv(z_0 + t_per * Vz))) == true)
			{
				t_per = -1.0;
			}
		}
		else
		{
			t_per = -1.0;
		}*/

		//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "A" << endl;
		
		//double tt1, tt2;
		//double a1, a2, a3;

		// Пытаемся определить время выхода точно
		if (true)
		{
			//double A, B, C;
			for (auto& i : now->Grans)
			{
				if (i->type == Axis)
				{
					continue;
				}
				A = i->aa;
				B = i->bb;
				C = i->cc;

				/*if (fabs(B) < 0.001)
				{
					tt1 = (-C / A - x_0) / Vx;
					if (std::fpclassify(tt1) == FP_NORMAL && tt1 > 0.0)
					{
						if (t_per > tt1)
						{
							if (now->belong(x_0 + (tt1 + min(ddt, 0.03 * tt1)) * Vx, sqrt(kv(y_0 + (tt1 + min(ddt, 0.03 * tt1)) * Vy) + kv(z_0 + (tt1 + min(ddt, 0.03 * tt1)) * Vz))) == false)
							{
								t_per = tt1;
								next_mb = i->Sosed;
							}
						}
					}
					continue;
				}*/
				/*cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "1 - " << A << " " << B << " " << C << endl;
				cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "2 - " << x_0 << " " << y_0 << " " << z_0 << endl;
				cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "3 - " << Vx << " " << Vy << " " << Vz << endl;*/
				a1 = 4.0 * kv(A * Vx * (C + A * x_0) - kv(B) * (Vy * y_0 + Vz * z_0)) - 4.0 * (kv(A * Vx) - //
					B * B * (kv(Vy) + kv(Vz))) * (kv(C + A * x_0) - B * B * (kv(y_0) + kv(z_0)));
				a2 = -A * Vx * (C + A * x_0) + B * B * (Vy * y_0 + Vz * z_0);
				a3 = A * A * Vx * Vx - B * B * (Vy * Vy + Vz * Vz);
				//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "4 - " << a1 << " " << a2 << " " << a3 << endl;
				if (fabs(a1) < 0.0000001)
				{
					a1 = 0.0;
				}
				if (a1 >= 0.0)
				{
					tt1 = (a2 + 0.5 * sqrt(a1)) / a3;
					tt2 = (a2 - 0.5 * sqrt(a1)) / a3;
					//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< A << " " << B << " " << C << endl;
					if (std::fpclassify(tt1) == FP_NORMAL && tt1 > 0.0)
					{
						if (t_per > tt1)
						{
							if (now->belong(x_0 + (tt1 + min(ddt, 0.005 * tt1)) * Vx, sqrt(kv(y_0 + (tt1 + min(ddt, 0.005 * tt1)) * Vy) + kv(z_0 + (tt1 + min(ddt, 0.01 * tt1)) * Vz))) == false)
							{
								t_per = tt1;
								next_mb = i->Sosed;
							}
						}
					}
					if (std::fpclassify(tt2) == FP_NORMAL && tt2 > 0.0)
					{
						if (t_per > tt2)
						{
							if (now->belong(x_0 + (tt2 + min(ddt, 0.005 * tt2)) * Vx, sqrt(kv(y_0 + (tt2 + min(ddt, 0.005 * tt2)) * Vy) + kv(z_0 + (tt2 + min(ddt, 0.01 * tt2)) * Vz))) == false)
							{
								t_per = tt2;
								next_mb = i->Sosed;
							}
						}
					}
				}
			}

			a = t_per - min(ddt, 0.005 * t_per);
			b = t_per + min(ddt, 0.005 * t_per);
		}

		hand = false;
		if (now->belong(x_0 + a * Vx, sqrt(kv(y_0 + a * Vy) + kv(z_0 + a * Vz))) == false)
		{
			hand = true;
		}
		else if (now->belong(x_0 + b * Vx, sqrt(kv(y_0 + b * Vy) + kv(z_0 + b * Vz))) == true)
		{
			hand = true;
		}
		else if (t_per > 100000000.0)
		{
			//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "ERROR  8568  " << t_per << endl;
			hand = true;
		}

		//double k;

		// Находим время  в ячейке
		lk = 0;
		if (hand)
		{
			a = 0.0;
			b = dt;
			t_per = 100000000000.0;
			for (auto& i : now->Grans)
			{
				a1 = -kv(i->a * Vx) * (kv(i->a * Vx) - kv(Vy) - kv(Vz)) * kv(Vz * y_0 - Vy * z_0);
				a2 = (kv(i->a * Vx) - kv(Vy) - kv(Vz)) * (Vz * z_0 + Vy * y_0);
				a3 = (kv(Vy) + kv(Vz)) * (-kv(i->a * Vx) + kv(Vy) + kv(Vz));
				if (a1 >= 0)
				{
					tt1 = (sqrt(a1) + a2) / a3;
					tt2 = (-sqrt(a1) + a2) / a3;
					if (std::fpclassify(tt1) == FP_NORMAL && tt1 > 0.0)
					{
						if (now->belong(x_0 + tt1 * Vx, sqrt(kv(y_0 + tt1 * Vy) + kv(z_0 + tt1 * Vz))) == false)
						{
							if (t_per > tt1)
							{
								t_per = tt1;
							}
						}
					}
					if (std::fpclassify(tt2) == FP_NORMAL && tt2 > 0.0)
					{
						if (now->belong(x_0 + tt2 * Vx, sqrt(kv(y_0 + tt2 * Vy) + kv(z_0 + tt2 * Vz))) == false)
						{
							if (t_per > tt2)
							{
								t_per = tt2;
							}
						}
					}
				}
			}

			while (now->belong(x_0 + b * Vx, sqrt(kv(y_0 + b * Vy) + kv(z_0 + b * Vz))) == true) // Слишком маленький шаг
			{
				lk++;
				if (lk > 1000)
				{
					cout << "Setka.cpp    " << "Fly_ex_Korol   " << "8055 dt = " << dt << " " << x_0 << " " << sqrt(kv(y_0) + kv(z_0)) << endl;
					cout << "Setka.cpp    " << "Fly_ex_Korol   " << normV << " " << Vx << " " << Vy << " " << Vz << endl;
					exit(-1);
				}
				b = b * 1.1;
			}

			if (t_per < b)
			{
				b = t_per;
			}

			k = (a + b) / 2.0;
			while (fabs(a - b) * normV > now->L / 300.0)
			{
				if (now->belong(x_0 + k * Vx, sqrt(kv(y_0 + k * Vy) + kv(z_0 + k * Vz))) == true)
				{
					a = k;
				}
				else
				{
					b = k;
				}
				k = (a + b) / 2.0;
			}
		}

		time = b; // Время нахождения в ячейке


		slay = false;
		// double t1, t2;
		//Нужно проверить пересечение со сферой
		D = kv(2.0 * Vx * x_0 + 2.0 * Vy * y_0 + 2.0 * Vz * z_0) - 4.0 * kvv(Vx, Vy, Vz) * (kvv(x_0, y_0, z_0) - kv(Rmax_));
		if (D >= 0.0)
		{
			t1 = (-2.0 * (Vx * x_0 + Vy * y_0 + Vz * z_0) - sqrt(D)) / (2.0 * kv(normV));
			t2 = (-2.0 * (Vx * x_0 + Vy * y_0 + Vz * z_0) + sqrt(D)) / (2.0 * kv(normV));
			if (t1 > 0.0 && t1 < time)
			{
				if (x_0 + Vx * t1 >= 0.0)
				{
					time = t1;
					a = time * 0.999;
					b = time * 1.001;
					slay = true;
				}
			}
			if (t2 > 0.0 && t2 < time)
			{
				if (x_0 + Vx * t2 >= 0.0)  // Для правой полусферы
				{
					time = t2;
					a = time * 0.999;
					b = time * 1.001;
					slay = true;
				}
			}
		}
		// Теперь пересечение с правой границей
		if (fabs(Vx) > 0.000001)
		{
			t1 = -x_0 / Vx;
			if (t1 > 0.0 && t1 <= time && kvv(0.0, y_0 + t1 * Vy, z_0 + t1 * Vz) >= kv(Rmax_))
			{
				time = t1;
				a = time * 0.999;
				b = time * 1.001;
				slay = true;
			}
		}


		// теперь время time, a, b уже определены
		// a - последний момент времени в ячейке, b - первый момент за пределами ячейки


		//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< a << " " << b << endl;
		//double uz, uz_M, uz_E;								// средние скорости в интегралах
		cp = sqrt(now->par[0].p / now->par[0].ro);
		vx = now->par[0].u;							// Скорости плазмы в ячейке
		vy = now->par[0].v;
		ro = now->par[0].ro;

		t_ex = 0.0;									// время до перезарядки
		// time - время нахождения атома в ячейке
		t2 = time;									// время после перезарядки (будет ниже)

		x_ex = x_0; 
		y_ex = y_0;
		z_ex = z_0;		    // координаты перезарядки
		//double u1, u2, u3;

		y_start = sqrt(kv(y_0) + kv(z_0));

		l = sqrt(kvv(time * Vx, time * Vy, time * Vz));                    // Расстояние до выхода из ячейки

		alpha = 0.0;								    // угол 
		u = 0.0;                                     // модуль относительной скорости атома и плазмы
		nu_ex = 0.0;								    // частота перезарядки
		//mu_ex = 0.0;									// вес перезаряженного атома
		mu2 = mu;									// вес не-перезаряженного атома
		mu3 = mu;
		mu_ex = mu;

		I = I_do;
		double drob = 3.0;
		//if (y_start < 40.0 / RR_)
		//{
		//	drob = 30.0;  // 30.0
		//}
		//else if (y_start < 60.0 / RR_)
		//{
		//	drob = 20.0;  // 20.0
		//}
		//else if (y_start < 80.0 / RR_)
		//{
		//	drob = 15.0;
		//}

		double vy_ = vy;
		nu_ex = 0.0;
		for (double tt = 0.0; tt < time * 0.99; tt += time / drob)
		{
			vy_ = vy;
			double yy = sqrt(kv(y_0 + (tt + time / (2.0 * drob)) * Vy) + kv(z_0 + (tt + time / (2.0 * drob)) * Vz));

			alpha = polar_angle(y_0 + (tt + time / (2.0 * drob)) * Vy, z_0 + (tt + time / (2.0 * drob)) * Vz);
			u = sqrt(kvv(Vx - vx, Vy - vy_ * cos(alpha), Vz - vy_ * sin(alpha)));

			if (true)//(u / cp > 7.0)
			{
				uz = Velosity_1(u, cp);
				nu_ex += (ro * uz * sigma(uz)) / Kn_;
			}
			else
			{
				nu_ex += (ro * MK.int_1(u, cp)) / Kn_;        // Пробуем вычислять интеграллы численно
			}
			//kappa += (nu_ex * time / drob);
		}

		nu_ex = nu_ex / drob;

		//alpha = polar_angle(y_0 + time / 2.0 * Vy, z_0 + time / 2.0 * Vz);
		//u = sqrt(kvv(Vx - vx, Vy - vy * cos(alpha), Vz - vy * sin(alpha)));
		//if (u / cp > 7.0)
		//{
		//	uz = Velosity_1(u, cp);
		//	nu_ex = (ro * uz * sigma(uz)) / Kn_;
		//	//cout << x_0 << " " << y_start << " " << u / cp << " " << u << " " << cp <<  endl;
		//}
		//else
		//{
		//	//uz = Velosity_1(u, cp);
		//	nu_ex = (ro * MK.int_1(u, cp)) / Kn_;
		//	//double sig = sqrt(kvv(Vx, Vy, Vz)) / ((1.0 / Kn_) * ro * uz * sigma(uz));
		//}

		sig = sqrt(kvv(Vx, Vy, Vz)) / (nu_ex);
		//sig = sqrt(kvv(Vx, Vy, Vz)) / (kappa);
		if (charge_x_1)
		{
			if (x_0 <= 1.0)
			{
				I += l / sig;
			}
		}
		else
		{
			I += l / sig;
		}

		/*if (now->type == C_2 || now->type == C_3 || now->type == C_1)
		{
			cout << "Ploho v c2, c3, c1 " << endl;
		}*/

		if (I < KSI)  // Не произошла перезарядка
		{
			I_do = I;
		}
		else  // Перезарядка была
		{
			ksi = (KSI - I_do) * sig;
			t_ex = ksi / sqrt(kvv(Vx, Vy, Vz));
			I_do = 0.0;
			alpha = polar_angle(y_0 + 0.5 * t_ex * Vy, z_0 + 0.5 * t_ex * Vz);  // Гарантирует расчёт угла в середине пути

			x_ex = x_0 + t_ex * Vx;  // Координаты перезарядки
			y_ex = y_0 + t_ex * Vy;
			z_ex = z_0 + t_ex * Vz;
			if (now->belong(x_ex, sqrt(kv(y_ex) + kv(z_ex))) == false)
			{
				x_ex = x_0 + a * Vx;  // Координаты перезарядки
				y_ex = y_0 + a * Vy;
				z_ex = z_0 + a * Vz;
				t_ex = a;
				alpha = polar_angle(y_0 + 0.5 * (t_ex)*Vy, z_0 + 0.5 * (t_ex)*Vz);
			}


			// Считаем источники для частицы до перезарядки
			u = sqrt(kvv(Vx - vx, Vy - vy * cos(alpha), Vz - vy * sin(alpha)));
			u1 = vx - Vx;
			u2 = vy * cos(alpha) - Vy;
			u3 = vy * sin(alpha) - Vz;
			double skalar = Vx * u1 + Vy * u2 + Vz * u3;

			now->mut.lock();

			//now->par[0].w_m[area] += mu / max(sin(alpha), 0.3 * Sinus[i_alp]);
			now->par[0].F_n += t_ex * mu;
			now->par[0].F_u += t_ex * Vx * mu;
			now->par[0].F_T += t_ex * kvv(Vx, Vy, Vz) * mu;
			now->par[0].H_n[area] += t_ex * mu;
			now->par[0].H_u[area] += t_ex * Vx * mu;
			now->par[0].H_T[area] += t_ex * kvv(Vx, Vy, Vz) * mu;

			area2 = 0;
			// определим область в которой находится атом сейчас (это параметр самой ячейки)
			if (now->type == C_5)
			{
				area2 = 3;
			}
			else if (now->type == C_4)
			{
				area2 = 2;
			}
			else if (now->type == C_1 || now->type == C_centr)
			{
				area2 = 0;
			}
			else
			{
				area2 = 1;
			}


			if (u / cp > 7.0)
			{
				uz = Velosity_1(u, cp);
				uz_M = Velosity_2(u, cp) / (uz * kv(cp) * cp * pi * sqrtpi_);
				uz_E = Velosity_3(u, cp);

				now->par[0].I_u += -mu_ex * uz_M * u1 / u;
				now->par[0].I_v += -mu_ex * uz_M * (u2 * cos(alpha) + u3 * sin(alpha)) / u;
				now->par[0].I_T += mu_ex * (-0.25 * (3.0 * kv(cp) + 2.0 * kv(u)) * (uz_E / uz) - //
					uz_M * skalar / u);
			}
			else
			{
				k1 = MK.int_1(u, cp);
				k2 = MK.int_2(u, cp);
				k3 = MK.int_3(u, cp);
				now->par[0].I_u += mu_ex * (k2 / k1) * u1 / u;
				now->par[0].I_v += mu_ex * (k2 / k1) * (u2 * cos(alpha) + u3 * sin(alpha)) / u;
				now->par[0].I_T += mu_ex * (-0.5 * k3 / k1 + k2 / k1 * skalar / u);
			}

			now->par[0].F_v += t_ex * (Vy * cos(alpha) + Vz * sin(alpha)) * mu;
			now->par[0].H_v[area] += t_ex * (Vy * cos(alpha) + Vz * sin(alpha)) * mu;
			Vrr = (Vy * cos(alpha) + Vz * sin(alpha));
			now->par[0].H_uu[area] += t_ex * Vx * Vx * mu;
			now->par[0].H_uv[area] += t_ex * Vrr * Vx * mu;
			now->par[0].H_vv[area] += t_ex * Vrr * Vrr * mu;
			now->par[0].H_uuu[area] += t_ex * Vx * Vx * Vx * mu;

			now->mut.unlock();


			r = sqrt(kvv(x_ex, y_ex, z_ex));
			rr = sqrt(kvv(0.0, y_ex, z_ex));
			//double Ur, Uthe, Uphi, Vr, Vthe, Vphi;

			/*if (r < 1.0)
			{
				cout << "Zaletel   " << mu << " " << x_ex << " " << rr << endl;
			}*/


			spherical_skorost(y_ex, z_ex, x_ex, vy * cos(alpha), vy * sin(alpha), vx, Ur, Uphi, Uthe);
			spherical_skorost(y_ex, z_ex, x_ex, Vy, Vz, Vx, Vr, Vphi, Vthe);

			if (r > Rmax_ && x_ex > 0.001 / RR_)
			{
				return;
			}



			if (true)    // Без расщепления
			{
				num = 1;
				kj = true;
				mu3 = mu;

				alpha = polar_angle(y_ex, z_ex);
				//double aa, bb, cc;
				//vector <double> Wr(1);
				//vector <double> Wthe(1);
				//vector <double> Wphi(1);
				//vector <double> mu_(1);

				MK.Change_Velosity4(sens, Ur / cp, Uthe / cp, Uphi / cp, Vr / cp, Vthe / cp, Vphi / cp, Wr, Wthe, Wphi, mu_, cp, r, 0);
				Wr[0] = Wr[0] * cp;
				Wthe[0] = Wthe[0] * cp;
				Wphi[0] = Wphi[0] * cp;

				dekard_skorost(y_ex, z_ex, x_ex, Wr[0], Wphi[0], Wthe[0], bb, cc, aa);


				if (true)//(Wr[0] > -0.4 && (now->type != C_3 && now->type != C_2 && now->type != C_1))
				{
					KSI = -log(1.0 - sens->MakeRandom());

					area = area2;
					x_0 = x_ex;
					y_0 = y_ex;
					z_0 = z_ex;
					Vx = aa;
					Vy = bb;
					Vz = cc;
					mu = mu3;

					goto my_start;

					 //Fly_exchenge_Imit_Korol_2(MK, sens, x_ex, y_ex, z_ex, aa, bb, cc, now, mu3, KSI, I_do, area2, mu_start);


			
				}
				else
				{
					Fly_exchenge_Imit_Korol_PUI(MK, sens, x_ex, y_ex, z_ex, aa, bb, cc, now, mu3, area2, false, mu_start, I_ - 1, J_ - 1, true);
				}
				return;
			}

		}


		//cout << "Setka.cpp    " << "B1" << endl;
		alpha = polar_angle(y_ex + 0.5 * t2 * Vy, z_ex + 0.5 * t2 * Vz);
		//cout << "Setka.cpp    " << "B2" << endl;
		now->mut.lock();
		now->par[0].F_n += t2 * mu2;
		now->par[0].F_u += t2 * Vx * mu2;
		now->par[0].F_v += t2 * (Vy * cos(alpha) + Vz * sin(alpha)) * mu2;
		now->par[0].F_T += t2 * kvv(Vx, Vy, Vz) * mu2;
		now->par[0].H_n[area] += t2 * mu2;
		now->par[0].H_u[area] += t2 * Vx * mu2;
		now->par[0].H_v[area] += t2 * (Vy * cos(alpha) + Vz * sin(alpha)) * mu2;
		now->par[0].H_T[area] += t2 * kvv(Vx, Vy, Vz) * mu2;

		Vrr = (Vy * cos(alpha) + Vz * sin(alpha));
		now->par[0].H_uu[area] += t2 * Vx * Vx * mu2;
		now->par[0].H_uv[area] += t2 * Vrr * Vx * mu2;
		now->par[0].H_vv[area] += t2 * Vrr * Vrr * mu2;
		now->par[0].H_uuu[area] += t2 * Vx * Vx * Vx * mu2;
		now->mut.unlock();

		X = x_ex + t2 * Vx;
		Y = y_ex + t2 * Vy;
		Z = z_ex + t2 * Vz;


		next = nullptr;             // В какую ячейку попадаем следующую
		//double xk, yk;
		xk = x_0 + b * Vx;
		yk = sqrt(kv(y_0 + b * Vy) + kv(z_0 + b * Vz));
		r_k = sqrt(kv(xk) + kv(yk));
		r_0 = sqrt(kv(x_0) + kv(y_0) + kv(z_0));

		//f_way << xk << " " << yk << " " << to_I << endl;
		//f_num++;

		if (r_0 <= Rmax_ && r_k >= Rmax_ && x_0 > 0.0)
		{
			//cout << "Setka.cpp    " << x_0 << " " << sqrt(kvv(0.0, y_0, z_0)) << " " << xk << " " << yk << endl;
			next = nullptr;
			return;
		}
		else if (r_k >= Rmax_ - 0.001 / RR_ && xk >= -0.001 / RR_)
		{
			//cout << "Setka.cpp    " << x_0 << " " << sqrt(kvv(0.0, y_0, z_0)) << " " << xk << " " << yk << endl;
			next = nullptr;
			return;
		}
		else
		{
			for (auto& i : now->Grans)
			{
				if (i->Sosed != nullptr)
				{
					if (i->Sosed->belong(xk, yk) == true)
					{
						next = i->Sosed;
					}
				}
			}
			if (next == nullptr)
			{
				for (auto& k : now->Grans)
				{
					if (k->Sosed != nullptr)
					{
						for (auto& i : k->Sosed->Grans)
						{
							if (i->Sosed != nullptr)
							{
								if (i->Sosed->belong(xk, yk) == true)
								{
									next = i->Sosed;
								}
							}
						}
					}
				}
			}
			if (next == nullptr)
			{
				if (now->belong(xk, yk) == true)
				{
					cout << "Ostalsya v svoey " << endl;
					next = now;
				}
			}
		}

		if (next == nullptr)
		{
			for (auto& i : this->All_Cells_zero)
			{
				if (i->belong(xk, yk) == true)
				{
					next = i;
				}
			}
		}

		if (next != nullptr)
		{
			x_0 = X;
			y_0 = Y;
			z_0 = Z;
			now = next;
			mu = mu2;

			goto my_start;
			//Fly_exchenge_Imit_Korol_2(MK, sens, X, Y, Z, Vx, Vy, Vz, next, mu2, KSI, I_do, area, mu_start);
		}
		else
		{
			if (xk > -1440.0 / RR_ && xk < 1200.0 / RR_ && yk < 200.0 / RR_) // Проверка, если траектория пропадёт не на границе а в нуле скорее всего
			{
				cout << "Setka.cpp  15246  propal  " << xk << " " << yk << " " << Vx << " " << Vy << " " << Vz << " " << time << endl;
				cout << x_0 << " " << y_0 << endl;
				double cv1, cv2;
				now->Get_Center(cv1, cv2);
				cout << cv1 << " " << cv2 << endl;
				for (auto& i : this->All_Cells)
				{
					if (i->belong(xk, yk) == true)
					{
						cout << "Nashol " << endl;
						next = i;

						x_0 = X;
						y_0 = Y;
						z_0 = Z;
						now = next;
						mu = mu2;

						goto my_start;

						//Fly_exchenge_Imit_Korol_2(MK, sens, X, Y, Z, Vx, Vy, Vz, next, mu2, KSI, I_do, area, mu_start);

					}
				}
			}
		}

		//cout << "Setka.cpp    " << "Vihod" << endl;
		return;
	}
}

void Setka::Fly_exchenge_Imit_Korol_PUI(MKmethod& MK, Sensor* sens, double x_0, double y_0, double z_0, double Vx, double Vy, double Vz, Cell* now, double mu, //
	int area, bool ExCh, const double& mu_start, int to_I = 0, int to_J = 0, bool georaschep = true, int zon_stat, bool pui__)
	// Новая версия полёта нейтрального атома. Функция без лишнего кода и насыщена комментариями
	// Функция использует и геометрическое и физическое расщепление!!!
	// Описание переменных:
	// sens - датчик случайных чисел
	// x_0 ... -  координаты начала полёта нейтрального атома в данной ячейке
	// Vx ... -   скорость летящего нейтрального атома
	// now -      текущая ячейка
	// mu -       вес нейтрального атома
	// area -     какой атом летит
	// mu_start - начальный вес атома при инициализации (чтобы понимать, сколько от исходного осталось)
	// to_I - в какую область предназначался атом
{
	Cell* next_mb = nullptr;

	// Алгоритм вырубания лишних траектори (траектории вырубаются при пересечении новой зоны, далее 
	// вес лмбо увеличивается, либо частица уничтожается)
	// Не уверен, что to_I бывает больше, чем now->zona, кажется не должно такого быть


	double normV = sqrt(kvv(Vx, Vy, Vz));
	double dt = now->L / normV;  // Время на один мини-шаг в ячейке (планируется сделать много шагов в ячейке)
	double ddt = 0.001 * dt;
	// Нахождение времени до попадания в следующую ячейку

	double a = 0.0, b = dt;
	double t_per = 100000000000.0;

	//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "A" << endl;
	double tt1, tt2;
	double a1, a2, a3;

	// Пытаемся определить время выхода точно
	if (true)
	{
		double A, B, C;
		for (auto& i : now->Grans)
		{
			if (i->type == Axis)
			{
				continue;
			}
			A = i->aa;
			B = i->bb;
			C = i->cc;

			a1 = 4.0 * kv(A * Vx * (C + A * x_0) - kv(B) * (Vy * y_0 + Vz * z_0)) - 4.0 * (kv(A * Vx) - //
				B * B * (kv(Vy) + kv(Vz))) * (kv(C + A * x_0) - B * B * (kv(y_0) + kv(z_0)));
			a2 = -A * Vx * (C + A * x_0) + B * B * (Vy * y_0 + Vz * z_0);
			a3 = A * A * Vx * Vx - B * B * (Vy * Vy + Vz * Vz);

			if (fabs(a1) < 0.0000001)
			{
				a1 = 0.0;
			}
			if (a1 >= 0.0)
			{
				tt1 = (a2 + 0.5 * sqrt(a1)) / a3;
				tt2 = (a2 - 0.5 * sqrt(a1)) / a3;
				//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< A << " " << B << " " << C << endl;
				if (std::fpclassify(tt1) == FP_NORMAL && tt1 > 0.0)
				{
					if (t_per > tt1)
					{
						if (now->belong(x_0 + (tt1 + min(ddt, 0.005 * tt1)) * Vx, sqrt(kv(y_0 + (tt1 + min(ddt, 0.005 * tt1)) * Vy) + kv(z_0 + (tt1 + min(ddt, 0.01 * tt1)) * Vz))) == false)
						{
							t_per = tt1;
							next_mb = i->Sosed;
						}
					}
				}
				if (std::fpclassify(tt2) == FP_NORMAL && tt2 > 0.0)
				{
					if (t_per > tt2)
					{
						if (now->belong(x_0 + (tt2 + min(ddt, 0.005 * tt2)) * Vx, sqrt(kv(y_0 + (tt2 + min(ddt, 0.005 * tt2)) * Vy) + kv(z_0 + (tt2 + min(ddt, 0.01 * tt2)) * Vz))) == false)
						{
							t_per = tt2;
							next_mb = i->Sosed;
						}
					}
				}
			}
		}

		a = t_per - min(ddt, 0.005 * t_per);
		b = t_per + min(ddt, 0.005 * t_per);
	}

	bool hand = false;
	if (now->belong(x_0 + a * Vx, sqrt(kv(y_0 + a * Vy) + kv(z_0 + a * Vz))) == false)
	{
		hand = true;
	}
	else if (now->belong(x_0 + b * Vx, sqrt(kv(y_0 + b * Vy) + kv(z_0 + b * Vz))) == true)
	{
		hand = true;
	}
	else if (t_per > 100000000.0)
	{
		//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "ERROR  8568  " << t_per << endl;
		hand = true;
	}

	double k;

	// Находим время  в ячейке
	int lk = 0;
	if (hand)
	{
		a = 0.0;
		b = dt;
		t_per = 100000000000.0;
		for (auto& i : now->Grans)
		{
			a1 = -kv(i->a * Vx) * (kv(i->a * Vx) - kv(Vy) - kv(Vz)) * kv(Vz * y_0 - Vy * z_0);
			a2 = (kv(i->a * Vx) - kv(Vy) - kv(Vz)) * (Vz * z_0 + Vy * y_0);
			a3 = (kv(Vy) + kv(Vz)) * (-kv(i->a * Vx) + kv(Vy) + kv(Vz));
			if (a1 >= 0)
			{
				tt1 = (sqrt(a1) + a2) / a3;
				tt2 = (-sqrt(a1) + a2) / a3;
				if (std::fpclassify(tt1) == FP_NORMAL && tt1 > 0.0)
				{
					if (now->belong(x_0 + tt1 * Vx, sqrt(kv(y_0 + tt1 * Vy) + kv(z_0 + tt1 * Vz))) == false)
					{
						if (t_per > tt1)
						{
							t_per = tt1;
						}
					}
				}
				if (std::fpclassify(tt2) == FP_NORMAL && tt2 > 0.0)
				{
					if (now->belong(x_0 + tt2 * Vx, sqrt(kv(y_0 + tt2 * Vy) + kv(z_0 + tt2 * Vz))) == false)
					{
						if (t_per > tt2)
						{
							t_per = tt2;
						}
					}
				}
			}
		}

		while (now->belong(x_0 + b * Vx, sqrt(kv(y_0 + b * Vy) + kv(z_0 + b * Vz))) == true) // Слишком маленький шаг
		{
			lk++;
			if (lk > 1000)
			{
				cout << "Setka.cpp    " << "Fly_ex_Korol   " << "8055 dt = " << dt << " " << x_0 << " " << sqrt(kv(y_0) + kv(z_0)) << endl;
				cout << "Setka.cpp    " << "Fly_ex_Korol   " << normV << " " << Vx << " " << Vy << " " << Vz << endl;
				exit(-1);
			}
			b = b * 1.1;
		}

		if (t_per < b)
		{
			b = t_per;
		}

		k = (a + b) / 2.0;
		while (fabs(a - b) * normV > now->L / 300.0)
		{
			if (now->belong(x_0 + k * Vx, sqrt(kv(y_0 + k * Vy) + kv(z_0 + k * Vz))) == true)
			{
				a = k;
			}
			else
			{
				b = k;
			}
			k = (a + b) / 2.0;
		}
	}

	double time = b; // Время нахождения в ячейке


	bool slay = false;
	double t1, t2;
	//Нужно проверить пересечение со сферой
	double D = kv(2.0 * Vx * x_0 + 2.0 * Vy * y_0 + 2.0 * Vz * z_0) - 4.0 * kvv(Vx, Vy, Vz) * (kvv(x_0, y_0, z_0) - kv(Rmax_));
	if (D >= 0.0)
	{
		t1 = (-2.0 * (Vx * x_0 + Vy * y_0 + Vz * z_0) - sqrt(D)) / (2.0 * kv(normV));
		t2 = (-2.0 * (Vx * x_0 + Vy * y_0 + Vz * z_0) + sqrt(D)) / (2.0 * kv(normV));
		if (t1 > 0.0 && t1 < time)
		{
			if (x_0 + Vx * t1 >= 0.0)
			{
				time = t1;
				a = time * 0.999;
				b = time * 1.001;
				slay = true;
			}
		}
		if (t2 > 0.0 && t2 < time)
		{
			if (x_0 + Vx * t2 >= 0.0)  // Для правой полусферы
			{
				time = t2;
				a = time * 0.999;
				b = time * 1.001;
				slay = true;
			}
		}
	}
	// Теперь пересечение с правой границей
	if (fabs(Vx) > 0.000001)
	{
		t1 = -x_0 / Vx;
		if (t1 > 0.0 && t1 <= time && kvv(0.0, y_0 + t1 * Vy, z_0 + t1 * Vz) >= kv(Rmax_))
		{
			time = t1;
			a = time * 0.999;
			b = time * 1.001;
			slay = true;
		}
	}


	// теперь время time, a, b уже определены
	// a - последний момент времени в ячейке, b - первый момент за пределами ячейки


	//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< a << " " << b << endl;
	double uz, uz_M, uz_E;								// средние скорости в интегралах
	double cp = sqrt(now->par[0].p / now->par[0].ro);

	if (std::fpclassify(cp) != FP_NORMAL)
	{
		cout << cp << "    ERROR  cp ihwuyegfyguy234223424" << endl;
		cout << now->par[0].p << " " << now->par[0].ro << endl;
		cout << now->contour[0]->x << " " << now->contour[0]->y << endl;
		exit(-1);
	}

	double vx = now->par[0].u;							// Скорости плазмы в ячейке
	double vy = now->par[0].v;
	double ro = now->par[0].ro;


	double t_ex = 0.0;									// время до перезарядки
	// time - время нахождения атома в ячейке
	t2 = time;									// время после перезарядки (будет ниже)

	double x_ex = x_0, y_ex = y_0, z_ex = z_0;		    // координаты перезарядки
	double u1, u2, u3;

	double y_start = sqrt(kv(y_0) + kv(z_0));

	double l = sqrt(kvv(time * Vx, time * Vy, time * Vz));                    // Расстояние до выхода из ячейки

	double alpha = 0.0;								    // угол 
	double u = 0.0;                                     // модуль относительной скорости атома и плазмы
	double nu_ex = 0.0;								    // частота перезарядки
	double nu_ex_pui = 0.0;
	double mu_ex = 0.0;									// вес перезаряженного атома
	double mu_ex_p = 0.0;
	vector <double> mu_ex_pui;
	mu_ex_pui.resize(now->i_pui);
	double mu2 = mu;									// вес не-перезаряженного атома

	double drob = 5.0;     // Для того чтобы точнее учесть проворот скорости во время полёта

	double kappa_sum = 0.0;

	if (now->type != C_centr)//(ExCh == false) // Если перезарядки в этой ячейке ещё не было
	{
		double kappa = 0.0;
		vector <double> kappa_pui;
		kappa_pui.resize(now->i_pui);
		for (int i = 0; i < now->i_pui; i++)
		{
			kappa_pui[i] = 0.0;
		}
		//drob = min(fabs(polar_angle(y_0, z_0) - polar_angle(y_0 + time * Vy, z_0 + time * Vz)) / 0.017 + 5.0, 80.0);

		//if (y_start < 40.0/RR_)
		//{
		//	drob = 30.0;  // 30.0
		//}
		//else if (y_start < 60.0 / RR_)
		//{
		//	drob = 20.0;  // 20.0
		//}
		//else if (y_start < 80.0 / RR_)
		//{
		//	drob = 15.0;
		//}
		//else if (y_start < 120.0 / RR_)
		//{
		//	drob = 10.0;
		//}


		if (true)  // Если нужно дробить траекторию
		{
			for (double tt = 0.0; tt < time * 0.98; tt += time / drob)
			{
				double yy = sqrt(kv(y_0 + (tt + time / (2.0 * drob)) * Vy) + kv(z_0 + (tt + time / (2.0 * drob)) * Vz));

				alpha = polar_angle(y_0 + (tt + time / (2.0 * drob)) * Vy, z_0 + (tt + time / (2.0 * drob)) * Vz);
				u = sqrt(kvv(Vx - vx, Vy - vy * cos(alpha), Vz - vy * sin(alpha)));

				if (std::fpclassify(u) != FP_NORMAL && std::fpclassify(u) != FP_ZERO)
				{
					cout << u << "    ERROR  uuu weferwfw43423ewefdef" << endl;
					cout << Vx << " " << Vy << " " << Vz << " " << vx << " " << vy << " " << alpha << endl;
					exit(-1);
				}

				if (now->pui_ == true)
				{
					for (int i = 0; i < now->i_pui; i++)
					{
						nu_ex_pui = now->get_nu_pui(u, i) / Kn_;   // Считаем частоту для перезаряженного на пикапах атома
						kappa_pui[i] += (nu_ex_pui * time / drob);
					}
				}

				if (u / cp > 7.0)
				{
					uz = Velosity_1(u, cp);
					nu_ex = (ro * uz * sigma(uz)) / Kn_;
				}
				else
				{
					nu_ex = (ro * MK.int_1(u, cp)) / Kn_;        // Пробуем вычислять интеграллы численно
				}
				kappa += (nu_ex * time / drob);
			}
		}

		kappa_sum = kappa;
		for (int i = 0; i < now->i_pui; i++)
		{
			kappa_sum = kappa_sum + kappa_pui[i];
		}



		t_ex = -(time / kappa_sum) * log(1.0 - sens->MakeRandom() * (1.0 - exp(-kappa_sum))); // Время до перезарядки
		t2 = time - t_ex;   // Время сколько лететь после того, как атом перезарядился
		mu_ex = mu * (1.0 - exp(-kappa_sum)); // вес перезаряженного атома
		mu2 = mu - mu_ex;                 // вес оставшегося неперезаряженного атома
		mu_ex_p = (kappa / kappa_sum) * mu_ex;
		for (int i = 0; i < now->i_pui; i++)
		{
			mu_ex_pui[i] = (kappa_pui[i] / kappa_sum) * mu_ex;
		}

		x_ex = x_0 + t_ex * Vx;  // Координаты перезарядки
		y_ex = y_0 + t_ex * Vy;
		z_ex = z_0 + t_ex * Vz;

		alpha = polar_angle(y_ex, z_ex);

		// проверка, если перезарядка произошла за пределами ячейки (нужно немного сдвигать её в этом случае)
		if (now->belong(x_ex, sqrt(kv(y_ex) + kv(z_ex))) == false)
		{
			x_ex = x_0 + a * Vx;  // Координаты перезарядки
			y_ex = y_0 + a * Vy;
			z_ex = z_0 + a * Vz;
			t_ex = a;
			t2 = time - t_ex;
			alpha = polar_angle(y_ex, z_ex);
		}

		double r = sqrt(kvv(x_ex, y_ex, z_ex));
		double all = polar_angle(x_ex, sqrt(kv(y_ex) + kv(z_ex)));

		// для сбора статистики
		if (func_stat)
		{
			if (r < 100.0 / RR_ && r > 60.0 / RR_)
			{
				bool sec = false;
				bool sec2 = false;
				double tt;
				double a1 = Vx * x_0 + Vy * y_0 + Vz * z_0;
				double a2 = kvv(Vx, Vy, Vz);
				double a3 = kvv(x_0, y_0, z_0) - kv(R_stat / RR_);
				if (4.0 * a1 * a1 - 4.0 * a2 * a3 >= 0.0)
				{
					double t1 = (-2.0 * a1 + sqrt(4.0 * a1 * a1 - 4.0 * a2 * a3)) / (2.0 * a2);
					double t2 = (-2.0 * a1 - sqrt(4.0 * a1 * a1 - 4.0 * a2 * a3)) / (2.0 * a2);
					if (t1 > 0.000000001 && t1 < time)
					{
						sec = true;
						tt = t1;
						if (t1 < t_ex)
						{
							sec2 = true; // пролетел до перезярядки
						}
						else
						{
							sec2 = false;
						}
					}
					if (t2 > 0.000000001 && t2 < time)
					{
						sec = true;
						tt = t2;
						if (t2 < t_ex)
						{
							sec2 = true; // пролетел до перезярядки
						}
						else
						{
							sec2 = false;
						}
					}
				}

				if (sec == true)
				{
					double xx, yy, zz;
					xx = x_0 + tt * Vx;
					yy = y_0 + tt * Vy;
					zz = z_0 + tt * Vz;
					double Vr, Vphi, Vthe;
					double phi_ = polar_angle(yy, zz);  //polar_angle(xx, sqrt(kv(yy) + kv(zz)));
					spherical_skorost(yy, zz, xx, Vy, Vz, Vx, Vr, Vphi, Vthe);

					double Vn = fabs((Vx * xx + Vy * yy + Vz * zz) / (R_stat / RR_));
					//double Vn = sqrt(kvv(Vx, Vy, Vz));
					double Vrr = (Vy * cos(phi_) + Vz * sin(phi_));
					double mumu = 0.0;
					if (sec2)
					{
						mumu = mu / Vn;
					}
					else
					{
						mumu = mu2 / Vn;
					}

					mut_mom.lock();
					this->mu_mom[area][(int)(phi_ / (pi_ / Al_stat))] += mumu;
					this->Vx_mom[area][(int)(phi_ / (pi_ / Al_stat))] += mumu * Vx;
					this->Vy_mom[area][(int)(phi_ / (pi_ / Al_stat))] += mumu * Vrr;
					this->Vxx_mom[area][(int)(phi_ / (pi_ / Al_stat))] += mumu * Vx * Vx;
					this->Vyy_mom[area][(int)(phi_ / (pi_ / Al_stat))] += mumu * Vrr * Vrr;
					this->Vxy_mom[area][(int)(phi_ / (pi_ / Al_stat))] += mumu * Vx * Vrr;
					this->Vxxx_mom[area][(int)(phi_ / (pi_ / Al_stat))] += mumu * Vx * Vx * Vx;
					mut_mom.unlock();

					//if (fabs(mu) > pred(area, now->zona, all) * 0.3)
					//{
					//	if (phi_ < 0.17453293)  // Угол в 10 градусов
					//	{
					//		mut_stat.lock();
					//		this->V_r_stat[this->number_stat] = Vr;
					//		this->V_t_stat[this->number_stat] = Vthe;
					//		this->V_p_stat[this->number_stat] = Vphi;
					//		this->phi_stat[this->number_stat] = phi_;
					//		this->num_stat[this->number_stat] = area;
					//		if (sec2)
					//		{
					//			this->mu_stat[this->number_stat] = mu / Vn;
					//		}
					//		else
					//		{
					//			this->mu_stat[this->number_stat] = mu2 / Vn;
					//		}
					//		this->number_stat++;
					//		mut_stat.unlock();
					//	}
					//}
				}
			}
		}

		double rr = sqrt(kvv(0.0, y_ex, z_ex));


		// Считаем источники для частицы до перезарядки

		u = sqrt(kvv(Vx - vx, Vy - vy * cos(alpha), Vz - vy * sin(alpha)));
		u1 = vx - Vx;
		u2 = vy * cos(alpha) - Vy;
		u3 = vy * sin(alpha) - Vz;
		double skalar = Vx * u1 + Vy * u2 + Vz * u3;

		double skalar2 = -Vx * u1 / u + -Vy * u2 / u + -Vz * u3 / u;

		// определим область перезарядки (по геометрическим зонам)
		int i_z = 0, i_alp = 0;
		i_z = geo_zones(r);
		i_alp = alpha_zones(x_ex, rr);

		int area2 = 0;
		// определим область в которой находится атом сейчас (это параметр самой ячейки)
		if (now->type == C_5)
		{
			area2 = 3;
		}
		else if (now->type == C_4)
		{
			area2 = 2;
		}
		else if (now->type == C_1 || now->type == C_centr)
		{
			area2 = 0;
		}
		else
		{
			area2 = 1;
		}


		now->mut.lock();

		//now->par[0].w_m[area] += mu / max(sin(alpha), 0.3 * Sinus[i_alp]);
		now->par[0].F_n += t_ex * mu;
		now->par[0].F_u += t_ex * Vx * mu;
		now->par[0].F_T += t_ex * kvv(Vx, Vy, Vz) * mu;

		int ghdf = min((int)(u * n_S / max_S), n_S - 1);
		now->S_p[ghdf] += t_ex * mu * (kappa_sum/time);
		now->S_m[ghdf] += t_ex * mu;

		int dfg = area;
		if (pui__ == true)
		{
			if (area == 0)
			{
				dfg = 4;
			}
			if (area == 1)
			{
				dfg = 5;
			}
		}

		now->par[0].H_n[dfg] += t_ex * mu;
		now->par[0].H_u[dfg] += t_ex * Vx * mu;
		now->par[0].H_T[dfg] += t_ex * kvv(Vx, Vy, Vz) * mu;

		double Vrr = (Vy * cos(alpha) + Vz * sin(alpha));
		now->par[0].H_uu[dfg] += t_ex * Vx * Vx * mu;
		now->par[0].H_uv[dfg] += t_ex * Vrr * Vx * mu;
		now->par[0].H_vv[dfg] += t_ex * Vrr * Vrr * mu;
		now->par[0].H_uuu[dfg] += t_ex * Vx * Vx * Vx * mu;

		if (now->pui_ == true)
		{
			for (int i = 0; i < now->i_pui; i++)
			{
				double k1 = now->get_nu_pui(u, i);
				double k2 = now->get_nu_pui2(u, i);
				double k3 = now->get_nu_pui3(u, i);
				now->par[0].II_u += -mu_ex_pui[i] * (k2 / k1) * u1 / u;
				now->par[0].II_v += -mu_ex_pui[i] * (k2 / k1) * (u2 * cos(alpha) + u3 * sin(alpha)) / u;
				now->par[0].II_T += mu_ex_pui[i] * (-0.5 * k3 / k1 + k2 / k1 * skalar2);
			}
		}

		if (u / cp > 7.0)
		{
			uz = Velosity_1(u, cp);
			uz_M = Velosity_2(u, cp) / (uz * kv(cp) * cp * pi * sqrtpi_);
			uz_E = Velosity_3(u, cp);

			now->par[0].I_u += -mu_ex_p * uz_M * u1 / u;
			now->par[0].I_v += -mu_ex_p * uz_M * (u2 * cos(alpha) + u3 * sin(alpha)) / u;
			now->par[0].I_T += mu_ex_p * (-0.25 * (3.0 * kv(cp) + 2.0 * kv(u)) * (uz_E / uz) - //
				uz_M * skalar / u);


		}
		else
		{
			double k1 = MK.int_1(u, cp);
			double k2 = MK.int_2(u, cp);
			double k3 = MK.int_3(u, cp);
			now->par[0].I_u += mu_ex_p * (k2 / k1) * u1 / u;
			now->par[0].I_v += mu_ex_p * (k2 / k1) * (u2 * cos(alpha) + u3 * sin(alpha)) / u;
			now->par[0].I_T += mu_ex_p * (-0.5 * k3 / k1 + k2 / k1 * skalar / u);
		}

		//drob = min(fabs(polar_angle(y_0, z_0) - polar_angle(y_0 + t_ex * Vy, z_0 + t_ex * Vz)) / 0.017 + 5.0, 80.0);

		if (true)
		{
			for (double tt = 0.0; tt < t_ex * 0.99; tt += t_ex / drob)
			{
				alpha = polar_angle(y_0 + (tt + t_ex / (2.0 * drob)) * Vy, z_0 + (tt + t_ex / (2.0 * drob)) * Vz);

				now->par[0].F_v += (t_ex / drob) * (Vy * cos(alpha) + Vz * sin(alpha)) * mu;
				now->par[0].H_v[dfg] += (t_ex / drob) * (Vy * cos(alpha) + Vz * sin(alpha)) * mu;
			}
		}
		else
		{
			now->par[0].F_v += t_ex * (Vy * cos(alpha) + Vz * sin(alpha)) * mu;
			now->par[0].H_v[dfg] += t_ex * (Vy * cos(alpha) + Vz * sin(alpha)) * mu;
		}

		//now->par[0].I_u += -mu_ex * uz_M * (sigma2(uz_M, cp) / sigma2(uz, cp)) * u1 / u;
		//now->par[0].I_v += -mu_ex * uz_M * (sigma2(uz_M, cp) / sigma2(uz, cp)) * (u2 * cos(alpha) + u3 * sin(alpha)) / u;
		//now->par[0].I_T += mu_ex * (-0.25 * (3.0 * kv(cp) + 2.0 * kv(u)) * (uz_E / uz) * (sigma2(uz_E, cp) / sigma2(uz, cp)) - //
		//	uz_M * (sigma2(uz_M, cp) / sigma2(uz, cp)) * skalar / u);

		/*now->par[0].II_u += -mu_ex * u1;
		now->par[0].II_v += -mu_ex * (u2 * cos(alpha) + u3 * sin(alpha));
		now->par[0].II_T += mu_ex * 0.5 * (kvv(Vx, Vy, Vz) - kvv(vx, vy, 0.0));*/

		now->mut.unlock();


		double Ur, Uthe, Uphi, Vr, Vthe, Vphi;
		spherical_skorost(y_ex, z_ex, x_ex, vy * cos(alpha), vy * sin(alpha), vx, Ur, Uphi, Uthe);
		spherical_skorost(y_ex, z_ex, x_ex, Vy, Vz, Vx, Vr, Vphi, Vthe);
		double mu3;


		int stat_zone_ = now->zona;

		if (area != area2)
		{
			stat_zone_ = -1;
		}

		// перезаряжаем с пикапами

		if (now->pui_ == true)
		{
			for (int i = 0; i < now->i_pui; i++)
			{
				bool kj = true;
				mu3 = mu_ex_pui[i];
				double alpha = polar_angle(y_ex, z_ex);
				double aa, bb, cc;

				now->Change_Velosity_PUI(sens, Vx, Vy, Vz, vx, vy * sin(alpha), vy * cos(alpha), aa, bb, cc, this->Nw, this->wmin, this->wmax, i);


				double peregel, time_do_peregel;
				int ii_z;
				int ii_alp;

				if ((aa * x_ex + bb * y_ex + cc * z_ex) > 0.0)
				{
					time_do_peregel = (-aa * x_ex - bb * y_ex - cc * z_ex) / kvv(aa, bb, cc);
					peregel = sqrt(kvv(x_ex + aa * time_do_peregel, y_ex + bb * time_do_peregel, z_ex + cc * time_do_peregel));
					ii_z = this->geo_zones(peregel);                     // Номер зоны перегелия атома
					ii_alp = alpha_zones(x_ex + aa * time_do_peregel, sqrt(kvv(y_ex + bb * time_do_peregel, z_ex + cc * time_do_peregel, 0.0)));
				}
				else
				{
					ii_z = i_z;
					ii_alp = i_alp;
				}

				if (kj == true)
				{
					if (fabs(mu3) >= pred(area2, ii_z, all))
					{
						kj = true;
					}
					else
					{
						if (fabs(mu3) >= sens->MakeRandom() * pred(area2, ii_z, all))
						{
							mu3 = pred(area2, ii_z, all) * sign(mu3);
							kj = true;
						}
						else
						{
							kj = false;
						}
					}
				}

				if (kj == true)
				{
					Fly_exchenge_Imit_Korol_PUI(MK, sens, x_ex, y_ex, z_ex, aa, bb, cc, now, mu3, area2, false, mu_start, ii_z, ii_alp, true, stat_zone_, true);
				}
			}

		}

		// Перезаряжаем с плазмой
		if (area2 == 0 || Ur / cp > 3.0)   // Без геометрического расщепления
		{
		aa:
			bool kj = true;
			mu3 = mu_ex_p;
			double alpha = polar_angle(y_ex, z_ex);
			double aa, bb, cc;  // Декартовы скорости после перезарядки
			vector <double> Wr(1);
			vector <double> Wthe(1);
			vector <double> Wphi(1);
			vector <double> mu_(1);
			MK.Change_Velosity4(sens, Ur / cp, Uthe / cp, Uphi / cp, Vr / cp, Vthe / cp, Vphi / cp, Wr, Wthe, Wphi, mu_, cp, r, 0);
			Wr[0] = Wr[0] * cp;
			Wthe[0] = Wthe[0] * cp;
			Wphi[0] = Wphi[0] * cp;
			dekard_skorost(y_ex, z_ex, x_ex, Wr[0], Wphi[0], Wthe[0], bb, cc, aa);
			double peregel, time_do_peregel;
			int ii_z;
			int ii_alp;

			if (Wr[0] < 0.0)
			{
				//wwt = sqrt(kv(Wphi[0]) + kv(Wthe[0]));
				//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "Letit vniz  " << x_ex * RR_ << " " << rr * RR_ << " " << mu3 << " " << area << " " << area2 << endl;
				//goto aa;
				time_do_peregel = (-aa * x_ex - bb * y_ex - cc * z_ex) / kvv(aa, bb, cc);
				peregel = sqrt(kvv(x_ex + aa * time_do_peregel, y_ex + bb * time_do_peregel, z_ex + cc * time_do_peregel));
				ii_z = this->geo_zones(peregel);                     // Номер зоны перегелия атома
				ii_alp = alpha_zones(x_ex + aa * time_do_peregel, sqrt(kvv(y_ex + bb * time_do_peregel, z_ex + cc * time_do_peregel, 0.0)));
				//kj = false;   // Убиваем траекторию, иначе может испортить статистику
			}
			else
			{
				ii_z = i_z;
				ii_alp = i_alp;
			}

			if (kj == true)
			{
				if (fabs(mu3) >= pred(area2, ii_z, all))
				{
					kj = true;
				}
				else
				{
					if (fabs(mu3) >= sens->MakeRandom() * pred(area2, ii_z, all))
					{
						mu3 = pred(area2, ii_z, all) * sign(mu3);
						kj = true;
					}
					else
					{
						kj = false;
					}
				}
			}

			if (kj == true)
			{
				/*now->mut.lock();
				now->par[0].II_u += mu3 * (Vx - aa);
				now->par[0].II_v += mu3 * ((Vy - bb) * cos(alpha) + (Vz - cc) * sin(alpha));
				now->par[0].II_T += mu3 * 0.5 * (kvv(Vx, Vy, Vz) - kvv(aa, bb, cc));
				now->mut.unlock();*/
				Fly_exchenge_Imit_Korol_PUI(MK, sens, x_ex, y_ex, z_ex, aa, bb, cc, now, mu3, area2, false, mu_start, ii_z, ii_alp, true, stat_zone_, false);
			}
			//this->Fly_exchenge_Imit(MK, sens, x_ex, y_ex, z_ex, aa, bb, cc, now, mu3, -log(1.0 - sens->MakeRandom()),//
			//	0.0, area2, mu_start, to_I, 0);

		}
		else  // делаем геометрическое расщепление
		{
			int I = geo_zones(r, 1.2);  // Число дополнительных траекторий
			if (I > to_I)  // Для того чтобы не расщеплять атом, который летел в to_I область на атомы из более крупных областей
			{
				I = to_I;
			}

			vector <double> Wr(I + 1);
			vector <double> Wthe(I + 1);
			vector <double> Wphi(I + 1);
			vector <double> mu_(I + 1);

			//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< " C " << endl;
			bool bbb = MK.Change_Velosity4(sens, Ur / cp, Uthe / cp, Uphi / cp, Vr / cp, Vthe / cp, Vphi / cp, Wr, Wthe, Wphi, mu_, cp, r, I);
			//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< " D " << endl;

			for (int i = 0; i <= I; i++)
			{
				Wr[i] = Wr[i] * cp;
				Wthe[i] = Wthe[i] * cp;
				Wphi[i] = Wphi[i] * cp;
			}

			double aa, bb, cc, time_do_peregel, peregel;
			int ii_z, ii_alp;
			bool kj = true;
			//now->par[0].num_atoms++;
			for (int i = 0; i < I; i++)
			{
				dekard_skorost(y_ex, z_ex, x_ex, Wr[i], Wphi[i], Wthe[i], bb, cc, aa);
				mu3 = mu_[i] * mu_ex_p;
				time_do_peregel = (-aa * x_ex - bb * y_ex - cc * z_ex) / kvv(aa, bb, cc);
				peregel = sqrt(kvv(x_ex + aa * time_do_peregel, y_ex + bb * time_do_peregel, z_ex + cc * time_do_peregel));
				ii_z = this->geo_zones(peregel);                     // Номер зоны перегелия атома
				ii_alp = alpha_zones(x_ex + aa * time_do_peregel, sqrt(kvv(y_ex + bb * time_do_peregel, z_ex + cc * time_do_peregel, 0.0)));

				if (fabs(mu3) >= pred(area2, ii_z, all))
				{
					kj = true;
				}
				else
				{
					if (fabs(mu3) >= sens->MakeRandom() * pred(area2, ii_z, all))
					{
						mu3 = pred(area2, ii_z, all);
						kj = true;
					}
					else
					{
						kj = false;
					}
				}
				if (kj == true)
				{
					/*now->mut.lock();
					now->par[0].II_u += mu3 * (Vx - aa);
					now->par[0].II_v += mu3 * ((Vy - bb) * cos(alpha) + (Vz - cc) * sin(alpha));
					now->par[0].II_T += mu3 * 0.5 * (kvv(Vx, Vy, Vz) - kvv(aa, bb, cc));
					now->mut.unlock();*/
					Fly_exchenge_Imit_Korol_PUI(MK, sens, x_ex, y_ex, z_ex, aa, bb, cc, now, mu3, area2, false, mu_start, ii_z, ii_alp, true, stat_zone_, false);
				}
			}

			if (bbb == true)  // Если вообще нужно запускать основной атом
			{
				mu3 = mu_[I] * mu_ex_p;

				// Нужно определить номер основного атома (от его перегелия)
				dekard_skorost(y_ex, z_ex, x_ex, Wr[I], Wphi[I], Wthe[I], bb, cc, aa);
				if (Wr[I] < 0.0)
				{
					time_do_peregel = (-aa * x_ex - bb * y_ex - cc * z_ex) / kvv(aa, bb, cc);
					peregel = sqrt(kvv(x_ex + aa * time_do_peregel, y_ex + bb * time_do_peregel, z_ex + cc * time_do_peregel));
					ii_z = this->geo_zones(peregel);                     // Номер зоны перегелия атома
					ii_alp = alpha_zones(x_ex + aa * time_do_peregel, sqrt(kvv(y_ex + bb * time_do_peregel, z_ex + cc * time_do_peregel, 0.0)));
				}
				else
				{
					ii_z = i_z;
					ii_alp = i_alp;
				}

				if (fabs(mu3) >= pred(area2, ii_z, all))
				{
					kj = true;
				}
				else
				{
					if (fabs(mu3) >= sens->MakeRandom() * pred(area2, ii_z, all))
					{
						mu3 = pred(area2, ii_z, all) * sign(mu3);
						kj = true;
					}
					else
					{
						kj = false;
					}
				}

				if (kj == true)
				{
					/*now->mut.lock();
					now->par[0].II_u += mu3 * (Vx - aa);
					now->par[0].II_v += mu3 * ((Vy - bb) * cos(alpha) + (Vz - cc) * sin(alpha));
					now->par[0].II_T += mu3 * 0.5 * (kvv(Vx, Vy, Vz) - kvv(aa, bb, cc));
					now->mut.unlock();*/
					Fly_exchenge_Imit_Korol_PUI(MK, sens, x_ex, y_ex, z_ex, aa, bb, cc, now, mu3, area2, true, mu_start, ii_z, ii_alp, false, stat_zone_, false);
				}
			}

		}


		// Здесь можно добавить рулетку для оставшегося неперезаряженного атома
		if (fabs(mu2) < pred(area, to_I, all))
		{
			if (fabs(mu2) >= sens->MakeRandom() * pred(area, to_I, all))
			{
				mu2 = pred(area, to_I, all) * sign(mu2);
			}
			else
			{
				return;
			}
		}

	}


	// этот блок кода будет работать и для частиц, которые не попали в предыдущий цикл
	// здесь t2, _ex будут иметь исходные значения, как при инициализации
	int dfg = area;
	if (pui__ == true)
	{
		if (area == 0)
		{
			dfg = 4;
		}
		if (area == 1)
		{
			dfg = 5;
		}
	}

	alpha = polar_angle(y_ex + 0.5 * t2 * Vy, z_ex + 0.5 * t2 * Vz);
	now->mut.lock();
	now->par[0].F_n += t2 * mu2;
	now->par[0].F_u += t2 * Vx * mu2;
	now->par[0].F_T += t2 * kvv(Vx, Vy, Vz) * mu2;
	now->par[0].H_n[dfg] += t2 * mu2;
	now->par[0].H_u[dfg] += t2 * Vx * mu2;
	now->par[0].H_T[dfg] += t2 * kvv(Vx, Vy, Vz) * mu2;


	int ghdf = min((int)(u * n_S / max_S), n_S - 1);
	now->S_p[ghdf] += t2 * mu2 * (kappa_sum / time);
	now->S_m[ghdf] += t2 * mu2;

	double Vrr = (Vy * cos(alpha) + Vz * sin(alpha));
	now->par[0].H_uu[dfg] += t2 * Vx * Vx * mu2;
	now->par[0].H_uv[dfg] += t2 * Vrr * Vx * mu2;
	now->par[0].H_vv[dfg] += t2 * Vrr * Vrr * mu2;
	now->par[0].H_uuu[dfg] += t2 * Vx * Vx * Vx * mu2;

	//drob = min(fabs(polar_angle(y_ex, z_ex) - polar_angle(y_ex + t2 * Vy, z_ex + t2 * Vz)) / 0.017 + 5.0, 80.0);

	if (true)
	{
		for (double tt = 0.0; tt < t2 * 0.99; tt += t2 / drob)  // было 20 и норм работало
		{
			alpha = polar_angle(y_ex + (tt + t2 / (2.0 * drob)) * Vy, z_ex + (tt + t2 / (2.0 * drob)) * Vz);

			now->par[0].F_v += (t2 / drob) * (Vy * cos(alpha) + Vz * sin(alpha)) * mu2;
			now->par[0].H_v[dfg] += (t2 / drob) * (Vy * cos(alpha) + Vz * sin(alpha)) * mu2;
		}
	}
	else
	{
		now->par[0].F_v += t2 * (Vy * cos(alpha) + Vz * sin(alpha)) * mu2;
		now->par[0].H_v[dfg] += t2 * (Vy * cos(alpha) + Vz * sin(alpha)) * mu2;
	}
	now->mut.unlock();

	double mu3 = mu2;
	/*double X = x_ex + t2 * Vx;
	double Y = y_ex + t2 * Vy;
	double Z = z_ex + t2 * Vz;*/
	double X = x_0 + b * Vx;
	double Y = y_0 + b * Vy;
	double Z = z_0 + b * Vz;

	Cell* next = nullptr;             // В какую ячейку попадаем следующую
	double xk, yk;
	xk = x_0 + b * Vx;
	yk = sqrt(kv(y_0 + b * Vy) + kv(z_0 + b * Vz));
	double r_k = sqrt(kv(xk) + kv(yk));
	double r_0 = sqrt(kv(x_0) + kv(y_0) + kv(z_0));

	/*if (yk < 5.0)
	{
		cout << "Setka.cpp    " << "Fly_ex_Korol   "<< b << "  " << xk << endl;
	}*/

	if (slay == true)  //Если надо убить траекторию
	{
		return;
	}

	if (xk < Left_ + 2.0 / RR_)  // ЭТОГО НЕ БЫЛО
	{
		return;
	}

	if (yk > R5_ - 2.0 / RR_)  // ЭТОГО НЕ БЫЛО
	{
		return;
	}

	if (r_k >= Rmax_ - 3.0 / RR_ && xk >= -0.5 / RR_)
	{
		return;
		next = nullptr;
	}

	if (next_mb != nullptr)
	{
		if (next_mb->belong(xk, yk) == true)
		{
			next = next_mb;
		}
	}

	if (next == nullptr)
	{
		for (auto& i : now->Grans)
		{
			if (i->Sosed != nullptr)
			{
				if (i->Sosed->belong(xk, yk) == true)
				{
					next = i->Sosed;
				}
			}
		}
		if (next == nullptr)
		{
			for (auto& k : now->Grans)
			{
				if (k->Sosed != nullptr)
				{
					for (auto& i : k->Sosed->Grans)
					{
						if (i->Sosed != nullptr)
						{
							if (i->Sosed->belong(xk, yk) == true)
							{
								next = i->Sosed;
							}
						}
					}
				}
			}
		}
	}

	if (next == nullptr)
	{
		for (auto& i : this->All_Cells_zero)
		{
			if (i->belong(xk, yk) == true)
			{
				next = i;
			}
		}
	}

	if (next == nullptr)
	{
		this->Smut.lock();
		cout << "Setka.cpp    " << "Fly_ex_Korol   " << "Long   find  " << endl;
		cout << "Setka.cpp    " << "Fly_ex_Korol  end " << xk << " " << yk << endl;
		cout << "Setka.cpp    " << "Fly_ex_Korol  start " << x_0 << " " << y_start << endl;
		this->Smut.unlock();
		//exit(-10);
		//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< x_0 << " " << sqrt(kvv(0.0, y_0, z_0)) << endl;
		//t_per = 100000000.0;
		//if (true)
		//{
		//	double A, B, C;
		//	for (auto& i : now->Grans)
		//	{
		//		double xx, yy;
		//		i->Get_Center(xx, yy);
		//		cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "centr = " << xx << " " << yy << endl;
		//		cout << "Setka.cpp    " << "Fly_ex_Korol   "<< t_per << endl;
		//		if (i->type == Axis)
		//		{
		//			continue;
		//		}
		//		A = i->aa;
		//		B = i->bb;
		//		C = i->cc;
		//		if (fabs(B) < 0.00001)
		//		{
		//			tt1 = (-C / A - x_0) / Vx;
		//			cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "TT1 = " << tt1 << endl;
		//			if (std::fpclassify(tt1) == FP_NORMAL && tt1 > 0.0)
		//			{
		//				if (t_per > tt1)
		//				{
		//					if (now->belong(x_0 + 1.02 * tt1 * Vx, sqrt(kv(y_0 + 1.02 * tt1 * Vy) + kv(z_0 + 1.02 * tt1 * Vz))) == false)
		//					{
		//						t_per = tt1;
		//						next_mb = i->Sosed;
		//					}
		//				}
		//			}
		//			continue;
		//		}
		//		cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "1 - " << A << " " << B << " " << C << endl;
		//		//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "2 - " << x_0 << " " << y_0 << " " << z_0 << endl;
		//		//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "3 - " << Vx << " " << Vy << " " << Vz << endl;
		//		a1 = 4.0 * kv(A * Vx * (C + A * x_0) - kv(B) * (Vy * y_0 + Vz * z_0)) - 4.0 * (kv(A * Vx) - //
		//			B * B * (kv(Vy) + kv(Vz))) * (kv(C + A * x_0) - B * B * (kv(y_0) + kv(z_0)));
		//		a2 = -A * Vx * (C + A * x_0) + B * B * (Vy * y_0 + Vz * z_0);
		//		a3 = A * A * Vx * Vx - B * B * (Vy * Vy + Vz * Vz);
		//		cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "4 - " << a1 << " " << a2 << " " << a3 << endl;
		//		if (fabs(a1) < 0.000001)
		//		{
		//			a1 = 0.0;
		//		}
		//		if (a1 >= 0.0)
		//		{
		//			tt1 = (a2 + 0.5 * sqrt(a1)) / a3;
		//			tt2 = (a2 - 0.5 * sqrt(a1)) / a3;
		//			//cout << "Setka.cpp    " << "Fly_ex_Korol   "<< A << " " << B << " " << C << endl;
		//			cout << "Setka.cpp    " << "Fly_ex_Korol   "<< tt1 << " " << tt2 << endl;
		//			if (std::fpclassify(tt1) == FP_NORMAL && tt1 > 0.0)
		//			{
		//				if (t_per > tt1)
		//				{
		//					if (now->belong(x_0 + 1.03 * tt1 * Vx, sqrt(kv(y_0 + 1.03 * tt1 * Vy) + kv(z_0 + 1.03 * tt1 * Vz))) == false)
		//					{
		//						cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "A " << endl;
		//						t_per = tt1;
		//					}
		//				}
		//			}
		//			if (std::fpclassify(tt2) == FP_NORMAL && tt2 > 0.0)
		//			{
		//				if (t_per > tt2)
		//				{
		//					if (now->belong(x_0 + 1.03 * tt2 * Vx, sqrt(kv(y_0 + 1.03 * tt2 * Vy) + kv(z_0 + 1.03 * tt2 * Vz))) == false)
		//					{
		//						cout << "Setka.cpp    " << "Fly_ex_Korol   "<< "B " << endl;
		//						t_per = tt2;
		//					}
		//				}
		//			}
		//		}
		//	}
		//	cout << "Setka.cpp    " << "Fly_ex_Korol   "<< t_per << endl;
		//	xk = x_0 + t_per * Vx;
		//	yk = sqrt(kv(y_0 + t_per * Vy) + kv(z_0 + t_per * Vz));
		//	cout << "Setka.cpp    " << "Fly_ex_Korol   "<< xk << " " << yk << endl;

		//	a = t_per * 0.99;
		//	b = t_per * 1.01;
		//}




		//exit(-1);
		for (auto& i : this->All_Cells)
		{
			if (i->belong(xk, yk) == true)
			{
				next = i;
				break;
			}
		}
	}

	if (next != nullptr)
	{
		// Для подсчёта статистики весов
		if (mu_statistic)
		{
			if (now->zona != zon_stat)
			{
				Mu_stat[area][now->zona] += mu;
				I_stat[area][now->zona] ++;
				zon_stat = now->zona;
			}
		}

		Fly_exchenge_Imit_Korol_PUI(MK, sens, X, Y, Z, Vx, Vy, Vz, next, mu3, area, false, mu_start, to_I, to_J, georaschep, zon_stat, pui__);
	}
	else
	{

		if (xk > -1000.0 / RR_ && xk < 1000.0 / RR_ && yk < 1600.0 / RR_)
		{
			cout << "Setka.cpp    " << "Fly_ex_Korol   " << "ERROR  8670  poteryal atom" << endl;
			cout << "Setka.cpp    " << "Fly_ex_Korol   " << "SDFSDF    " << xk << " " << yk << endl;
		}
	}

	return;
}

void Setka::Fly_exchenge_Imit(MKmethod& MK, Sensor* sens, double x_0, double y_0, double z_0, double Vx, double Vy, double Vz, Cell* now, double mu, double KSI, //
	double I_do, int area, const double& mu_start, int to_I = 0, int iii = 0)
{
	//cout << "Setka.cpp    " << "Start  " << iii <<  endl;
	double normV = sqrt(kvv(Vx, Vy, Vz));
	double dt = now->L / (normV); 

	/*double Vr, Vthe, Vphi;
	spherical_skorost(y_0, z_0, x_0, Vy, Vz, Vx, Vr, Vphi, Vthe);
	if (Vr > 0.0 && to_I < now->zona)
	{
		to_I = now->zona;
		if (mu < Mu[area][now->zona])
		{
			if (mu >= sens->MakeRandom() * Mu[area][now->zona])
			{
				mu = Mu[area][now->zona];
			}
			else
			{
				return;
			}
		}
	}*/


	//cout << "Setka.cpp    " << "V = " << Vx << " " << Vy << " " << Vz << endl;
	//cout << "Setka.cpp    " << "A" << endl;
	double a = 0.0, b = dt;
	double t_per = 1000000000.0;

	/*double t_per = (-Vy * y_0 - Vz * z_0) / (kv(Vy) + kv(Vz));
	if (t_per > 0.0000001)
	{
		if (now->belong(x_0 + t_per * Vx, sqrt(kv(y_0 + t_per * Vy) + kv(z_0 + t_per * Vz))) == true)
		{
			t_per = -1.0;
		}
	}
	else
	{
		t_per = -1.0;
	}*/

	double tt1, tt2;
	double a1, a2, a3;
	double k;

	int lk = 0;
	if (true)
	{
		a = 0.0;
		b = dt;
		t_per = 1000000000.0;
		for (auto& i : now->Grans)
		{
			a1 = -kv(i->a * Vx) * (kv(i->a * Vx) - kv(Vy) - kv(Vz)) * kv(Vz * y_0 - Vy * z_0);
			a2 = (kv(i->a * Vx) - kv(Vy) - kv(Vz)) * (Vz * z_0 + Vy * y_0);
			a3 = (kv(Vy) + kv(Vz)) * (-kv(i->a * Vx) + kv(Vy) + kv(Vz));
			if (a1 >= 0)
			{
				tt1 = (sqrt(a1) + a2) / a3;
				tt2 = (-sqrt(a1) + a2) / a3;
				if (std::fpclassify(tt1) == FP_NORMAL && tt1 > 0.0)
				{
					if (now->belong(x_0 + tt1 * Vx, sqrt(kv(y_0 + tt1 * Vy) + kv(z_0 + tt1 * Vz))) == false)
					{
						if (t_per > tt1)
						{
							t_per = tt1;
						}
					}
				}
				if (std::fpclassify(tt2) == FP_NORMAL && tt2 > 0.0)
				{
					if (now->belong(x_0 + tt2 * Vx, sqrt(kv(y_0 + tt2 * Vy) + kv(z_0 + tt2 * Vz))) == false)
					{
						if (t_per > tt2)
						{
							t_per = tt2;
						}
					}
				}
			}
		}

		while (now->belong(x_0 + b * Vx, sqrt(kv(y_0 + b * Vy) + kv(z_0 + b * Vz))) == true) 
		{
			lk++;
			if (lk > 1000)
			{
				cout << "Setka.cpp    " << "8055 dt = " << dt << " " << x_0 << " " << sqrt(kv(y_0) + kv(z_0)) << endl;
				cout << "Setka.cpp    " << normV << " " << Vx << " " << Vy << " " << Vz << endl;
				exit(-1);
			}
			b = b * 1.1;
		}

		if (t_per < b)
		{
			b = t_per;
		}

		k = (a + b) / 2.0;
		while (fabs(a - b) * normV > now->L / 200.0)
		{
			if (now->belong(x_0 + k * Vx, sqrt(kv(y_0 + k * Vy) + kv(z_0 + k * Vz))) == true)
			{
				a = k;
			}
			else
			{
				b = k;
			}
			k = (a + b) / 2.0;
		}
	}

	double time = b; 

	/*if (iii == 219)
	{
		cout << "Setka.cpp    " << "a, b, time, dt,  a, now->L  =   " << a << "   " << b << " " << time << " " << dt << " " << now->L << " " << endl;
	}*/
	//cout << "Setka.cpp    " << a << " " << b << " " << time << " " << dt << " " << now->L << " " << endl;

	double uz, uz_M, uz_E, t1;// y, z, r;
	double cp = sqrt(now->par[0].p / now->par[0].ro);
	bool change_was = false;
	double vx = now->par[0].u;
	double vy = now->par[0].v;
	double ro = now->par[0].ro;

	double t_ex = 0.0;
	double mu2 = mu;
	double mu3 = 0.0;
	double nu_ex, u, alpha;
	double t2 = time;

	double x_ex = x_0, y_ex = y_0, z_ex = z_0;
	double u1, u2, u3;

	int nk = 100;
	double sig = 0.0;
	double l = 0.0;
	double I = I_do;
	l = sqrt(kvv(time / nk * Vx, time / nk * Vy, time / nk * Vz));
	for (double tt = 0.0; tt < time * 0.9999; tt += time / nk)
	{
		alpha = polar_angle(y_0 + (tt + time / (2.0 * nk)) * Vy, z_0 + (tt + time / (2.0 * nk)) * Vz);
		u = sqrt(kvv(Vx - vx, Vy - vy * cos(alpha), Vz - vy * sin(alpha)));
		uz = Velosity_1(u, cp);
		sig = sqrt(kvv(Vx, Vy, Vz)) / ((1.0 / Kn_) * ro * uz * sigma(uz));
		I += l / sig;
	}

	//double l = sqrt(kvv(time * Vx, time * Vy, time * Vz));
	//alpha = polar_angle(y_0, z_0);
	//u = sqrt(kvv(Vx - vx, Vy - vy * cos(alpha), Vz - vy * sin(alpha)));
	//uz = Velosity_1(u, cp);
	////double IS = MK.Int_cp_1(u / cp);
	////sig = sqrt(kvv(Vx, Vy, Vz)) / ((1.0 / Kn_) * ro * IS);
	//sig = sqrt(kvv(Vx, Vy, Vz)) / ((1.0 /Kn_) * ro * uz * sigma(uz));
	//double I = I_do + l / sig;



	if (I < KSI)  
	{
		I_do = I;
	}
	else  
	{
		//cout << "Setka.cpp    " << "C" << endl;

		double ksi = (KSI - I_do) * sig;
		t_ex = ksi / sqrt(kvv(Vx, Vy, Vz));
		I_do = 0.0;
		//KSI = -log(1.0 - sens->MakeRandom());
		alpha = polar_angle(y_0 + 0.5 * t_ex * Vy, z_0 + 0.5 * t_ex * Vz);  //      

		x_ex = x_0 + t_ex * Vx;  //  
		y_ex = y_0 + t_ex * Vy;
		z_ex = z_0 + t_ex * Vz;
		if (now->belong(x_ex, sqrt(kv(y_ex) + kv(z_ex))) == false)
		{
			//cout << "Setka.cpp    " << "TUT" << endl;
			x_ex = x_0 + a * Vx;  //  
			y_ex = y_0 + a * Vy;
			z_ex = z_0 + a * Vz;
			t_ex = a;
			alpha = polar_angle(y_0 + 0.5 * (t_ex)*Vy, z_0 + 0.5 * (t_ex)*Vz);
		}


		//      
		u = sqrt(kvv(Vx - vx, Vy - vy * cos(alpha), Vz - vy * sin(alpha)));
		uz = Velosity_1(u, cp);
		uz_M = Velosity_2(u, cp) / (uz * kv(cp) * cp * pi * sqrtpi_);
		uz_E = Velosity_3(u, cp);
		u1 = vx - Vx;
		u2 = vy * cos(alpha) - Vy;
		u3 = vy * sin(alpha) - Vz;
		double skalar = Vx * u1 + Vy * u2 + Vz * u3;

		now->mut.lock();
		now->par[0].F_n += t_ex * mu;
		now->par[0].F_u += t_ex * Vx * mu;
		now->par[0].F_v += t_ex * (Vy * cos(alpha) + Vz * sin(alpha)) * mu;
		now->par[0].F_T += t_ex * kvv(Vx, Vy, Vz) * mu;
		now->par[0].H_n[area] += t_ex * mu;
		now->par[0].H_u[area] += t_ex * Vx * mu;
		now->par[0].H_v[area] += t_ex * (Vy * cos(alpha) + Vz * sin(alpha)) * mu;
		now->par[0].H_T[area] += t_ex * kvv(Vx, Vy, Vz) * mu;

		now->par[0].I_u += -mu * uz_M * u1 / u;
		now->par[0].I_v += -mu * uz_M * (u2 * cos(alpha) + u3 * sin(alpha)) / u;
		now->par[0].I_T += mu * (-0.25 * (3.0 * kv(cp) + 2.0 * kv(u)) * (uz_E / uz) - //
			uz_M * skalar / u);

		//now->par[0].I_u += -mu * uz_M * (sigma2(uz_M, cp) / sigma2(uz, cp)) * u1 / u;
		//now->par[0].I_v += -mu * uz_M * (sigma2(uz_M, cp) / sigma2(uz, cp)) * (u2 * cos(alpha) + u3 * sin(alpha)) / u;
		//now->par[0].I_T += mu * (-0.25 * (3.0 * kv(cp) + 2.0 * kv(u)) * (uz_E / uz) * (sigma2(uz_E, cp) / sigma2(uz, cp)) - //
		//	uz_M * (sigma2(uz_M, cp) / sigma2(uz, cp)) * skalar / u);

		/*now->par[0].II_u += -mu * u1;
		now->par[0].II_v += -mu * (u2 * cos(alpha) + u3 * sin(alpha));
		now->par[0].II_T += mu * 0.5 * (kvv(Vx, Vy, Vz) - kvv(vx, vy, 0.0));*/
		now->mut.unlock();

		//cout << "Setka.cpp    " << "SDDD" << endl;
		int area2;
		if (now->type == C_5)
		{
			area2 = 3;
		}
		else if (now->type == C_4)
		{
			area2 = 2;
		}
		else if (now->type == C_1 || now->type == C_centr)
		{
			area2 = 0;
		}
		else
		{
			area2 = 1;
		}

		//cout << "Setka.cpp    " << "SDDD 2" << endl;
		double r = sqrt(kvv(x_ex, y_ex, z_ex));
		double rr = sqrt(kvv(0.0, y_ex, z_ex));
		double Ur, Uthe, Uphi, Vr, Vthe, Vphi;
		//cout << "Setka.cpp    " << "E" << endl;
		spherical_skorost(y_ex, z_ex, x_ex, vy * cos(alpha), vy * sin(alpha), vx, Ur, Uphi, Uthe);
		spherical_skorost(y_ex, z_ex, x_ex, Vy, Vz, Vx, Vr, Vphi, Vthe);
		//double uu = sqrt(kvv(Vr - Ur, Vthe - Uthe, Vphi - Uphi)) / cp;
		//double UUr = (Ur - Vr) / cp;
		//double UUa = sqrt(kvv(0.0, Vthe - Uthe, Vphi - Uphi)) / cp;

		//f_way << x_ex << " " << sqrt(kvv(z_ex, y_ex, 0.0)) << " " << 15 << endl;
		//f_num++;


		if (r > Rmax_ && x_ex > 0.001 / RR_)
		{
			return;
		}

		int i_z = geo_zones(r);

		//  uu > 5    
		//  (UUr >= -5.0 && UUr <= 0.0)


		if (true)//(area2 == 0)   //  
		{
			int num = 1;
			bool kj = true;
			mu3 = mu;

			double alpha = polar_angle(y_ex, z_ex);
			double aa, bb, cc;
			vector <double> Wr(1);
			vector <double> Wthe(1);
			vector <double> Wphi(1);
			vector <double> mu_(1);

			MK.Change_Velosity(sens, Ur / cp, Uthe / cp, Uphi / cp, Vr / cp, Vthe / cp, Vphi / cp, Wr, Wthe, Wphi, mu_, cp, r, 0);
			Wr[0] = Wr[0] * cp;
			Wthe[0] = Wthe[0] * cp;
			Wphi[0] = Wphi[0] * cp;

			/*if (area2 == 1 && mu2 > 0.01)
			{
				cout << "Setka.cpp    " << mu3 << " " << area << " " << Vr << " " << Wr[0] << " " << x_ex << " " << sqrt(kv(y_ex) + kv(z_ex)) << endl;
			}*/

			double wwt, peregel;
			int ii_z;
			bool rulet = true;

			if (Wr[0] < 0.0)
			{
				//cout << "Setka.cpp    " << Wr[0] << " " << area << " " << mu << endl;
				rulet = false;
				wwt = sqrt(kv(Wphi[0]) + kv(Wthe[0]));
				peregel = r * wwt / (sqrt(kv(wwt) + kv(Wr[0])));
				ii_z = this->geo_zones(peregel);                     //    
			}
			else
			{
				rulet = true;
				ii_z = i_z;
				//if (ii_z == 2 || ii_z == 1)  //        
				//{
				//	ii_z = 0;
				//}
			}

			if (true)
			{
				double ko = 1.0;

				/*if (area2 == 0 || area2 == 1)
				{
					ko = 0.001;
				}*/

				if (mu3 >= Mu[area2][ii_z] * ko)
				{
					kj = true;
				}
				else
				{
					if (mu3 >= sens->MakeRandom() * Mu[area2][ii_z] * ko)
					{
						mu3 = Mu[area2][ii_z] * ko;
						kj = true;
					}
					else
					{
						kj = false;
					}
				}
			}

			if (kj == true)
			{
				dekard_skorost(y_ex, z_ex, x_ex, Wr[0], Wphi[0], Wthe[0], bb, cc, aa);
				KSI = -log(1.0 - sens->MakeRandom());
				Fly_exchenge_Imit(MK, sens, x_ex, y_ex, z_ex, aa, bb, cc, now, mu3, KSI, I_do, area2, mu_start, ii_z, iii);
			}

			return;
		}
		else   //  
		{
			double alpha = polar_angle(y_ex, z_ex);

			int I = geo_zones(r, 1.5);  //   

			if (I > to_I + 1)   //   ,         ,   
			{
				I = to_I + 1;
			}

			if (I > 1)
			{
				I = 1;
			}

			int num = 1;
			/*if (rr < 50.0)
			{
				num = 3;
			}*/
			mu = mu / num;
			for (int kkk = 0; kkk < num; kkk++)
			{
				vector <double> Wr(I + 1);
				vector <double> Wthe(I + 1);
				vector <double> Wphi(I + 1);
				vector <double> mu_(I + 1);

				//bool bbb = true;
				bool bbb = MK.Change_Velosity(sens, Ur / cp, Uthe / cp, Uphi / cp, Vr / cp, Vthe / cp, Vphi / cp, Wr, Wthe, Wphi, mu_, cp, r, I, x_ex, y_ex, z_ex);
				//Change_Velosity_Split(sens, Ur / cp, Uthe / cp, Uphi / cp, Vr / cp, Vthe / cp, Vphi / cp, Wr, Wthe, Wphi, mu_, cp, r, I, x_ex, y_ex, z_ex);
				double aa, bb, cc;
				for (int i = 0; i <= I; i++)
				{
					Wr[i] = Wr[i] * cp;
					Wthe[i] = Wthe[i] * cp;
					Wphi[i] = Wphi[i] * cp;
				}


				//         
				/*if (area == 3 && I >= 5 && bbb == true)
				{
					cout << "Setka.cpp    " << "to_I = " << to_I << endl;
					for (int i = 0; i < I; i++)
					{
						dekard_skorost(y_ex, z_ex, x_ex, Wr[i], Wphi[i], Wthe[i], bb, cc, aa);
						for (int ij = 0; ij < 27000; ij++)
						{
							cout << "Setka.cpp    " << x_ex + ij * aa << " " << sqrt(kvv(y_ex + ij * bb, z_ex + ij * cc, 0.0)) << " " << i << endl;
						}
					}
					if (bbb == true)
					{
						dekard_skorost(y_ex, z_ex, x_ex, Wr[I], Wphi[I], Wthe[I], bb, cc, aa);
						for (int ij = 0; ij < 30000; ij++)
						{
							cout << "Setka.cpp    " << x_ex + ij * aa << " " << sqrt(kvv(y_ex + ij * bb, z_ex + ij * cc, 0.0)) << " " << I << endl;
						}
					}
					exit(-1);
				}*/

				bool kj = true;
				double ko = 1.0;

				/*if (area2 == 0 || area2 == 1)
				{
					ko = 0.001;
				}*/

				/*if (area2 == 1)
				{
					for (int i = 0; i < I; i++)
					{
						cout << "Setka.cpp    " << mu_[i] << endl;
					}
					exit(-1);
				}*/

				for (int i = 0; i < I; i++)
				{
					kj = true;
					mu3 = mu_[i] * mu;

					if (true)//(i != to_I)
					{
						if (fabs(mu3) >= Mu[area2][i] * ko / num)
						{
							kj = true;
						}
						else
						{
							if (fabs(mu3) >= sens->MakeRandom() * Mu[area2][i] * ko / num)
							{
								mu3 = Mu[area2][i] * ko * sign(mu3) / num;
								kj = true;
							}
							else
							{
								kj = false;
							}
						}
					}

					if (kj == true)
					{
						dekard_skorost(y_ex, z_ex, x_ex, Wr[i], Wphi[i], Wthe[i], bb, cc, aa);
						KSI = -log(1.0 - sens->MakeRandom());
						Fly_exchenge_Imit(MK, sens, x_ex, y_ex, z_ex, aa, bb, cc, now, mu3, KSI, 0.0, area2, mu_start / num, i, iii);
						//Fly_exchenge_Imit_Korol(MK, sens, x_ex, y_ex, z_ex, aa, bb, cc, now, mu3, area2, true, mu_start, i);
					}
				}

				double wwt, peregel;
				int ii_z;

				if (Wr[I] < 0.0)
				{
					wwt = sqrt(kv(Wphi[I]) + kv(Wthe[I]));
					peregel = r * wwt / (sqrt(kv(wwt) + kv(Wr[I])));
					ii_z = this->geo_zones(peregel);                     //    
				}
				else
				{
					ii_z = i_z;
				}

				kj = true;
				if (bbb == true)
				{

					mu3 = mu_[I] * mu;
					if (true)
					{
						if (fabs(mu3) >= Mu[area2][ii_z] * mu_start * ko / num)
						{
							kj = true;
						}
						else
						{
							if (fabs(mu3) >= sens->MakeRandom() * Mu[area2][ii_z] * mu_start * ko / num)
							{
								mu3 = Mu[area2][ii_z] * mu_start * ko * sign(mu3) / num;
								kj = true;
							}
							else
							{
								kj = false;
							}
						}
					}


					if (kj == true)
					{
						dekard_skorost(y_ex, z_ex, x_ex, Wr[I], Wphi[I], Wthe[I], bb, cc, aa);
						KSI = -log(1.0 - sens->MakeRandom());
						Fly_exchenge_Imit(MK, sens, x_ex, y_ex, z_ex, aa, bb, cc, now, mu3, KSI, 0.0, area2, mu_start / num, ii_z, iii);
						//Fly_exchenge_Imit_Korol(MK, sens, x_ex, y_ex, z_ex, aa, bb, cc, now, mu3, area2, true, mu_start, to_I);
					}
				}
			}

			//dekard_skorost(y_ex, z_ex, x_ex, Wr[I], Wphi[I], Wthe[I], bb, cc, aa);
			//KSI = -log(1.0 - sens->MakeRandom());
			//Fly_exchenge_Imit(sens, x_ex, y_ex, z_ex, aa, bb, cc, now, mu_[I] * mu, KSI, I_do, area, mu_start, to_I, iii);
			return;
		}
	}


	//cout << "Setka.cpp    " << "B1" << endl;
	alpha = polar_angle(y_ex + 0.5 * t2 * Vy, z_ex + 0.5 * t2 * Vz);
	//cout << "Setka.cpp    " << "B2" << endl;
	now->mut.lock();
	now->par[0].F_n += t2 * mu2;
	now->par[0].F_u += t2 * Vx * mu2;
	now->par[0].F_v += t2 * (Vy * cos(alpha) + Vz * sin(alpha)) * mu2;
	now->par[0].F_T += t2 * kvv(Vx, Vy, Vz) * mu2;
	now->par[0].H_n[area] += t2 * mu2;
	now->par[0].H_u[area] += t2 * Vx * mu2;
	now->par[0].H_v[area] += t2 * (Vy * cos(alpha) + Vz * sin(alpha)) * mu2;
	now->par[0].H_T[area] += t2 * kvv(Vx, Vy, Vz) * mu2;
	now->mut.unlock();

	double X = x_ex + t2 * Vx;
	double Y = y_ex + t2 * Vy;
	double Z = z_ex + t2 * Vz;


	Cell* next = nullptr;             //     
	double xk, yk;
	xk = x_0 + b * Vx;
	yk = sqrt(kv(y_0 + b * Vy) + kv(z_0 + b * Vz));
	double r_k = sqrt(kv(xk) + kv(yk));
	double r_0 = sqrt(kv(x_0) + kv(y_0) + kv(z_0));

	//f_way << xk << " " << yk << " " << to_I << endl;
	//f_num++;

	if (r_0 <= Rmax_ && r_k >= Rmax_ && x_0 > 0.0)
	{
		//cout << "Setka.cpp    " << x_0 << " " << sqrt(kvv(0.0, y_0, z_0)) << " " << xk << " " << yk << endl;
		next = nullptr;
		return;
	}
	else if (r_k >= Rmax_ - 0.001 / RR_ && xk >= -0.001 / RR_)
	{
		//cout << "Setka.cpp    " << x_0 << " " << sqrt(kvv(0.0, y_0, z_0)) << " " << xk << " " << yk << endl;
		next = nullptr;
		return;
	}
	else
	{
		for (auto& i : now->Grans)
		{
			if (i->Sosed != nullptr)
			{
				if (i->Sosed->belong(xk, yk) == true)
				{
					next = i->Sosed;
				}
			}
		}
		if (next == nullptr)
		{
			for (auto& k : now->Grans)
			{
				if (k->Sosed != nullptr)
				{
					for (auto& i : k->Sosed->Grans)
					{
						if (i->Sosed != nullptr)
						{
							if (i->Sosed->belong(xk, yk) == true)
							{
								next = i->Sosed;
							}
						}
					}
				}
			}
		}
	}

	if (next == nullptr)
	{
		for (auto& i : this->All_Cells_zero)
		{
			if (i->belong(xk, yk) == true)
			{
				next = i;
			}
		}
	}

	if (next != nullptr)
	{
		Fly_exchenge_Imit(MK, sens, X, Y, Z, Vx, Vy, Vz, next, mu2, KSI, I_do, area, mu_start, to_I, iii);
	}
	else
	{
		if (xk > -1000.0 / RR_ && xk < 1000.0 / RR_ && yk < 200.0 / RR_) // ,           
		{
			cout << "Setka.cpp 17048 f  " << xk << " " << yk << endl;
			cout << Vx << " " << Vy << " " << Vz << endl;
		}
	}

	//cout << "Setka.cpp    " << "Vihod" << endl;
	return;

}

//
//void Setka::Fly_exchenge_Split(Sensor* sens, double x_0, double y_0, double z_0, double Vx, double Vy, double Vz, Cell* now, double mu,//
//	const double& mu_0, bool ExCh, int zone)
//{
//	//cout << "Setka.cpp    " << "NEW" << endl;
//	//cout << "Setka.cpp    "  << x_0 << " " << sqrt(kvv(y_0, z_0, 0.0)) << " " << mu_0 << endl;
//	double f1, f2;
//	now->Get_Center(f1, f2);
//	//cout << "Setka.cpp    " << "Cent = " << f1 << " " << f2 << endl;
//	double normV = sqrt(kvv(Vx, Vy, Vz));
//	double dt = geo_accur * now->L / normV;
//	//cout << "Setka.cpp    " << "dt = " << dt << endl;
//
//
//	// Нахождение времени до попадания в следующую ячейку
//	//cout << "Setka.cpp    " << "V = " << Vx << " " << Vy << " " << Vz << endl;
//	while (now->belong(x_0 + 100.0 * dt * Vx, sqrt(kv(y_0 + 100.0 * dt * Vy) + kv(z_0 + 100.0 * dt * Vz))) == false) // Слишком большой шаг
//	{
//		dt = dt / 2.0;
//		//cout << "Setka.cpp    " << "dt = " << dt << endl;
//	}
//
//	int a = 100;
//	int b = 1000;
//	while (now->belong(x_0 + (1.0 * b) * dt * Vx, sqrt(kv(y_0 + (1.0 * b) * dt * Vy) + kv(z_0 + (1.0 * b) * dt * Vz))) == true)
//	{
//		b = b + 50;
//		//cout << "Setka.cpp    " << "b = " <<  x_0 + (1.0 * b) * dt * Vx << " " << sqrt(kv(y_0 + (1.0 * b) * dt * Vy) + kv(z_0 + (1.0 * b) * dt * Vz)) << " " << b <<  endl;
//	}
//
//	int k = (a + b) / 2;
//	while (k != a)
//	{
//		if (now->belong(x_0 + (1.0 * k) * dt * Vx, sqrt(kv(y_0 + (1.0 * k) * dt * Vy) + kv(z_0 + (1.0 * k) * dt * Vz))) == true)
//		{
//			a = k;
//		}
//		else
//		{
//			b = k;
//		}
//		k = (a + b) / 2;
//		//cout << "Setka.cpp    " << a << " " << b << endl;
//	}
//	
//	Smut.lock();
//	//cout << "Setka.cpp    " << b << endl;
//	mmu1 += 1.0 * b;
//	mn1++;
//	Smut.unlock();
//
//	if (b != a + 1)  // Можно потом удалить проверку
//	{
//		cout << "Setka.cpp    " << "ERRORIUEHUYCEUVDCC W" << endl;
//		exit(-1);
//	}
//
//	double time = dt * (b); // Время нахождения в ячейке
//
//	double uz, uz_M, uz_E, t1;// y, z, r;
//	double cp = sqrt(now->par[0].p / now->par[0].ro);
//	bool change_was = false;
//	double vx = now->par[0].u;
//	double vy = now->par[0].v;
//	double ro = now->par[0].ro;
//
//	double t_ex = 0.0;
//	double mu2 = mu;
//	double nu_ex, u, alpha, mu_ex;
//	double t2 = time;
//
//	double x_ex = x_0, y_ex = y_0, z_ex = z_0;
//	double u1, u2, u3;
//
//	if (ExCh == false) // Если перезарядки в этой ячейке ещё не было
//	{
//		alpha = polar_angle(y_0, z_0);
//		u = sqrt(kvv(Vx - vx, Vy - vy * cos(alpha), Vz - vy * sin(alpha)));
//		uz = Velosity_1(u, cp);
//		nu_ex = (ro * uz * sigma(uz));
//		double kappa = (nu_ex * time) / Kn_;
//		t_ex = -(time / kappa) * log(1.0 - sens->MakeRandom() * (1.0 - exp(-kappa))); // Время до перезарядки
//		t2 = t2 - t_ex;
//		mu_ex = mu * (1.0 - exp(-kappa));
//		mu2 = mu - mu_ex;
//		alpha = polar_angle(y_0 + 0.5 * (t_ex)*Vy, z_0 + 0.5 * (t_ex)*Vz);
//
//
//		x_ex = x_0 + t_ex * Vx;  // Координаты перезарядки
//		y_ex = y_0 + t_ex * Vy;
//		z_ex = z_0 + t_ex * Vz;
//		if (now->belong(x_ex, sqrt(kv(y_ex) + kv(z_ex))) == false)
//		{
//			x_ex = x_0 + (1.0 * a * dt) * Vx;  // Координаты перезарядки
//			y_ex = y_0 + (1.0 * a * dt) * Vy;
//			z_ex = z_0 + (1.0 * a * dt) * Vz;
//			t_ex = (1.0 * a * dt);
//			t2 = time - t_ex;
//			alpha = polar_angle(y_0 + 0.5 * (t_ex)*Vy, z_0 + 0.5 * (t_ex)*Vz);
//		}
//
//		// Считаем источники для частицы до перезарядки
//		u = sqrt(kvv(Vx - vx, Vy - vy * cos(alpha), Vz - vy * sin(alpha)));
//		uz = Velosity_1(u, cp);
//		uz_M = Velosity_2(u, cp) / (uz * kv(cp) * cp * pi * sqrtpi_);
//		uz_E = Velosity_3(u, cp);
//		u1 = vx - Vx;
//		u2 = vy * cos(alpha) - Vy;
//		u3 = vy * sin(alpha) - Vz;
//		double skalar = Vx * u1 + Vy * u2 + Vz * u3;
//
//
//		now->mut.lock();
//		now->par[0].F_n += t_ex * mu;
//		now->par[0].F_u += t_ex * Vx * mu;
//		now->par[0].F_v += t_ex * (Vy * cos(alpha) + Vz * sin(alpha)) * mu;
//		now->par[0].F_T += t_ex * kvv(Vx, Vy, Vz) * mu;
//		now->par[0].I_u += mu * t_ex * uz_M * uz * sigma(uz_M) * u1 / u;
//		now->par[0].I_v += mu * t_ex * uz_M * uz * sigma(uz_M) * (u2 * cos(alpha) + u3 * sin(alpha)) / u;
//		now->par[0].I_T += mu * t_ex * (0.5 * (3.0 * kv(cp) + 2.0 * kv(u)) * uz_E * sigma(uz_E) + 2.0 * uz_M * uz * sigma(uz_M) * skalar / u);
//		now->mut.unlock();
//
//
//		double alpha = polar_angle(y_ex, z_ex);
//
//		double Ur, Uthe, Uphi, Vr, Vthe, Vphi;
//		double r = sqrt(kvv(x_ex, y_ex, z_ex));
//
//		spherical_skorost(y_ex, z_ex, x_ex, vy* cos(alpha), vy* sin(alpha), vx, Ur, Uphi, Uthe);
//		spherical_skorost(y_ex, z_ex, x_ex, Vy, Vz, Vx, Vr, Vphi, Vthe);
//
//		int I = 0;
//		//if (Ur <= 0)
//		//{
//		//	if (r > Ri[I_ - 2])  // Выбираем, сколько расщеплений атома будет.
//		//	{
//		//		I = 4;
//		//	}
//		//	else if (r > Ri[I_ - 3])
//		//	{
//		//		I = 3;
//		//	}
//		//	else if (r > Ri[I_ - 4])
//		//	{
//		//		I = 2;
//		//	}
//		//	else if (r > Ri[I_ - 5])
//		//	{
//		//		I = 1;
//		//	}
//		//	else
//		//	{
//		//		I = 0;
//		//	}
//		//}
//		//cout << "Setka.cpp    " << "I = " << I << " " << Ur << endl;
//		vector <double> Wr(I + 1);
//		vector <double> Wthe(I + 1);
//		vector <double> Wphi(I + 1);
//		vector <double> mu_(I + 1);
//
//		//cout << "Setka.cpp    " << 1 << endl;
//		Change_Velosity_Split(sens, Ur / cp, Uthe / cp, Uphi / cp, Vr / cp, Vthe / cp, Vphi / cp, Wr, Wthe, Wphi, mu_, cp, r, I);
//		//cout << "Setka.cpp    " << 2 << endl;
//
//		for (int i = 0; i <= I; i++)
//		{
//			Wr[i] = Wr[i] * cp;
//			Wthe[i] = Wthe[i] * cp;
//			Wphi[i] = Wphi[i] * cp;
//			//cout << "Setka.cpp    " << mu_[i] << endl;
//		}
//		//exit(-1);
//
//		double aa, bb, cc;
//		bool bol = true;
//		for (int i = 0; i < I; i++)
//		{
//			if (mu_[i] >= 0.1)
//			{
//				bol = true;
//			}
//			else
//			{
//				if (mu_[i] >= sens->MakeRandom() * 0.1)
//				{
//					mu_[i] =0.1;
//					bol = true;
//				}
//				else
//				{
//					bol = false;
//				}
//			}
//			if (bol == true)
//			{
//				dekard_skorost(y_ex, z_ex, x_ex, Wr[i], Wphi[i], Wthe[i], bb, cc, aa);
//				/*for (int k = 0; k < 30000; k++)
//				{
//					cout << "Setka.cpp    " << x_ex + 0.3 * k * aa << " " << sqrt(kv(y_ex + 0.3 * k * bb) + kv(z_ex + 0.3 * k * cc)) << endl;
//				}*/
//				Fly_exchenge_Split(sens, x_ex, y_ex, z_ex, aa, bb, cc, now, mu_[i] * mu_ex, mu_ex, true, i);
//			}
//		}
//		//cout << "Setka.cpp    " << 3 << endl;
//
//		// У основного атома ограничение на вес больше
//		if (mu_[I] >= Mu[I_ - 1])
//		{
//			bol = true;
//		}
//		else
//		{
//			if (mu_[I] >= sens->MakeRandom() * Mu[I_ - 1])
//			{
//				mu_[I] = Mu[I_ - 1];
//				bol = true;
//			}
//			else
//			{
//				bol = false;
//			}
//		}
//		//cout << "Setka.cpp    " << 4 << endl;
//		if (bol == true)
//		{
//			dekard_skorost(y_ex, z_ex, x_ex, Wr[I], Wphi[I], Wthe[I], bb, cc, aa);
//			/*for (int k = 0; k < 30000; k++)
//			{
//				cout << "Setka.cpp    " << x_ex + 0.3 * k * aa << " " << sqrt(kv(y_ex + 0.3 * k * bb) + kv(z_ex + 0.3 * k * cc)) << endl;
//			}*/
//			Fly_exchenge_Split(sens, x_ex, y_ex, z_ex, aa, bb, cc, now, mu_[I] * mu_ex, mu_ex, true, I_ - 1);
//		}
//		//exit(-1);
//
//		/*mmu5 += mu_[I];
//		mn5++;
//
//		if (I >= 1)
//		{
//			mmu1 += mu_[0];
//			mn1++;
//		}
//		if (I >= 2)
//		{
//			mmu2 += mu_[1];
//			mn2++;
//		}
//		if (I >= 3)
//		{
//			mmu3 += mu_[2];
//			mn3++;
//		}
//		if (I >= 4)
//		{
//			mmu4 += mu_[3];
//			mn4++;
//		}*/
//		
//	}
//
//	/*if (mu2 < weight_ * mu_0 && ExCh == false)
//	{
//		if (mu2 >= sens->MakeRandom() * weight_ * mu_0)
//		{
//			mu2 = weight_ * mu_0;
//		}
//		else
//		{
//			return;
//		}
//	}*/
//	bool bol = true;
//	//cout << "Setka.cpp    " << 5 << endl;
//
//	if (mu2 < Mu[zone] && ExCh == false)
//	{
//		if (mu2 >= sens->MakeRandom() * Mu[zone])
//		{
//			mu2 = Mu[zone];
//		}
//		else
//		{
//			return;
//		}
//	}
//	//cout << "Setka.cpp    " << 6 << endl;
//	alpha = polar_angle(y_ex + 0.5 * t2 * Vy, z_ex + 0.5 * t2 * Vz);
//	u = sqrt(kvv(Vx - vx, Vy - vy * cos(alpha), Vz - vy * sin(alpha)));
//	uz = Velosity_1(u, cp);
//	uz_M = Velosity_2(u, cp) / (uz * kv(cp) * cp * pi * sqrtpi_);
//	uz_E = Velosity_3(u, cp);
//	u1 = vx - Vx;
//	u2 = vy * cos(alpha) - Vy;
//	u3 = vy * sin(alpha) - Vz;
//	double skalar = Vx * u1 + Vy * u2 + Vz * u3;
//
//	now->mut.lock();
//	now->par[0].F_n += t2 * mu2;
//	now->par[0].F_u += t2 * Vx * mu2;
//	now->par[0].F_v += t2 * (Vy * cos(alpha) + Vz * sin(alpha)) * mu2;
//	now->par[0].F_T += t2 * kvv(Vx, Vy, Vz) * mu2;
//	now->par[0].I_u += mu2 * t2 * uz_M * uz * sigma(uz_M) * u1 / u;
//	now->par[0].I_v += mu2 * t2 * uz_M * uz * sigma(uz_M) * (u2 * cos(alpha) + u3 * sin(alpha)) / u;
//	now->par[0].I_T += mu2 * t2 * (0.5 * (3.0 * kv(cp) + 2.0 * kv(u)) * uz_E * sigma(uz_E) + 2.0 * uz_M * uz * sigma(uz_M) * skalar / u);
//	now->mut.unlock();
//
//	double X = x_ex + t2 * Vx;
//	double Y = y_ex + t2 * Vy;
//	double Z = z_ex + t2 * Vz;
//
//	//cout << "Setka.cpp    " << 7 << endl;
//	Cell* next = nullptr;             // В какую ячейку попадаем следующую
//	double xk, yk;
//	xk = x_0 + (1.0 * b) * dt * Vx;
//	yk = sqrt(kv(y_0 + (1.0 * b) * dt * Vy) + kv(z_0 + (1.0 * b) * dt * Vz));
//	if (xk < Left_ || (sqrt(kv(xk) + kv(yk)) > R5_ && xk > 0) || (xk < 0 && yk > R5_))
//	{
//		next = nullptr;
//	}
//	else
//	{
//		for (auto& i : now->Grans)
//		{
//			if (i->Sosed != nullptr)
//			{
//				if (i->Sosed->belong(xk, yk) == true)
//				{
//					next = i->Sosed;
//				}
//			}
//		}
//		if (next == nullptr)
//		{
//			for (auto& k : now->Grans)
//			{
//				if (k->Sosed != nullptr)
//				{
//					for (auto& i : k->Sosed->Grans)
//					{
//						if (i->Sosed != nullptr)
//						{
//							if (i->Sosed->belong(xk, yk) == true)
//							{
//								next = i->Sosed;
//							}
//						}
//					}
//				}
//			}
//		}
//		/*if (next == nullptr)
//		{
//			cout << "Setka.cpp    " << xk << " " << yk << " " << "no sosed" << endl;
//		}*/
//	}
//	//cout << "Setka.cpp    " << 8 << endl;
//	if (next != nullptr)
//	{
//		//cout << "Setka.cpp    " << "Next cell" << endl;
//		Fly_exchenge_Split(sens, X, Y, Z, Vx, Vy, Vz, next, mu2, mu_0, false, zone);
//		//cout << "Setka.cpp    " << "Vozvrat next" << endl;
//	}
//	return;
//	
//}

void Setka::Change_Velosity(Sensor* s, const double& Ur, const double& Uthe, const double& Uphi, //
	const double& Vr, const double& Vthe, const double& Vphi, double& X, double& Y, double& Z, const double& cp)
{
	double x = sqrt(kvv(Vr - Ur, Vthe - Uthe, Vphi - Uphi));
	double p4 = 0.5 * sqrtpi_ * x / (1.0 + 0.5 * sqrtpi_ * x);
	double ksi1, ksi2, ksi3, ksi4, ksi5, ksi6;
	double om1, om2, om3, lo;
	double y1, y2, y3;
	double v1, v2, v3, u1, u2, u3;
	double uu, yy, vv, D, ko;
	do
	{
		ksi1 = s->MakeRandom();
		ksi2 = s->MakeRandom();
		ksi3 = s->MakeRandom();
		ksi4 = s->MakeRandom();
		ksi5 = s->MakeRandom();
		ksi6 = s->MakeRandom();
		//cout << "Setka.cpp    " << "sd " << endl;
		if (p4 < ksi1)
		{
			om1 = 1.0 - 2.0 * ksi4;
			om2 = sqrt(1.0 - kv(om1)) * cos(2.0 * pi_ * ksi5);
			om3 = sqrt(1.0 - kv(om1)) * sin(2.0 * pi_ * ksi5);
			// Более экономичный алгоритм
			/*do
			{
				om2 = 1.0 - 2.0 * s->MakeRandom();
				om3 = 1.0 - 2.0 * s->MakeRandom();
				D = kv(om2) + kv(om3);
			} while (D > 1);
			ko = sqrt((1.0 - kv(om1)) / D);
			om2 = om2 * ko;
			om3 = om3 * ko;*/

			lo = sqrt(-log(ksi2 * ksi3));
			y1 = lo * om1;
			y2 = lo * om2;
			y3 = lo * om3;
		}
		else
		{
			y1 = sqrt(-log(ksi2)) * cos(pi_ * ksi3);
			y2 = sqrt(-log(ksi4)) * cos(2.0 * pi_ * ksi5);
			y3 = sqrt(-log(ksi4)) * sin(2.0 * pi_ * ksi5);
		}
		v1 = y1 + Ur;
		v2 = y2 + Uthe;
		v3 = y3 + Uphi;
		u1 = Vr - v1;
		u2 = Vthe - v2;
		u3 = Vphi - v3;
		uu = sqrt(kvv(u1, u2, u3));
		yy = sqrt(kvv(y1, y2, y3));
	} while (((uu * sigma2(uu, cp)) / (sigma2(x, cp) * (x + yy))) <= ksi6);
	//cout << "Setka.cpp    " << v2 << endl;
	X = v1;
	Y = v2;
	Z = v3;
}

void Setka::Change_Velosity_Split(Sensor* s, const double& Ur, const double& Uthe, const double& Uphi, //
	const double& Vr, const double& Vthe, const double& Vphi, vector <double>& Wr, vector <double>& Wthe,//
	vector <double>& Wphi, vector <double>& mu_, const double& cp, const double& r, int I, const double& x_ex, const double& y_ex, const double& z_ex)
	// I -  сколько дополнительных атомов вылетает (кроме основного)
{
	// I - в какой зоне мы сейчас находимся?
	// r - координата точки перезарядки
	double e1 = sqrtpi_ * kv(Ur);
	double e2 = 2.0 * fabs(Ur);
	double e3 = 0.5 * sqrtpi_;
	double p1 = e1 / (e1 + e2 + e3);
	double p2 = e2 / (e1 + e2 + e3);
	double ksi1, ksi2, ksi3, ksi4, ksi5, ksi6, ksi7, ksi8;
	double z = 0.0, h = 0.0;
	//cout << "Setka.cpp    " << "A " << endl;
	for (int i = 0; i < I; i++)
	{
		do
		{
			ksi1 = s->MakeRandom();
			if (p1 > ksi1)
			{
				ksi2 = s->MakeRandom();
				ksi3 = s->MakeRandom();
				z = sqrt(-log(ksi2)) * cos(ksi3 * pi_);
			}
			else if (p1 + p2 > ksi1)
			{
				ksi2 = s->MakeRandom();
				if (ksi2 <= 0.5)
				{
					z = -sqrt(-log(2.0 * ksi2));
				}
				else
				{
					z = sqrt(-log(2.0 * (1.0 - ksi2)));
				}
			}
			else
			{
				ksi2 = s->MakeRandom();
				ksi3 = s->MakeRandom();
				ksi4 = s->MakeRandom();
				ksi5 = s->MakeRandom();
				z = sign(ksi5 - 0.5) * sqrt(-log(ksi2) - log(ksi3) * //
					kv(cos(pi_ * ksi4)));
			}

			Wr[i] = Ur + z;
			h = kv(Ur + z) / kv(fabs(Ur) + fabs(z));
			ksi6 = s->MakeRandom();
		} while (h <= ksi6 || z >= -Ur);
	}
	//cout << "Setka.cpp    " << "B " << endl;
	//vector <double> w_(I + 1);
	vector <double> gam_(I + 1);
	vector <double> F_(I + 1);
	for (int i = 0; i < I; i++)
	{
		gam_[i] = 1.0 / (kv(r) / kv(Ri[i]) - 1.0);
		//cout << "Setka.cpp    " << gam_[i] << endl;
		//w_[i] = sqrt(kv(Wr[i]) * gam_[i]);
		F_[i] = 0.5 * gam_[i] * ((0.5 + kv(Ur)) *//
			(1.0 - sign(Ur) * erf(fabs(Ur))) - //
			Ur * exp(-kv(Ur)) / sqrtpi_);
	}
	//cout << "Setka.cpp    " << "C " << endl;
	double Val;
	double The;

	for (int i = 0; i < I; i++)
	{
		ksi7 = s->MakeRandom();
		ksi8 = s->MakeRandom();
		The = 2.0 * pi_ * ksi7;
		if (i > 0)
		{
			Val = sqrt(kv(w_c_v_s(r, Wr[i], i - 1)) + ksi8 * (kv(w_c_v_s(r, Wr[i], i)) - kv(w_c_v_s(r, Wr[i], i - 1))));
		}
		else
		{
			Val = sqrt(ksi8 * (kv(w_c_v_s(r, Wr[i], i))));
		}
		Wthe[i] = Val * cos(The);
		Wphi[i] = Val * sin(The);
	}
	//cout << "Setka.cpp    " << "D " << endl;
	double x = sqrt(kvv(Vr - Ur, Vthe - Uthe, Vphi - Uphi));
	double uu = exp(-kv(x)) / sqrtpi_ + (x + 1.0 / (2.0 * x)) * erf(x);

	for (int i = 0; i < I; i++)
	{
		double u = sqrt(kvv(Vr - Wr[i], Vthe - Wthe[i], Vphi - Wphi[i]));
		if (i > 0)
		{
			mu_[i] = (F_[i] - F_[i - 1]) * (u * sigma2(u, cp) / (uu * sigma2(uu, cp))) * exp(-kv(Wthe[i] - Uthe) - kv(Wphi[i]));
		}
		else
		{
			mu_[i] = F_[i] * (u * sigma2(u, cp) / (uu * sigma2(uu, cp))) * exp(-kv(Wthe[i] - Uthe) - kv(Wphi[i]));
		}
	}

	// Расчёт основного атома
	//x = sqrt(kvv(Vr - Ur, Vthe - Uthe, Vphi - Uphi));
	double p4 = 0.5 * sqrtpi_ * x / (1.0 + 0.5 * sqrtpi_ * x);
	double om1, om2, om3, lo;
	double y1, y2, y3;
	double v1, v2, v3, u1, u2, u3;
	double yy, vv, D, ko, uuu;
	//cout << "Setka.cpp    " << "E " << endl;
	do
	{
		ksi1 = s->MakeRandom();
		ksi2 = s->MakeRandom();
		ksi3 = s->MakeRandom();
		ksi4 = s->MakeRandom();
		ksi5 = s->MakeRandom();
		ksi6 = s->MakeRandom();
		//cout << "Setka.cpp    " << "sd " << endl;
		if (p4 < ksi1)
		{
			om1 = 1.0 - 2.0 * ksi4;
			om2 = sqrt(1.0 - kv(om1)) * cos(2.0 * pi_ * ksi5);
			om3 = sqrt(1.0 - kv(om1)) * sin(2.0 * pi_ * ksi5);
			// Более экономичный алгоритм
			/*do
			{
				om2 = 1.0 - 2.0 * s->MakeRandom();
				om3 = 1.0 - 2.0 * s->MakeRandom();
				D = kv(om2) + kv(om3);
			} while (D > 1);
			ko = sqrt((1.0 - kv(om1)) / D);
			om2 = om2 * ko;
			om3 = om3 * ko;*/

			lo = sqrt(-log(ksi2 * ksi3));
			y1 = lo * om1;
			y2 = lo * om2;
			y3 = lo * om3;
		}
		else
		{
			y1 = sqrt(-log(ksi2)) * cos(pi_ * ksi3);
			y2 = sqrt(-log(ksi4)) * cos(2.0 * pi_ * ksi5);
			y3 = sqrt(-log(ksi4)) * sin(2.0 * pi_ * ksi5);
		}
		v1 = y1 + Ur;
		v2 = y2 + Uthe;
		v3 = y3 + Uphi;
		u1 = Vr - v1;
		u2 = Vthe - v2;
		u3 = Vphi - v3;
		uuu = sqrt(kvv(u1, u2, u3));
		yy = sqrt(kvv(y1, y2, y3));
		h = ((uuu * sigma2(uuu, cp)) / (sigma2(x, cp) * (x + yy)));
	} while (h <= ksi6 || (v1 <= 0 && kvv(v2, v3, 0.0) < kv(w_c_v_s(r, v1, I - 1))));
	//cout << "Setka.cpp    " << v2 << endl;
	Wr[I] = v1;
	Wthe[I] = v2;
	Wphi[I] = v3;
	mu_[I] = 1.0;
	for (int i = 0; i < I; i++)
	{
		mu_[I] = mu_[I] - mu_[i];
	}

	/*if (mu_[I] < -5.0)
	{
		double aa, bb, cc;
		cout << "Setka.cpp    " << "8742  MINUS  I = " << mu_[I] << "  "  << I << "  r = " << r << endl;
		cout << "Setka.cpp    " << Ur << " " << Uthe << " " << Uphi << " " << Vr << " " << Vthe << " " << Vphi << endl;
		cout << "Setka.cpp    " << "x = " << x << "   uu = " << uu << endl;
		cout << "Setka.cpp    " << "erf(x) = " << erf(x) << endl;
		double u = sqrt(kvv(Vr - Wr[3], Vthe - Wthe[3], Vphi - Wphi[3]));
		cout << "Setka.cpp    " << " u = " << u << endl;
		cout << "Setka.cpp    " << (F_[3] - F_[2]) * (u * sigma2(u, cp) / (uu * sigma2(uu, cp))) * exp(-kv(Wthe[3] - Uthe) - kv(Wphi[3])) << endl;
		cout << "Setka.cpp    " << (u * sigma2(u, cp) / (uu * sigma2(uu, cp))) << endl;
		exit(-1);
		for (int i = 0; i < I; i++)
		{
			cout << "Setka.cpp    " << "mu = " << mu_[i] << endl;
		}
		for (int i = 0; i < I; i++)
		{
			cout << "Setka.cpp    " << "gam = " << gam_[i] << endl;
			cout << "Setka.cpp    " << "F = " << F_[i] << endl;
		}
		for (int i = 0; i <= I; i++)
		{
			cout << "Setka.cpp    " << "Wr = " << Wr[i] << "    Wthe = " << Wthe[i] << "    Wphi = " << Wphi[i] << endl;
		}
		for (int i = 0; i < I; i++)
		{
			dekard_skorost(y_ex, z_ex, x_ex, Wr[i], Wphi[i], Wthe[i], bb, cc, aa);
			for (int ij = 0; ij < 27000; ij++)
			{
				cout << "Setka.cpp    " << x_ex + ij * aa << " " << sqrt(kvv(y_ex + ij * bb, z_ex + ij * cc, 0.0)) << " " << i << endl;
			}
		}
		dekard_skorost(y_ex, z_ex, x_ex, Wr[I], Wphi[I], Wthe[I], bb, cc, aa);
		for (int ij = 0; ij < 30000; ij++)
		{
			cout << "Setka.cpp    " << x_ex + ij * aa << " " << sqrt(kvv(y_ex + ij * bb, z_ex + ij * cc, 0.0)) << " " << I << endl;
		}
		exit(-1);
	}*/
}

double Setka::w_c_v_s(const double& r, const double& v, int i)  // W для перезарядки
{
	double gamma = 1.0 / (kv(r) / kv(Ri[i]) - 1.0);
	return sqrt(gamma * kv(v));
}

void Setka::func_pogloshenie(void)
{
	// Считаем поглощение вдоль луча

	for (int ijk = 0; ijk <= 18; ijk++)
	{
		for (int sort = 0; sort < 4; sort++)
		{
			double alf = ijk * (pi_ / 2.0) / (9.0);
			double e1 = cos(alf);
			double e2 = sin(alf);
			double dl = pogRmax / 300;
			double x = 0.0, y = 0.0;
			double du = (pogVmax - pogVmin) / pogl_rad_;
			double u = 0.0;

			int b;

			double n = 1.0, cp = 1.0, ux = 0.0, uy = 0.0;

			double pogl[pogl_rad_];
			double pogl2[pogl_rad_];
			for (int i = 0; i < pogl_rad_; i++)
			{
				pogl[i] = 0.0;
				pogl2[i] = 0.0;
			}


			for (int i = 0; i < 300; i++)
			{
				x = e1 * (i * dl + dl / 2.0);
				y = e2 * (i * dl + dl / 2.0);

				auto C = Find_cell(b, x, y);
				if (b == 0)
				{
					cout << "ERRE 22423424343324  " << x << " " << y << endl;
					continue;
				}

				n = C->par[0].H_n[sort];
				cp = sqrt(C->par[0].H_T[sort]);
				ux = C->par[0].H_u[sort];
				uy = C->par[0].H_v[sort];

				if (cp <= 0.00000001) cp = 1.0;

				//cout << n << " " << ux << " " << uy << " " << cp << endl;

				for (int j = 0; j < pogl_rad_; j++)
				{
					u = pogVmin + (j + 0.5) * du;

					pogl[j] += dl * exp(-kv(u - ux * e1 - uy * e2) / kv(cp)) * 3.0 * n / (sqrtpi_ * cp);
					pogl2[j] += 3.0 * dl * C->pogloshenie[sort][j] / du;
				}
			}

			ofstream fout_cr;
			fout_cr.open(to_string(sort) + "_Pogloshenie_maxwell_" + to_string(ijk) + ".txt");

			for (int j = 0; j < pogl_rad_; j++)
			{
				u = pogVmin + (j + 0.5) * du;
				fout_cr << u << " " << pogl[j] << " " << pogl2[j] << endl;
				// << exp(-pogl[j] * 45.2717) << endl;
			}

			fout_cr.close();
		}
	}

}

double Setka::Velosity_1(const double& u, const double& cp)
{
	if (u < 0.00001)
	{
		return 2.0 * cp / sqrtpi_ + 2.0 * u * u / (3.0 * cp * sqrtpi_) - u * u * u * u / (15.0 * cp * cp * cp * sqrtpi_);
	}
	else
	{
		return  exp(-u * u / kv(cp)) * cp / sqrtpi_ + (u + kv(cp) / (2.0 * u)) * erf(u / cp);
	}
}

double Setka::Velosity_2(const double& u, const double& cp)  // Считает на совсем скорость, а только её числитель (см. статью)
{
	if (u < 0.00001)
	{
		return (8.0 / 3.0) * kv(cp) * kv(cp) * pi * u + (8.0 / 15.0) * kv(cp) * pi * u * u * u - (4.0 / 105.0) * pi * kv(u) * kv(u) * u;
	}
	else
	{
		return  cp * cp * cp * pi * (exp(-u * u / kv(cp)) * cp * u * 2.0 * (kv(cp) + 2.0 * kv(u)) +//
			sqrtpi_ * (4.0 * kv(u) * kv(u) + 4.0 * cp * cp * kv(u) - kv(cp) * kv(cp)) * erf(u / cp)) / (4.0 * u * u);
	}
}

double Setka::Velosity_3(const double& u, const double& cp)
{
	if (u < 0.00001)
	{
		return 8.0 * cp / (3.0 * sqrtpi_) + 8.0 * u * u / (9.0 * cp * sqrtpi_) - 44.0 * u * u * u * u / (135.0 * cp * cp * cp * sqrtpi_);
	}
	else
	{
		return  exp(-u * u / kv(cp)) * cp * (5.0 * kv(cp) + 2.0 * kv(u)) / (sqrtpi_ * (3.0 * kv(cp) + 2.0 * kv(u))) +//
			(4.0 * kv(u) * kv(u) + 12.0 * cp * cp * kv(u) + 3.0 * kv(cp) * kv(cp)) * erf(u / cp) / (2.0 * u * (3.0 * kv(cp) + 2.0 * kv(u)));
	}
}

void Setka::test_velosity(const double& a, const double& b)
{
	double c;
	c = sin(a) + cos(b) * sin(log(a * a + b * a + b * b)) + minmod(a,b);
}

void Setka::test_main(void)
{
	for (unsigned long int n = 0; n < 1000; ++n) {
		for (unsigned long int n2 = 0; n2 < 1000; ++n2) {
			for (unsigned long int n3 = 0; n3 < 1000; ++n3) {
				double a = 1.23423432;
				double b = 2.23424123;
				this->test_velosity(a, b);
			}
		}
	}
}

void Setka::Smooth_kvadr3(const double& x1, const double& y1, const double& z1, //
	const double& x2, const double& y2, const double& z2, const double& x3, const double& y3,//
	const double& z3, const double& x4, const double& y4, const double& z4, const double& x5,//
	const double& y5, const double& z5, double& xx, double& yy, double& zz)
{
	//Функция определяющая новвое положение точки x, y, z
		   //В соответствии с тремя другими точками, строя квадратичный сплайн по методу наименьших квадратов

	double ex[3], ey[3], c[3], xx2, yy2, xx3, yy3, xx4, yy4, xx5, nn;
	double a, b, cc, g2, g3, g4, h2, h4;


	ex[0] = x5 - x1;
	ex[1] = y5 - y1;
	ex[2] = z5 - z1;

	ey[0] = x3 - x1;
	ey[1] = y3 - y1;
	ey[2] = z3 - z1;

	nn = sqrt(kv(ex[0]) + kv(ex[1]) + kv(ex[2]));
	ex[0] /= nn;
	ex[1] /= nn;
	ex[2] /= nn;

	nn = sqrt(kv(ey[0]) + kv(ey[1]) + kv(ey[2]));
	ey[0] /= nn;
	ey[1] /= nn;
	ey[2] /= nn;


	//Если они не сонаправлены
	if (fabs(DOT_PRODUCT(ex, ey)) < 0.99)
	{
		c[0] = ex[1] * ey[2] - ex[2] * ey[1];
		c[1] = ex[2] * ey[0] - ex[0] * ey[2];
		c[2] = ex[0] * ey[1] - ex[1] * ey[0];

		ey[0] = c[1] * ex[2] - c[2] * ex[1];
		ey[1] = c[2] * ex[0] - c[0] * ex[2];
		ey[2] = c[0] * ex[1] - c[1] * ex[0];


		nn = sqrt( kv(ey[0]) + kv(ey[1]) + kv(ey[2]));
		ey[0] /= nn;
		ey[1] /= nn;
		ey[2] /= nn;


		c[0] = x2 - x1; c[1] = y2 - y1; c[2] = z2 - z1;
		xx2 = DOT_PRODUCT(ex, c);
		yy2 = DOT_PRODUCT(ey, c);

		c[0] = x3 - x1; c[1] = y3 - y1; c[2] = z3 - z1;
		xx3 = DOT_PRODUCT(ex, c);
		yy3 = DOT_PRODUCT(ey, c);

		c[0] = x4 - x1; c[1] = y4 - y1; c[2] = z4 - z1;
		xx4 = DOT_PRODUCT(ex, c);
		yy4 = DOT_PRODUCT(ey, c);

		c[0] = x5 - x1; c[1] = y5 - y1; c[2] = z5 - z1;
		xx5 = DOT_PRODUCT(ex, c);

		g2 = xx2 / xx5;
		g3 = xx3 / xx5;
		g4 = xx4 / xx5;

		if (fabs(y3) < 0.0000001) 
		{
			goto a11;
		}

		h2 = yy2 / yy3;
		h4 = yy4 / yy3;


		a = ((-1 + g4) * g4 * ((-1 + g4) * (1 + h2) + 3 * h4) + kv3(g3) * (h2 + g4 * h2 + h4 - 3 * g4 * h4) +
			kv3(g2) * (1 + g4 + g3 * (-3 + h4) + h4 - 3 * g4 * h4) + g3 * (-3 + h2 + h4 + g4 * (1 + g4 *
				(1 + g4 * (-3 + h2) - 2 * h4) + h4)) + g2 * (1 - 3 * h2 +
					g3 * (1 + kv(g4)) * (1 + h2) + h4 - 2 * g3 * kv(g4) * h4 + kv3(g3) * (-3 * h2 + h4) +
					g4 * (h2 + g4 * (g4 + h2 - 3 * g4 * h2 - 2 * h4) + h4) + kv(g3) *
					(-2 + h2 + g4 * (-2 + h2 + h4))) + kv(g3) * (3 - 2 * h2 - 2 * h4 +
						g4 * (-2 + h4 + g4 * (3 - 2 * h2 + 3 * h4))) + kv(g2) * (-2 + 3 * h2 +
							kv(g3) * (3 + 3 * h2 - 2 * h4) - 2 * h4 + g3 * (1 + g4 - 2 * h2 - 2 * g4 * h2 + g4 * h4) + g4 * (-2 * h2 + h4 +
								g4 * (-2 + 3 * h2 + 3 * h4)))) / (3 * kv(-1 + g4) * kv(g4) -
									2 * g3 * kv(-1 + g4) * g4 * (1 + g4) + 2 * kv3(g3) * (-3 + g4 + kv(g4) - 3 * kv3(g4)) +
									kv4(g3) * (3 + g4 * (-2 + 3 * g4)) + kv4(g2) * (3 + 3 * kv(g3) -
										2 * g3 * (1 + g4) + g4 * (-2 + 3 * g4)) + 2 * kv3(g2) * (-3 - 3 * kv3(g3) + g4 + kv(g4) -
											3 * kv3(g4) + kv(g3) * (1 + g4) + g3 * (1 + kv(g4))) + 2 * g2 * (-(kv4(g3) * (1 + g4)) -
												kv(-1 + g4) * g4 * (1 + g4) + kv3(g3) * (1 + kv(g4)) + kv(g3) * (1 + kv3(g4)) -
												g3 * (1 + kv4(g4))) + kv(g3) * (3 + g4 * (2 +
													g4 * (-6 + g4 * (2 + 3 * g4)))) + kv(g2) * (3 + 3 * kv4(g3) + 2 * kv3(g3) * (1 + g4) - 6 * kv(g3) * (1 + kv(g4)) +
														2 * g3 * (1 + kv3(g4)) + g4 * (2 + g4 * (-6 + g4 * (2 + 3 * g4)))));

		b = ((-1 + g4) * g4 * (1 + h2 - kv(g4) * (1 + h2) - 3 * h4) - g3 * (-3 + kv4(g4) * (-3 + h2) +
			h2 - kv(g4) * (-2 + h4) + h4) - kv4(g3) * (h2 + g4 * h2 + h4 - 3 * g4 * h4) - kv4(g2) * (1 + g4 + g3 * (-3 + h4) +
				h4 - 3 * g4 * h4) + kv3(g3) * (h2 + kv(g4) * h2 + h4 - 3 * kv(g4) * h4) + kv3(g2) * (1 + kv(g4) +
					kv(g3) * (-3 + h4) + h4 - 3 * kv(g4) * h4) + kv(g3) * (-3 + h2 + h4 + g4 * (1 - 2 * h4 + g4 * (1 + g4 * (-3 + h2) + h4))) + g2 *
			(-1 + 3 * h2 + kv4(g3) * (3 * h2 - h4) - h4 + kv(g4) * (-2 * h2 + kv(g4) * (-1 + 3 * h2) + h4) + kv(g3) * (1 - 2 * h2 +
				kv(g4) * (1 - 2 * h2 + h4))) + kv(g2) * (1 + g3 * (1 + kv(g4)) * (-2 + h2) - 3 * h2 + h4 + g3 * kv(g4) * h4 + kv3(g3) * (-3 * h2 + h4) +
					kv(g3) * (1 + g4 + h2 + g4 * h2 - 2 * g4 * h4) + g4 * (h2 - 2 * h4 + g4 * (g4 + h2 - 3 * g4 * h2 + h4)))) /
			(3 * kv(-1 + g4) * kv(g4) - 2 * g3 * kv(-1 + g4) * g4 * (1 + g4) + 2 * kv3(g3) * (-3 + g4 + kv(g4) - 3 * kv3(g4)) +
				kv4(g3) * (3 + g4 * (-2 + 3 * g4)) + kv4(g2) * (3 + 3 * kv(g3) - 2 * g3 * (1 + g4) + g4 * (-2 + 3 * g4)) + 2 * kv3(g2) *
				(-3 - 3 * kv3(g3) + g4 + kv(g4) - 3 * kv3(g4) + kv(g3) * (1 + g4) + g3 * (1 + kv(g4))) + 2 * g2 * (-(kv4(g3) * (1 + g4)) -
					kv(-1 + g4) * g4 * (1 + g4) + kv3(g3) * (1 + kv(g4)) + kv(g3) * (1 + kv3(g4)) - g3 * (1 + kv4(g4))) +
				kv(g3) * (3 + g4 * (2 + g4 * (-6 + g4 * (2 + 3 * g4)))) + kv(g2) * (3 + 3 * kv4(g3) + 2 * kv3(g3) * (1 + g4) -
					6 * kv(g3) * (1 + kv(g4)) + 2 * g3 * (1 + kv3(g4)) + g4 * (2 + g4 * (-6 + g4 * (2 + 3 * g4)))));

		cc = (kv(-1 + g4) * kv(g4) * (1 + h2) - g3 * (-1 + g4) * g4 * (-1 + kv(g4) - h4) + kv4(g3) * (h2 + kv(g4) * h2 + h4 - g4 * h4) +
			kv3(g3) * (-2 * (1 + kv3(g4)) * h2 + (-2 + g4 + kv(g4)) * h4) + kv(g3) * (h2 + h4 + g4 * (1 + g4 * (-2 + g4 + kv(g4) * h2 -
				2 * h4) + h4)) + kv4(g2) * (1 + kv(g4) + h4 + kv(g3) * h4 - g4 * h4 - g3 * (1 + g4 + g4 * h4)) + g2 * (kv(g3) * (1 + kv3(g4)) * (1 + h2) -
					g3 * (1 + kv4(g4)) * (1 + h2) + (-1 + g4) * g4 * (h2 - kv(g4) * h2 + h4) - kv4(g3) * (h2 + g4 * h2 + g4 * h4) + kv3(g3) *
					(h2 + kv(g4) * h2 + kv(g4) * h4)) + kv3(g2) * (-2 * (1 + kv3(g4)) - 2 * kv3(g3) * h4 + (-2 + g4 + kv(g4)) * h4 + kv(g3) * (1 + g4 + g4 * h4) +
						g3 * (1 + kv(g4) * (1 + h4))) + kv(g2) * (1 + g3 * (1 + kv3(g4)) * (1 + h2) + h4 + kv4(g3) * h4 + kv3(g3) * (h2 + g4 * h2 + g4 * h4) -
							2 * kv(g3) * (1 + h2 + kv(g4) * (1 + h2 + h4)) + g4 * (h2 + h4 + g4 * (g4 * (g4 + h2) - 2 * (h2 + h4))))) /
			(3 * kv(-1 + g4) * kv(g4) - 2 * g3 * kv(-1 + g4) * g4 * (1 + g4) + 2 * kv3(g3) * (-3 + g4 + kv(g4) - 3 * kv3(g4)) +
				kv4(g3) * (3 + g4 * (-2 + 3 * g4)) + kv4(g2) * (3 + 3 * kv(g3) - 2 * g3 * (1 + g4) + g4 * (-2 + 3 * g4)) + 2 * kv3(g2) * (-3 - 3 * kv3(g3) + g4 + kv(g4) -
					3 * kv3(g4) + kv(g3) * (1 + g4) + g3 * (1 + kv(g4))) + 2 * g2 * (-(kv4(g3) * (1 + g4)) - kv(-1 + g4) * g4 * (1 + g4) +
						kv3(g3) * (1 + kv(g4)) + kv(g3) * (1 + kv3(g4)) - g3 * (1 + kv4(g4))) + kv(g3) * (3 + g4 * (2 +
							g4 * (-6 + g4 * (2 + 3 * g4)))) + kv(g2) * (3 + 3 * kv4(g3) + 2 * kv3(g3) * (1 + g4) - 6 * kv(g3) * (1 + kv(g4)) +
								2 * g3 * (1 + kv3(g4)) + g4 * (2 + g4 * (-6 + g4 * (2 + 3 * g4)))));


			yy3 = (a * kv(xx3 / xx5) + b * (xx3 / xx5) + cc) * yy3;


		xx = x1 + xx3 * ex[1] + yy3 * ey[1];
		yy = y1 + xx3 * ex[2] + yy3 * ey[2];
		zz = z1 + xx3 * ex[3] + yy3 * ey[3];

		return;
		}
		else
		{
			xx = x3;
			yy = y3;
			zz = z3;

			return;
		}

				a11:
				xx = x3;
				yy = y3;
				zz = z3;

				return;

}
