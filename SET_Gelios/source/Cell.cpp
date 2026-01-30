#include "Cell.h"
#include <math.h>


Cell::Cell(Point* A, Point* B, Point* C, Point* D)
{
	this->contour.push_back(A);
	this->contour.push_back(B);
	this->contour.push_back(C);
	this->contour.push_back(D);
	this->number = -1;
	this->type = C_no;
	this->Initial();
}

Cell::Cell(Point* A, Point* B, Point* C)
{
	this->contour.push_back(A);
	this->contour.push_back(B);
	this->contour.push_back(C);
	this->number = -1;
	this->type = C_no;
	this->Initial();
}

Cell::Cell(void)
{
	this->number = -1;
	this->type = C_no;
	this->Initial();
}

void Cell::Initial(void)
{
	this->par[0].num_atoms[0] = 0;
	this->par[0].num_atoms[1] = 0;
	this->par[0].num_atoms[2] = 0;
	this->par[0].num_atoms[3] = 0;
	this->axis_ = false;
	this->y_ax = 0.0;

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < pogl_rad_; j++)
		{
			pogloshenie[i][j] = 0.0;
		}
	}
	

}

void Cell::Get_Center2(double& x, double& y)
// Среднее арифметическое вершин
{
	if (this->contour.size() == 4)
	{
		x = 0.0;
		y = 0.0;
		for (auto& i : this->contour)
		{
			x += i->x;
			y += i->y;
		}
		x = x / (1.0 * this->contour.size());
		y = y / (1.0 * this->contour.size());
		/*if (this->contour[0]->y < 0.0001)
		{
			y = y * (4.0 / 3.0);
		}*/
	}
	else
	{
		x = 0.0;
		y = 0.0;
		for (auto& i : this->contour)
		{
			x += i->x;
			y += i->y;
		}
		x = x / (1.1 * this->contour.size());
		y = y / (1.1 * this->contour.size());
	}
}

void Cell::Get_Center(double& x, double& y)
// Центр масс - барицентр многоугольника
{
	double A = 0.0;
	double Gx = 0.0;
	double Gy = 0.0;
	int k = 0;
	for (int i = 0; i < this->contour.size(); i++)
	{
		auto a = this->contour[i];
		k = i + 1;
		if (k == this->contour.size()) k = 0;
		auto b = this->contour[k];
		A = A + (a->x * b->y - a->y * b->x);
		Gx = Gx + (a->x + b->x) * (a->x * b->y - a->y * b->x);
		Gy = Gy + (a->y + b->y) * (a->x * b->y - a->y * b->x);
	}
	A = A * 0.5;
	Gx = Gx / (6.0 * A);
	Gy = Gy / (6.0 * A);
	x = Gx;
	y = Gy;
}

void Cell::Get_Center_posle(double& x, double& y)
{
	double A = 0.0;
	double Gx = 0.0;
	double Gy = 0.0;
	int k = 0;
	for (int i = 0; i < this->contour.size(); i++)
	{
		auto a = this->contour[i];
		k = i + 1;
		if (k == this->contour.size()) k = 0;
		auto b = this->contour[k];
		A = A + (a->x2 * b->y2 - a->y2 * b->x2);
		Gx = Gx + (a->x2 + b->x2) * (a->x2 * b->y2 - a->y2 * b->x2);
		Gy = Gy + (a->y2 + b->y2) * (a->x2 * b->y2 - a->y2 * b->x2);
	}
	A = A * 0.5;
	Gx = Gx / (6.0 * A);
	Gy = Gy / (6.0 * A);
	x = Gx;
	y = Gy;
}

void Cell::Get_Center_posle2(double& x, double& y)
{
	if (this->contour.size() == 4)
	{
		x = 0.0;
		y = 0.0;
		for (auto& i : this->contour)
		{
			x += i->x2;
			y += i->y2;
		}
		x = x / (1.0 * this->contour.size());
		y = y / (1.0 * this->contour.size());

		/*if (this->contour[0]->y2 < 0.0001)
		{
			y = y * (4.0 / 3.0);
		}*/
	}
	else
	{
		x = 0.0;
		y = 0.0;
		for (auto& i : this->contour)
		{
			x += i->x2;
			y += i->y2;
		}
		x = x / (1.35 * this->contour.size());
		y = y / (1.35 * this->contour.size());
	}
}

double Cell::Get_Volume(void)
{
	if (this->contour.size() == 4)
	{
		double a1 = this->contour[0]->x - this->contour[2]->x;
		double a2 = this->contour[0]->y - this->contour[2]->y;
		double a = sqrt(kv(a1) + kv(a2));

		double b1 = this->contour[1]->x - this->contour[3]->x;
		double b2 = this->contour[1]->y - this->contour[3]->y;
		double b = sqrt(kv(b1) + kv(b2));

		a1 = a1 / a;
		a2 = a2 / a;
		b1 = b1 / b;
		b2 = b2 / b;
		double c = a1 * b1 + a2 * b2;
		double s = sqrt(1 - kv(c));
		return 0.5 * a * b * s;
	}
	else
	{
		double x = 0.0;
		double y = 0.0;
		
		this->Get_Center(x, y);

		double a1, a2, a3, p;
		double S = 0.0;
		for (auto& i : this->Grans)
		{
			a1 = sqrt(kv(i->A->x - x) + kv(i->A->y - y));
			a2 = sqrt(kv(i->B->x - x) + kv(i->B->y - y));
			a3 = sqrt(kv(i->A->x - i->B->x) + kv(i->A->y - i->B->y));
			p = (a1 + a2 + a3) / 2.0;
			S = S +  sqrt(p * (p - a1) * (p - a2) * (p - a3));
		}
		return S;
	}
	return 0.0;
}

double Cell::Get_Volume_rotate(const double& angle)
{
	double A = 0.0;
	double G = 0.0;
	int k = 0;
	for(int i = 0; i < this->contour.size(); i++)
	{
		auto a = this->contour[i];
		k = i + 1;
		if (k == this->contour.size()) k = 0;
		auto b = this->contour[k];
		A = A + (a->x * b->y - a->y * b->x);
		G = G + (a->y + b->y) * (a->x * b->y - a->y * b->x);
	}
	A = A * 0.5;

	return fabs(A) * (1.0 / (6.0 * A)) * G * pi_ * angle / 180.0;   // В градусах можно подавать угол
	
}

double Cell::Get_Volume_posle(void)
{
	if (this->contour.size() == 4)
	{
		double a1 = this->contour[0]->x2 - this->contour[2]->x2;
		double a2 = this->contour[0]->y2 - this->contour[2]->y2;
		double a = sqrt(kv(a1) + kv(a2));

		double b1 = this->contour[1]->x2 - this->contour[3]->x2;
		double b2 = this->contour[1]->y2 - this->contour[3]->y2;
		double b = sqrt(kv(b1) + kv(b2));

		a1 = a1 / a;
		a2 = a2 / a;
		b1 = b1 / b;
		b2 = b2 / b;
		double c = a1 * b1 + a2 * b2;
		double s = sqrt(1 - kv(c));
		return 0.5 * a * b * s;
	}
	else
	{
		double x = 0.0;
		double y = 0.0;

		for (auto& i : this->contour)
		{
			x += i->x2;
			y += i->y2;
		}
		x = x / (1.0 * this->contour.size());
		y = y / (1.0 * this->contour.size());

		double a1, a2, a3, p;
		double S = 0.0;
		for (auto& i : this->Grans)
		{
			a1 = sqrt(kv(i->A->x2 - x) + kv(i->A->y2 - y));
			a2 = sqrt(kv(i->B->x2 - x) + kv(i->B->y2 - y));
			a3 = sqrt(kv(i->A->x2 - i->B->x2) + kv(i->A->y2 - i->B->y2));
			p = (a1 + a2 + a3) / 2.0;
			S = S + sqrt(p * (p - a1) * (p - a2) * (p - a3));
		}
		return S;
	}
	return 0.0;
}

double Cell::Get_Volume_posle_rotate(const double& angle)
{
	double A = 0.0;
	double G = 0.0;
	int k = 0;
	for (int i = 0; i < this->contour.size(); i++)
	{
		auto a = this->contour[i];
		k = i + 1;
		if (k == this->contour.size()) k = 0;
		auto b = this->contour[k];
		A = A + (a->x2 * b->y2 - a->y2 * b->x2);
		G = G + (a->y2 + b->y2) * (a->x2 * b->y2 - a->y2 * b->x2);
	}

	A = A * 0.5;

	return fabs(A) * (1.0 / (6.0 * A)) * G * pi_ * angle / 180.0;

}



bool Cell::belong(const double& x, const double& y)
{
	// Попробуем быструю проверку
	if (x < this->x_min)
	{
		return false;
	}
	if (x > this->x_max)
	{
		return false;
	}
	if (y < this->y_min - 0.001/ RR_)
	{
		return false;
	}
	if (y > this->y_max)
	{
		return false;
	}
	// Если это не сработало, то нужно проверять точнее
	
	
	for (auto& i : this->Grans)
	{
		if (i->belong(x, y) == false)
		{
			return false;
		}
	}
	return true;
}

void Cell::renew(void)
{
	double d = 10000000000000000.0;
	for (auto& i : this->Grans)
	{
		d = min(d, i->Get_lenght());
	}
	d = min(d, sqrt(kv(this->contour[0]->x - this->contour[2]->x) + kv(this->contour[0]->y - this->contour[2]->y)));
	this->L = d;

	// Теперь найдём "коробку", внутри которой лежит ячейка
	this->x_min = this->contour[0]->x;
	this->x_max = this->contour[0]->x;
	this->y_min = this->contour[0]->y;
	this->y_max = this->contour[0]->y;

	for (auto& i : this->contour)
	{
		if (x_min > i->x)
		{
			x_min = i->x;
		}
		if (x_max < i->x)
		{
			x_max = i->x;
		}
		if (y_min > i->y)
		{
			y_min = i->y;
		}
		if (y_max < i->y)
		{
			y_max = i->y;
		}
	}
}

void Cell::Calc_Sourse(void)
{
	double U_M_H1, U_M_H2, U_M_H3, U_M_H4;
	double U_H1, U_H2, U_H3, U_H4;
	double sigma_H1, sigma_H2, sigma_H3, sigma_H4;
	double nu_H1, nu_H2, nu_H3, nu_H4;
	double q2_1, q2_2, q3;
	double u, v, ro, p, Q;

	double u_H1 = this->par[0].H_u[0], v_H1 = this->par[0].H_v[0], ro_H1 = this->par[0].H_n[0], p_H1 = 0.5 * this->par[0].H_T[0] * this->par[0].H_n[0];
	double u_H2 = this->par[0].H_u[1], v_H2 = this->par[0].H_v[1], ro_H2 = this->par[0].H_n[1], p_H2 = 0.5 * this->par[0].H_T[1] * this->par[0].H_n[1];
	double u_H3 = this->par[0].H_u[2], v_H3 = this->par[0].H_v[2], ro_H3 = this->par[0].H_n[2], p_H3 = 0.5 * this->par[0].H_T[2] * this->par[0].H_n[2];
	double u_H4 = this->par[0].H_u[3], v_H4 = this->par[0].H_v[3], ro_H4 = this->par[0].H_n[3], p_H4 = 0.5 * this->par[0].H_T[3] * this->par[0].H_n[3];


	u = this->par[0].u;
	v = this->par[0].v;
	ro = this->par[0].ro;
	p = this->par[0].p;
	Q = this->par[0].Q;
	

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

	this->par[0].I1_mc[0] = (n_p_LISM_ / Kn_) * (nu_H1 * (u_H1 - u));
	this->par[0].I1_mc[1] = (n_p_LISM_ / Kn_) * (nu_H2 * (u_H2 - u));
	this->par[0].I1_mc[2] = (n_p_LISM_ / Kn_) * (nu_H3 * (u_H3 - u));
	this->par[0].I1_mc[3] = (n_p_LISM_ / Kn_) * (nu_H4 * (u_H4 - u));

	this->par[0].I2_mc[0] = (n_p_LISM_ / Kn_) * (nu_H1 * (v_H1 - v));
	this->par[0].I2_mc[1] = (n_p_LISM_ / Kn_) * (nu_H2 * (v_H2 - v));
	this->par[0].I2_mc[2] = (n_p_LISM_ / Kn_) * (nu_H3 * (v_H3 - v));
	this->par[0].I2_mc[3] = (n_p_LISM_ / Kn_) * (nu_H4 * (v_H4 - v));

	
	this->par[0].I3_mc[0] = (n_p_LISM_ / Kn_) * (nu_H1 * ((kv(u_H1) + kv(v_H1) - kv(u) - kv(v)) / 2.0 + //
		(U_H1 / U_M_H1) * (2.0 * p_H1 / ro_H1 - p / ro)));
	this->par[0].I3_mc[1] = (n_p_LISM_ / Kn_) * (nu_H2 * ((kv(u_H2) + kv(v_H2) - kv(u) - kv(v)) / 2.0 + //
			(U_H2 / U_M_H2) * (2.0 * p_H2 / ro_H2 - p / ro)));
	this->par[0].I3_mc[2] = (n_p_LISM_ / Kn_) * (nu_H3 * ((kv(u_H3) + kv(v_H3) - kv(u) - kv(v)) / 2.0 + //
			(U_H3 / U_M_H3) * (2.0 * p_H3 / ro_H3 - p / ro)));
	this->par[0].I3_mc[3] = (n_p_LISM_ / Kn_) * (nu_H4 * ((kv(u_H4) + kv(v_H4) - kv(u) - kv(v)) / 2.0 + //
			(U_H4 / U_M_H4) * (2.0 * p_H4 / ro_H4 - p / ro)));

	// Теперь считаем Мультифлюидные источники

	u_H1 = this->par[0].u_H1, v_H1 = this->par[0].v_H1, ro_H1 = this->par[0].ro_H1, p_H1 = this->par[0].p_H1;
	u_H2 = this->par[0].u_H2, v_H2 = this->par[0].v_H2, ro_H2 = this->par[0].ro_H2, p_H2 = this->par[0].p_H2;
	u_H3 = this->par[0].u_H3, v_H3 = this->par[0].v_H3, ro_H3 = this->par[0].ro_H3, p_H3 = this->par[0].p_H3;
	u_H4 = this->par[0].u_H4, v_H4 = this->par[0].v_H4, ro_H4 = this->par[0].ro_H4, p_H4 = this->par[0].p_H4;

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

	this->par[0].I1_mf[0] = (n_p_LISM_ / Kn_) * (nu_H1 * (u_H1 - u));
	this->par[0].I1_mf[1] = (n_p_LISM_ / Kn_) * (nu_H2 * (u_H2 - u));
	this->par[0].I1_mf[2] = (n_p_LISM_ / Kn_) * (nu_H3 * (u_H3 - u));
	this->par[0].I1_mf[3] = (n_p_LISM_ / Kn_) * (nu_H4 * (u_H4 - u));

	this->par[0].I2_mf[0] = (n_p_LISM_ / Kn_) * (nu_H1 * (v_H1 - v));
	this->par[0].I2_mf[1] = (n_p_LISM_ / Kn_) * (nu_H2 * (v_H2 - v));
	this->par[0].I2_mf[2] = (n_p_LISM_ / Kn_) * (nu_H3 * (v_H3 - v));
	this->par[0].I2_mf[3] = (n_p_LISM_ / Kn_) * (nu_H4 * (v_H4 - v));


	this->par[0].I3_mf[0] = (n_p_LISM_ / Kn_) * (nu_H1 * ((kv(u_H1) + kv(v_H1) - kv(u) - kv(v)) / 2.0 + //
		(U_H1 / U_M_H1) * (2.0 * p_H1 / ro_H1 - p / ro)));
	this->par[0].I3_mf[1] = (n_p_LISM_ / Kn_) * (nu_H2 * ((kv(u_H2) + kv(v_H2) - kv(u) - kv(v)) / 2.0 + //
		(U_H2 / U_M_H2) * (2.0 * p_H2 / ro_H2 - p / ro)));
	this->par[0].I3_mf[2] = (n_p_LISM_ / Kn_) * (nu_H3 * ((kv(u_H3) + kv(v_H3) - kv(u) - kv(v)) / 2.0 + //
		(U_H3 / U_M_H3) * (2.0 * p_H3 / ro_H3 - p / ro)));
	this->par[0].I3_mf[3] = (n_p_LISM_ / Kn_) * (nu_H4 * ((kv(u_H4) + kv(v_H4) - kv(u) - kv(v)) / 2.0 + //
		(U_H4 / U_M_H4) * (2.0 * p_H4 / ro_H4 - p / ro)));
}


void Cell::Get_Sourse_MK1(double& q1, double& q2, double& q3, const double& u, const double& v, const double& ro, const double& p, bool interpol)
{
	/*q1 = 0.0;
	q2 = 0.0;
	q3 = 0.0;
	return;*/

	if (false)//(xx < -0.8 && xx > -2.0)
	{
		q1 = this->par[0].I_u;// +this->par[0].II_u;
		q2 = this->par[0].I_v;// +this->par[0].II_v;
		q3 = this->par[0].I_T;// +this->par[0].II_T;
		return;
	}

	double xx, yy;
	this->Get_Center(xx, yy);
	double rr = sqrt(kv(xx) + kv(yy));

	


	double r1, r2;
	Cell* K_do, * K_posle;


	double u_H[pop_atom];
	double v_H[pop_atom];
	double ro_H[pop_atom];
	double p_H[pop_atom];
	double T_H[pop_atom];
	double U_M_H[pop_atom];
	double U_H[pop_atom];
	double sigma_H[pop_atom];
	double nu_H[pop_atom];

	double ku = this->par[0].k_u;
	double kv = this->par[0].k_v;
	double kT = this->par[0].k_T;

	if ((rr > 10.0 / RR_ && rr < 1500.0 / RR_ && xx >= 0.0)&&(interpol == true))
	{
		K_do = this->Back;
		K_posle = this->Next;
		if (rr < this->r_istoch)
		{
			r2 = this->r_istoch;
			r1 = this->r_istoch;
			while (r1 > rr)
			{
				r1 = K_do->r_istoch;
				K_do = K_do->Back;
			}
			K_do = K_do->Next;
			K_posle = K_do->Next;
		}
		else
		{
			r2 = this->r_istoch;
			r1 = this->r_istoch;
			while (r2 < rr)
			{
				r2 = K_posle->r_istoch;
				K_posle = K_posle->Next;
			}
			K_posle = K_posle->Back;
			K_do = K_posle->Back;
		}

		r1 = K_do->r_istoch;
		r2 = K_posle->r_istoch;

		for (int i = 0; i < pop_atom; i++)
		{
			u_H[i] = linear(r1, K_do->par[0].H_u[i], r2, K_posle->par[0].H_u[i], rr);
			v_H[i] = linear(r1, K_do->par[0].H_v[i], r2, K_posle->par[0].H_v[i], rr);
			ro_H[i] = linear(r1, K_do->par[0].H_n[i], r2, K_posle->par[0].H_n[i], rr);
			T_H[i] = linear(r1, K_do->par[0].H_T[i], r2, K_posle->par[0].H_T[i], rr);
			p_H[i] = 0.5 * T_H[i] * ro_H[i];
		}


		/*ku = linear(r1, K_do->par[0].k_u, r2, K_posle->par[0].k_u, rr);
		kv = linear(r1, K_do->par[0].k_v, r2, K_posle->par[0].k_v, rr);
		kT = linear(r1, K_do->par[0].k_T, r2, K_posle->par[0].k_T, rr);*/


	}
	else
	{
		for (int i = 0; i < pop_atom; i++)
		{
			u_H[i] = this->par[0].H_u[i];
			v_H[i] = this->par[0].H_v[i];
			ro_H[i] = this->par[0].H_n[i];
			T_H[i] = this->par[0].H_T[i];
			p_H[i] = 0.5 * T_H[i] * ro_H[i];
		}
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

	this->par[0].M_u = 0.0;
	this->par[0].M_v = 0.0;
	this->par[0].M_T = 0.0;

	for (int i = 0; i < pop_atom; i++)
	{
		this->par[0].M_u = this->par[0].M_u + (n_p_LISM_ / Kn_) * nu_H[i] * (u_H[i] - u);
		this->par[0].M_v = this->par[0].M_v + (n_p_LISM_ / Kn_) * nu_H[i] * (v_H[i] - v);
		this->par[0].M_T = this->par[0].M_T + (n_p_LISM_ / Kn_) * nu_H[i] * ((kv(u_H[i]) + kv(v_H[i]) - kv(u) - kv(v)) / 2.0 + //
			(U_H[i] / U_M_H[i]) * (2.0 * p_H[i] / ro_H[i] - p / ro));
	}


	q1 = this->par[0].M_u * ku;// +this->par[0].II_u;
	q2 = this->par[0].M_v * kv;// +this->par[0].II_v;
	q3 = this->par[0].M_T * kT;// +this->par[0].II_T;
}

void Cell::Get_Sourse_MK2(double& q1, double& q2, double& q3, const double& u, const double& v, const double& ro, const double& p)
{
	double U_M_H1, U_M_H2, U_M_H3, U_M_H4;
	double U_H1, U_H2, U_H3, U_H4;
	double sigma_H1, sigma_H2, sigma_H3, sigma_H4;
	double nu_H1, nu_H2, nu_H3, nu_H4;

	double u_H1 = this->par[0].H_u2[0], v_H1 = this->par[0].H_v2[0], ro_H1 = this->par[0].H_n2[0], p_H1 = 0.5 * this->par[0].H_T2[0] * this->par[0].H_n2[0];
	double u_H2 = this->par[0].H_u2[1], v_H2 = this->par[0].H_v2[1], ro_H2 = this->par[0].H_n2[1], p_H2 = 0.5 * this->par[0].H_T2[1] * this->par[0].H_n2[1];
	double u_H3 = this->par[0].H_u2[2], v_H3 = this->par[0].H_v2[2], ro_H3 = this->par[0].H_n2[2], p_H3 = 0.5 * this->par[0].H_T2[2] * this->par[0].H_n2[2];
	double u_H4 = this->par[0].H_u2[3], v_H4 = this->par[0].H_v2[3], ro_H4 = this->par[0].H_n2[3], p_H4 = 0.5 * this->par[0].H_T2[3] * this->par[0].H_n2[3];


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

	this->par[0].M_u = (n_p_LISM_ / Kn_) * (nu_H1 * (u_H1 - u) + nu_H2 * (u_H2 - u) //
		+ nu_H3 * (u_H3 - u) + nu_H4 * (u_H4 - u));
	this->par[0].M_v = (n_p_LISM_ / Kn_) * (nu_H1 * (v_H1 - v) + nu_H2 * (v_H2 - v) //
		+ nu_H3 * (v_H3 - v) + nu_H4 * (v_H4 - v));
	this->par[0].M_T = (n_p_LISM_ / Kn_) * (nu_H1 * ((kv(u_H1) + kv(v_H1) - kv(u) - kv(v)) / 2.0 + //
		(U_H1 / U_M_H1) * (2.0 * p_H1 / ro_H1 - p / ro)) + //
		nu_H2 * ((kv(u_H2) + kv(v_H2) - kv(u) - kv(v)) / 2.0 + //
			(U_H2 / U_M_H2) * (2.0 * p_H2 / ro_H2 - p / ro)) + //
		nu_H3 * ((kv(u_H3) + kv(v_H3) - kv(u) - kv(v)) / 2.0 + //
			(U_H3 / U_M_H3) * (2.0 * p_H3 / ro_H3 - p / ro)) + //
		nu_H4 * ((kv(u_H4) + kv(v_H4) - kv(u) - kv(v)) / 2.0 + //
			(U_H4 / U_M_H4) * (2.0 * p_H4 / ro_H4 - p / ro)));


	q1 = this->par[0].M_u * this->par[0].k_u2;
	q2 = this->par[0].M_v * this->par[0].k_v2;
	q3 = this->par[0].M_T * this->par[0].k_T2;
}

void Cell::Get_Sourse_MK3(double& q1, double& q2, double& q3, const double& u, const double& v, const double& ro, const double& p)
{
	double U_M_H1, U_M_H2, U_M_H3, U_M_H4;
	double U_H1, U_H2, U_H3, U_H4;
	double sigma_H1, sigma_H2, sigma_H3, sigma_H4;
	double nu_H1, nu_H2, nu_H3, nu_H4;

	double u_H1 = this->par[0].H_u3[0], v_H1 = this->par[0].H_v3[0], ro_H1 = this->par[0].H_n3[0], p_H1 = 0.5 * this->par[0].H_T3[0] * this->par[0].H_n3[0];
	double u_H2 = this->par[0].H_u3[1], v_H2 = this->par[0].H_v3[1], ro_H2 = this->par[0].H_n3[1], p_H2 = 0.5 * this->par[0].H_T3[1] * this->par[0].H_n3[1];
	double u_H3 = this->par[0].H_u3[2], v_H3 = this->par[0].H_v3[2], ro_H3 = this->par[0].H_n3[2], p_H3 = 0.5 * this->par[0].H_T3[2] * this->par[0].H_n3[2];
	double u_H4 = this->par[0].H_u3[3], v_H4 = this->par[0].H_v3[3], ro_H4 = this->par[0].H_n3[3], p_H4 = 0.5 * this->par[0].H_T3[3] * this->par[0].H_n3[3];


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

	this->par[0].M_u = (n_p_LISM_ / Kn_) * (nu_H1 * (u_H1 - u) + nu_H2 * (u_H2 - u) //
		+ nu_H3 * (u_H3 - u) + nu_H4 * (u_H4 - u));
	this->par[0].M_v = (n_p_LISM_ / Kn_) * (nu_H1 * (v_H1 - v) + nu_H2 * (v_H2 - v) //
		+ nu_H3 * (v_H3 - v) + nu_H4 * (v_H4 - v));
	this->par[0].M_T = (n_p_LISM_ / Kn_) * (nu_H1 * ((kv(u_H1) + kv(v_H1) - kv(u) - kv(v)) / 2.0 + //
		(U_H1 / U_M_H1) * (2.0 * p_H1 / ro_H1 - p / ro)) + //
		nu_H2 * ((kv(u_H2) + kv(v_H2) - kv(u) - kv(v)) / 2.0 + //
			(U_H2 / U_M_H2) * (2.0 * p_H2 / ro_H2 - p / ro)) + //
		nu_H3 * ((kv(u_H3) + kv(v_H3) - kv(u) - kv(v)) / 2.0 + //
			(U_H3 / U_M_H3) * (2.0 * p_H3 / ro_H3 - p / ro)) + //
		nu_H4 * ((kv(u_H4) + kv(v_H4) - kv(u) - kv(v)) / 2.0 + //
			(U_H4 / U_M_H4) * (2.0 * p_H4 / ro_H4 - p / ro)));


	q1 = this->par[0].M_u * this->par[0].k_u3;
	q2 = this->par[0].M_v * this->par[0].k_v3;
	q3 = this->par[0].M_T * this->par[0].k_T3;
}



double Cell::get_nu_pui(const double& L, int i)
{
	if (L > L_Igor)
	{
		cout << "PROBLEM  Cell  833  nu_pui > 20   " << L << endl;
		exit(-1);
	}

	double dl = L_Igor / (k_Igor - 1);
	int k1 = (int)(L / dl);
	if (k1 == (k_Igor - 1))
	{
		return  this->nu_pui[i][k1];
	}
	else
	{
		int k2 = k1 + 1;
		double al = (L - k1 * dl) / dl;
		return this->nu_pui[i][k1] * (1.0 - al) + this->nu_pui[i][k2] * al;
	}

}

double Cell::get_nu_pui2(const double& L, int i)
{
	if (L > L_Igor)
	{
		cout << "PROBLEM  Cell  833  nu_pui > 20" << endl;
		exit(-1);
	}

	double dl = L_Igor / (k_Igor - 1);
	int k1 = (int)(L / dl);
	if (k1 == (k_Igor - 1))
	{
		return  this->nu2_pui[i][k1];
	}
	else
	{
		int k2 = k1 + 1;
		double al = (L - k1 * dl) / dl;
		return  this->nu2_pui[i][k1] * (1.0 - al) + this->nu2_pui[i][k2] * al;
	}

}

double Cell::get_nu_pui3(const double& L, int i)
{
	if (L > L_Igor)
	{
		cout << "PROBLEM  Cell  833  nu_pui > 20" << endl;
		exit(-1);
	}

	double dl = L_Igor / (k_Igor - 1);
	int k1 = (int)(L / dl);
	if (k1 == (k_Igor - 1))
	{
		return  this->nu3_pui[i][k1];
	}
	else
	{
		int k2 = k1 + 1;
		double al = (L - k1 * dl) / dl;
		return  this->nu3_pui[i][k1] * (1.0 - al) + this->nu3_pui[i][k2] * al;
	}

}



double Cell::get_fpui(const double& W, const double& Wmin, const double& Wmax)
{
	if (W < Wmin || W > Wmax)
	{
		return 0.0;
	}

	int N = this->fpui.size();

	//int i = (int)((N - 1) * log(W / Wmin) / log(Wmax / Wmin));
	int i = (int)((W - Wmin) / (Wmax - Wmin) * (N - 1));

	if (i >= N - 1)
	{
		return this->fpui[i];
	}
	else
	{
		//double WL = Wmin * pow(Wmax / Wmin, 1.0 * i / (N - 1));
		double WL = Wmin + (Wmax - Wmin) * i / (N - 1);
		//double WR = Wmin * pow(Wmax / Wmin, 1.0 * (i + 1) / (N - 1));
		double WR = Wmin + (Wmax - Wmin) * (i + 1) / (N - 1);

		return (this->fpui[i + 1] - this->fpui[i]) / (WR - WL) * (W - WL) + this->fpui[i];
	}
}

bool Cell::Change_Velosity_PUI(Sensor* sens, const double& Vh1, const double& Vh2, const double& Vh3, //
	const double& Vp1, const double& Vp2, const double& Vp3, double& W1, double& W2, double& W3, int Nw, //
	const double& wmin, const double& wmax, int ifg)
{
	// Число делений в распределениях ИГОРЯ   NW = 100 было

	double u = sqrt(kv(Vh1 - Vp1) + kv(Vh2 - Vp2) + kv(Vh3 - Vp3));
	double Wr;

	double Wmax3 = this->Wmax_sort[ifg] * this->Wmax_sort[ifg] * this->Wmax_sort[ifg];
	double Wmin3 = this->Wmin_sort[ifg] * this->Wmin_sort[ifg] * this->Wmin_sort[ifg];

	double Thetta = 0.0;
	double phi = 0.0;
	double h;
	double M = (this->Wmax_sort[ifg] + u) * sigma(fabs(this->Wmax_sort[ifg] - u)) * this->fpui_max_sort[ifg];
	double uu = 0.0;
	double coss = 0.0;

	// Находим Wr и Thetta

	do
	{
		//cout << "DO " << endl;
		Wr = pow(sens->MakeRandom() * (Wmax3 - Wmin3) + Wmin3, 1.0 / 3.0);
		coss = 1.0 - 2.0 * sens->MakeRandom();
		uu = sqrt(Wr * Wr + u * u - 2.0 * Wr * u * coss);
		h = uu * sigma(uu) * this->get_fpui(Wr, wmin, wmax) / M;
		//cout << h << "    " << Wr << "  " << uu << "    " << this->get_fpui(Wr, wmin, wmax) << endl;
		//cout << "   " << M << "   " << this->fpui_max << "  " << uu * sigma(uu) * this->get_fpui(Wr, wmin, wmax) << endl;
		//cout << this->number << endl;
		if (h > 1.0)
		{
			cout << " h > 1  iffiuehrkf    " << h << endl;
			exit(-1);
		}
	} while (h < sens->MakeRandom());
	//cout << "ENDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD" << endl;
	Thetta = acos(coss);

	phi = 2.0 * pi_ * sens->MakeRandom();

	double ez1, ez2, ez3;
	double ex1, ex2, ex3, ex;
	double ey1, ey2, ey3;

	ez1 = Vh1 - Vp1;
	ez2 = Vh2 - Vp2;
	ez3 = Vh3 - Vp3;
	ez1 = ez1 / u;
	ez2 = ez2 / u;
	ez3 = ez3 / u;

	ex = sqrt(kv(Vp1) + kv(Vp2) + kv(Vp3));

	// Находим какой-то базис
	if (fabs(ez1 * Vp1 + ez2 * Vp2 + ez3 * Vp3) / ex < 0.9)
	{
		Vector_product(ez1, ez2, ez3, Vp1, Vp2, Vp3, ex1, ex2, ex3);
		ex = sqrt(kv(ex1) + kv(ex2) + kv(ex3));
		ex1 = ex1 / ex;
		ex2 = ex2 / ex;
		ex3 = ex3 / ex;
		Vector_product(ez1, ez2, ez3, ex1, ex2, ex3, ey1, ey2, ey3);
	}
	else
	{
		if (fabs(ez1 - 1.0) > 0.001)
		{
			Vector_product(ez1, ez2, ez3, 1.0, 0.0, 0.0, ex1, ex2, ex3);
			ex = sqrt(kv(ex1) + kv(ex2) + kv(ex3));
			ex1 = ex1 / ex;
			ex2 = ex2 / ex;
			ex3 = ex3 / ex;
			Vector_product(ez1, ez2, ez3, ex1, ex2, ex3, ey1, ey2, ey3);
		}
		else if (fabs(ez2 - 1.0) > 0.001)
		{
			Vector_product(ez1, ez2, ez3, 0.0, 1.0, 0.0, ex1, ex2, ex3);
			ex = sqrt(kv(ex1) + kv(ex2) + kv(ex3));
			ex1 = ex1 / ex;
			ex2 = ex2 / ex;
			ex3 = ex3 / ex;
			Vector_product(ez1, ez2, ez3, ex1, ex2, ex3, ey1, ey2, ey3);
		}
		else if (fabs(ez3 - 1.0) > 0.001)
		{
			Vector_product(ez1, ez2, ez3, 0.0, 0.0, 1.0, ex1, ex2, ex3);
			ex = sqrt(kv(ex1) + kv(ex2) + kv(ex3));
			ex1 = ex1 / ex;
			ex2 = ex2 / ex;
			ex3 = ex3 / ex;
			Vector_product(ez1, ez2, ez3, ex1, ex2, ex3, ey1, ey2, ey3);
		}
		else
		{
			cout << "ERORORORORO  995 wijrihfberfr" << endl;
			cout << ez1 << " " << ez2 << " " << ez3 << endl;
			cout << Vp1 << " " << Vp2 << " " << Vp3 << endl;
			cout << Vh1 << " " << Vh2 << " " << Vh3 << endl;
			exit(-2);
		}
	}

	double V1, V2, V3;
	V1 = Wr * sin(Thetta) * cos(phi);
	V2 = Wr * sin(Thetta) * sin(phi);
	V3 = Wr * cos(Thetta);

	W1 = Vp1 + V1 * ex1 + V2 * ey1 + V3 * ez1;
	W2 = Vp2 + V1 * ex2 + V2 * ey2 + V3 * ez2;
	W3 = Vp3 + V1 * ex3 + V2 * ey3 + V3 * ez3;

	if (fpclassify(W1) == FP_NAN || fpclassify(W1) == FP_SUBNORMAL || fpclassify(W1) == FP_INFINITE)
	{
		cout << "NUNUNUNUN  1012" << endl;
		cout << Wr << " " << Thetta << " " << phi << endl;
		exit(-1);
	}

	return true;
}


bool Cell::Change_Velosity_PUI2(Sensor* sens, const double& Vh1, const double& Vh2, const double& Vh3, //
	const double& Vp1, const double& Vp2, const double& Vp3, double& W1, double& W2, double& W3, int Nw, //
	const double& wmin, const double& wmax, int ifg)
{
	// Число делений в распределениях ИГОРЯ   NW = 100 было

	double u = sqrt(kv(Vh1 - Vp1) + kv(Vh2 - Vp2) + kv(Vh3 - Vp3));
	double nu = this->get_nu_pui(u, ifg);
	int n1 = 45;                            // Понижает точность но увеличивает скорость счёта
	// Лучше чтобы n1 = тому n1, когда считали частоту в начале
	double dthe = pi_ / (n1 - 1);
	double L, R, d, f1, f2;
	double LL = u;
	double Int1 = 0.0;
	double Int_do = 0.0;
	double ksi = sens->MakeRandom();
	double Wr, Wt, Wp;

	// Находим Wr
	for (int ik = 0; ik < Nw - 1; ik++)
	{
		//L = wmin * pow(wmax / wmin, 1.0 * ik / (Nw - 1));
		//R = wmin * pow(wmax / wmin, 1.0 * (ik + 1) / (Nw - 1));
		L = wmin + (wmax - wmin) * ik / (Nw - 1);
		R = wmin + (wmax - wmin) * (ik + 1) / (Nw - 1);
		for (int ij = 0; ij < n1; ij++)
		{
			double the = dthe * ij;
			d = sqrt(L * L + LL * LL - 2.0 * L * LL * cos(the));
			f1 = this->fpui[ik] * d * sigma(d) * L * L * sin(the);

			d = sqrt(R * R + LL * LL - 2.0 * R * LL * cos(the));
			f2 = this->fpui[ik + 1] * d * sigma(d) * R * R * sin(the);
			Int1 = Int1 + 0.5 * (R - L) * (f2 + f1) * dthe / nu;
		}

		if (Int1 >= ksi)
		{
			double k = (Int1 - Int_do) / (R - L);
			Wr = L + (ksi - Int_do) / k;
			break;
		}

		if (ik == Nw - 2)
		{
			Wr = R;
		}

		Int_do = Int1;
	}

	// Находим  Tetta                 ПРОВЕРИТЬ
	double Thetta = 0.0;
	double phi = 0.0;
	double ksi2;
	do
	{
		ksi2 = 1.0 - 2.0 * sens->MakeRandom();
		Thetta = acos(ksi2);
		d = sqrt(kv(Wr) + kv(u) - 2.0 * Wr * u * ksi2);
	} while (d * sigma(d) / (fabs(Wr + u) * sigma(fabs(Wr - u))) < sens->MakeRandom());

	phi = 2.0 * pi_ * sens->MakeRandom();

	double ez1, ez2, ez3;
	double ex1, ex2, ex3, ex;
	double ey1, ey2, ey3;

	ez1 = Vh1 - Vp1;
	ez2 = Vh2 - Vp2;
	ez3 = Vh3 - Vp3;
	ez1 = ez1 / u;
	ez2 = ez2 / u;
	ez3 = ez3 / u;

	ex = sqrt(kv(Vp1) + kv(Vp2) + kv(Vp3));

	// Находим какой-то базис
	if (fabs(ez1 * Vp1 + ez2 * Vp2 + ez3 * Vp3) / ex < 0.9)
	{
		Vector_product(ez1, ez2, ez3, Vp1, Vp2, Vp3, ex1, ex2, ex3);
		ex = sqrt(kv(ex1) + kv(ex2) + kv(ex3));
		ex1 = ex1 / ex;
		ex2 = ex2 / ex;
		ex3 = ex3 / ex;
		Vector_product(ez1, ez2, ez3, ex1, ex2, ex3, ey1, ey2, ey3);
	}
	else
	{
		if (fabs(ez1 - 1.0) > 0.001)
		{
			Vector_product(ez1, ez2, ez3, 1.0, 0.0, 0.0, ex1, ex2, ex3);
			ex = sqrt(kv(ex1) + kv(ex2) + kv(ex3));
			ex1 = ex1 / ex;
			ex2 = ex2 / ex;
			ex3 = ex3 / ex;
			Vector_product(ez1, ez2, ez3, ex1, ex2, ex3, ey1, ey2, ey3);
		}
		else if (fabs(ez2 - 1.0) > 0.001)
		{
			Vector_product(ez1, ez2, ez3, 0.0, 1.0, 0.0, ex1, ex2, ex3);
			ex = sqrt(kv(ex1) + kv(ex2) + kv(ex3));
			ex1 = ex1 / ex;
			ex2 = ex2 / ex;
			ex3 = ex3 / ex;
			Vector_product(ez1, ez2, ez3, ex1, ex2, ex3, ey1, ey2, ey3);
		}
		else if (fabs(ez3 - 1.0) > 0.001)
		{
			Vector_product(ez1, ez2, ez3, 0.0, 0.0, 1.0, ex1, ex2, ex3);
			ex = sqrt(kv(ex1) + kv(ex2) + kv(ex3));
			ex1 = ex1 / ex;
			ex2 = ex2 / ex;
			ex3 = ex3 / ex;
			Vector_product(ez1, ez2, ez3, ex1, ex2, ex3, ey1, ey2, ey3);
		}
		else
		{
			cout << "ERORORORORO  995 wijrihfberfr" << endl;
			cout << ez1 << " " << ez2 << " " << ez3 << endl;
			cout << Vp1 << " " << Vp2 << " " << Vp3 << endl;
			cout << Vh1 << " " << Vh2 << " " << Vh3 << endl;
			exit(-2);
		}
	}

	double V1, V2, V3;
	V1 = Wr * sin(Thetta) * cos(phi);
	V2 = Wr * sin(Thetta) * sin(phi);
	V3 = Wr * cos(Thetta);

	W1 = Vp1 + V1 * ex1 + V2 * ey1 + V3 * ez1;
	W2 = Vp2 + V1 * ex2 + V2 * ey2 + V3 * ez2;
	W3 = Vp3 + V1 * ex3 + V2 * ey3 + V3 * ez3;

	if (fpclassify(W1) == FP_NAN || fpclassify(W1) == FP_SUBNORMAL || fpclassify(W1) == FP_INFINITE)
	{
		cout << "NUNUNUNUN  1012" << endl;
		cout << Wr << " " << Thetta << " " << phi << endl;
		cout << nu << endl;
		cout << Int1 << " " << Int_do << " " << ksi << endl;
		cout << L << " " << R << endl;
		for (int ijk = 0; ijk < 100; ijk++)
		{
			cout << this->fpui[ijk] << "  ";
		}
		exit(-1);
	}

	return true;
}


