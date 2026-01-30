#include "Gran.h"
#include <math.h>

Gran::Gran(Point* A, Point* B, Gran_type type)
{
	this->A = A;
	this->B = B;
	this->type = type;
	this->Gran_copy = nullptr;
	this->Master = nullptr;
	this->Sosed = nullptr;
	this->number = -1;
	this->main_gran = true;
	this->Sosed_down = nullptr;
	this->Sosed_up = nullptr;
	this->metod_HLLC = 1;
	if (fabs(A->x - B->x) > 0.00001 / RR_)
	{
		this->parallel = false;
		this->a = (A->y - B->y) / (A->x - B->x);
		this->b = (A->y * B->x - B->y * A->x) / (B->x - A->x);
	}
	else
	{
		this->parallel = true;
		this->a = 0.0;
		this->b = 0.0;
	}

	if (this->parallel == false)
	{
		double n1, n2;
		this->Get_normal(n1, n2);
		double x, y;
		x = A->x + n1;
		y = A->y + n2;
		if (y - this->a * x - this->b > 0)
		{
			this->koef = 1;
		}
		else
		{
			this->koef = -1;
		}
	}
	else
	{
		this->koef = 0;
	}
}

void Gran::renew(void)
{
	if (fabs(A->x - B->x) > 0.00001 /RR_)
	{
		this->parallel = false;
		this->a = (A->y - B->y) / (A->x - B->x);
		this->b = (A->y * B->x - B->y * A->x) / (B->x - A->x);
	}
	else
	{
		this->parallel = true;
		this->a = 0.0;
		this->b = 0.0;
	}

	if (this->parallel == false)
	{
		double n1, n2;
		this->Get_normal(n1, n2);
		double x, y;
		x = A->x + n1;
		y = A->y + n2;
		if (y - this->a * x - this->b > 0)
		{
			this->koef = 1;
		}
		else
		{
			this->koef = -1;
		}
	}
	else
	{
		this->koef = 0;
	}

	this->aa = A->y - B->y;
	this->bb = B->x - A->x;
	double nnn = sqrt(kvv(0.0, this->aa, this->bb));
	this->aa = this->aa / nnn;
	this->bb = this->bb / nnn;
	this->cc = -this->aa * A->x - this->bb * A->y;
}


bool Gran::belong(const double& x, const double& y)
{
	if (this->type == Axis)
	{
		return true;
	}

	if (this->parallel == false)
	{
		if (this->koef * (y - this->a * x - this->b) <= 0.0)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	else
	{
		double n1, n2;
		this->Get_normal(n1, n2);
		if ((x - this->A->x) / n1 < 0)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
}

double Gran::Get_lenght(void)
{
	return sqrt(kv(this->A->x - this->B->x) + kv(this->B->y - this->A->y));
}

bool Gran::belong_gran(const double& x, const double& y)
{
	if (this->parallel == false)
	{
		if (fabs(y - this->a * x - this->b) / (sqrt(1.0 + kv(this->a))) <= 0.0001 / RR_)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	else
	{
		if (fabs(x - this->A->x) <= 0.0001 / RR_ && max(this->A->y, this->B->y) >= y && min(this->A->y, this->B->y) <= y)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
}

void Gran::Get_Center(double& x, double& y)
{
	x = (this->A->x + this->B->x) / 2.0;
	y = (this->A->y + this->B->y) / 2.0;
}

void Gran::Get_Center_posle(double& x, double& y)
{
	x = (this->A->x2 + this->B->x2) / 2.0;
	y = (this->A->y2 + this->B->y2) / 2.0;
}

void Gran::Get_normal(double& n1, double& n2)
{
	double t1 = this->B->x - this->A->x;
	double t2 = this->B->y - this->A->y;
	double n = sqrt(kv(t1) + kv(t2));
	t1 = t1 / n;
	t2 = t2 / n;
	n1 = t2;
	n2 = -t1;
	return;
}

void Gran::Get_normal_posle(double& n1, double& n2)
{
	double t1 = this->B->x2 - this->A->x2;
	double t2 = this->B->y2 - this->A->y2;
	double n = sqrt(kv(t1) + kv(t2));
	t1 = t1 / n;
	t2 = t2 / n;
	n1 = t2;
	n2 = -t1;
	return;
}

double Gran::Get_square(void)
{
	return sqrt(kv(this->B->x - this->A->x) + kv(this->B->y - this->A->y));
}

double Gran::Get_square_rotate(const double& angle)
{
	return pi_ * sqrt(kv(this->B->x - this->A->x) + kv(this->B->y - this->A->y)) * 0.5 * (this->B->y + this->A->y) * angle/180.0;
}

void Gran::Get_par(Parametr& par, int i)  // Здесь задаются граничные условия
{

	if (this->type == Usualy)  // Не надо менять
	{
		par = this->Sosed->par[i];
	}
	else if (this->type == Extern)  // Не надо менять
	{
		par = this->Master->par[i];
		//par.p = par.p / 1000000.0;
		if (par.u >= Velosity_inf)                                                         // ОТСОС ЖИДКОСТИ!!!!!!!!!!!!!!!!!!!!!!!!!!
		{
			par.u = Velosity_inf;
		}
	}
	else if (this->type == Upper_wall)  // Не надо менять
	{
		par = this->Master->par[i];
	}
	else if (this->type == Axis)  // Не надо менять
	{
		par = this->Master->par[i];
		par.v = -par.v;
		par.v_H1 = -par.v_H1;
		par.v_H2 = -par.v_H2;
		par.v_H3 = -par.v_H3;
		par.v_H4 = -par.v_H4;
	}
	else if (this->type == Inner_sphere)
	{
		auto par2 = this->Master->par[i];   // this->Sosed->par[i];
		double x, y;
		this->Get_Center(x, y);
		double dist = sqrt(x * x + y * y);
		double ro = 1.0/(kv(chi_real) * dist * dist);
		double P_E = ro * chi_real * chi_real / (ggg * 10.0 * 10.0);
		double T_p = (P_E * pow(1.0 / dist, 2.0 * ggg)) / (2.0 * ro / (dist * dist));
		//par = { ro / (dist * dist), P_E * pow(1.0 / dist, 2.0 * ggg), chi_ * x / dist, chi_ * y / dist, ro / (dist * dist),//
		//0.0003, (0.0003 * chi_real * chi_real / (ggg * 5.0 * 5.0)) * pow(1.0 / dist, 2.0 * ggg), chi_real* x / dist, chi_real* y / dist,//
		//	par2.ro_H2, par2.p_H2, par2.u_H2, par2.v_H2 , //
		//	par2.ro_H3, par2.p_H3, par2.u_H3, par2.v_H3,//
		//	par2.ro_H4, par2.p_H4, par2.u_H4, par2.v_H4};

		//par = { ro / (dist * dist), P_E * pow(1.0 / dist, 2.0 * ggg), chi_ * x / dist, chi_ * y / dist, ro / (dist * dist),//
		//0.001666666, (0.001666666 * chi_real * chi_real / (ggg * 14.1344 * 14.1344)), chi_real* x / dist, chi_real* y / dist,//
		//	0.000001, 0.0000001, par2.u_H2, par2.v_H2 , //
		//	0.000001, 0.0000001, par2.u_H3, par2.v_H3,//
		//	0.000001, 0.0000001, par2.u_H4, par2.v_H4 };

		//par = { ro / kv(chi_real / chi_) / (dist * dist), P_E * pow(1.0 / dist, 2.0 * ggg), chi_ * x * (chi_real / chi_) / dist, chi_ * y * (chi_real / chi_) / dist, ro / kv(chi_real / chi_) / (dist * dist),//
		//0.001666666, (0.001666666 * chi_real * chi_real / (ggg * 14.1344 * 14.1344)), chi_real * x / dist, chi_real * y / dist,//
		//	0.00000001, 0.0000001, par2.u_H2, par2.v_H2 , //
		//	0.00000001, 0.0000001, par2.u_H3, par2.v_H3,//
		//	0.00000001, 0.0000001, par2.u_H4, par2.v_H4 };

		par = { ro, P_E, chi_real * x/dist, chi_real * y / dist, ro,//
		0.0000001, (0.0000001 * chi_real * chi_real / (ggg * 14.1344 * 14.1344)), chi_real * x / dist, chi_real * y / dist,//
			par2.ro_H2, par2.p_H2, par2.u_H2, par2.v_H2 , //
			par2.ro_H3, par2.p_H3, par2.u_H3, par2.v_H3,//
			par2.ro_H4, par2.p_H4, par2.u_H4, par2.v_H4 };

		//par = { ro * kv(chi_real / chi_), P_E, chi_real * x / dist / (chi_real / chi_), chi_real * y / dist / (chi_real / chi_), ro * (chi_real / chi_),//
		//0.0000001, (0.0000001 * chi_real * chi_real / (ggg * 14.1344 * 14.1344)), chi_real * x / dist, chi_real * y / dist,//
		//	par2.ro_H2, par2.p_H2, par2.u_H2, par2.v_H2 , //
		//	par2.ro_H3, par2.p_H3, par2.u_H3, par2.v_H3,//
		//	par2.ro_H4, par2.p_H4, par2.u_H4, par2.v_H4 };


		//cout << "Gran.cpp    " << par.ro << " " << dist << endl;
	}
	else if (this->type == Input)
	{
		auto par2 = this->Master->par[i];
		par = {1.0, 1.0, Velosity_inf, 0.0, 100.0, par2.ro_H1, par2.p_H1, par2.u_H1, par2.v_H1,//
		par2.ro_H2, par2.p_H2, par2.u_H2, par2.v_H2 , par2.ro_H3, par2.p_H3, par2.u_H3, par2.v_H3,//
		1.0, 0.5, Velosity_inf, 0.0};

		//par = { 4.0, 2.5, Velosity_inf, 0.0, 100.0, par2.ro_H1, par2.p_H1, par2.u_H1, par2.v_H1,//
		//par2.ro_H2, par2.p_H2, par2.u_H2, par2.v_H2 , par2.ro_H3, par2.p_H3, par2.u_H3, par2.v_H3,//
		//1.0, 0.5, Velosity_inf, 0.0 };
	}

	return;
}


void Gran::Get_par_TVD(Parametr& par, int i)  // Здесь задаются граничные условия
// Сносит значение в текущей ячейке на текущую гарнь
{
	if (this->type == Usualy)  // Не надо менять
	{
		if (this->Sosed_down == nullptr)
		{
			bool b = false;
			for (auto& i : this->Master->Grans)
			{
				if (i->type == Axis)
				{
					b = true;
					break;
				}
			}
			if (b == true)
			{
				double x, y, x2, y2, dist1, dist2, dist3;
				par = this->Master->par[i];
				this->Get_Center(x, y);
				//cout << "Gran.cpp    " << x << " " << y << endl;
				auto par2 = this->Master->par[i];
				auto par1 = this->Master->par[i];
				par1.v = -par1.v;
				par1.v_H1 = -par1.v_H1;
				par1.v_H2 = -par1.v_H2;
				par1.v_H3 = -par1.v_H3;
				par1.v_H4 = -par1.v_H4;
				auto par3 = this->Sosed->par[i];
				this->Master->Get_Center(x2, y2);
				//cout << "Gran.cpp    " << x2 << " " << y2 << endl;
				dist2 = sqrt(kv(x - x2) + kv(y - y2));
				this->Sosed->Get_Center(x2, y2);
				//cout << "Gran.cpp    " << x2 << " " << y2 << endl;
				dist3 = sqrt(kv(x - x2) + kv(y - y2));
				this->Master->Get_Center(x2, y2);
				y2 = -y2;
				//cout << "Gran.cpp    " << x2 << " " << y2 << endl;
				dist1 = sqrt(kv(x - x2) + kv(y - y2));
				par.v = linear(-dist1, par1.v, -dist2, par2.v, dist3, par3.v, 0.0);
				par.v_H1 = linear(-dist1, par1.v_H1, -dist2, par2.v_H1, dist3, par3.v_H1, 0.0);
				par.v_H2 = linear(-dist1, par1.v_H2, -dist2, par2.v_H2, dist3, par3.v_H2, 0.0);
				par.v_H3 = linear(-dist1, par1.v_H3, -dist2, par2.v_H3, dist3, par3.v_H3, 0.0);
				par.v_H4 = linear(-dist1, par1.v_H4, -dist2, par2.v_H4, dist3, par3.v_H4, 0.0);
				//cout << "Gran.cpp    " << " __________ " << endl;
			}
			else
			{
				par = this->Master->par[i];
			}
			return;
		}
		if ((this->A->type == P_Contact && this->B->type == P_Contact) ||//
			(this->A->type == P_Inner_shock && this->B->type == P_Inner_shock) ||//
			(this->A->type == P_Outer_shock && this->B->type == P_Outer_shock))
		{
			double x, y, x2, y2, dist1, dist2;
			this->Get_Center(x, y);
			auto par2 = this->Master->par[i];
			auto par1 = this->Sosed_down->par[i];
			this->Master->Get_Center(x2, y2);
			dist2 = sqrt(kv(x - x2) + kv(y - y2));
			this->Sosed_down->Get_Center(x2, y2);
			dist1 = sqrt(kv(x - x2) + kv(y - y2));

			par.ro = linear(-dist1, par1.ro, -dist2, par2.ro, 0.0);
			par.p = linear(-dist1, par1.p, -dist2, par2.p, 0.0);
			par.u = linear(-dist1, par1.u, -dist2, par2.u, 0.0);
			par.v = linear(-dist1, par1.v, -dist2, par2.v, 0.0);
			par.Q = linear(-dist1, par1.Q, -dist2, par2.Q, 0.0);
			if (par.ro <= 0.0)
			{
				par.ro = par2.ro;
			}
			if (par.p <= 0.0)
			{
				par.p = par2.p;
			}

			par.ro_H1 = linear(-dist1, par1.ro_H1, -dist2, par2.ro_H1, 0.0);
			par.p_H1 = linear(-dist1, par1.p_H1, -dist2, par2.p_H1, 0.0);
			par.u_H1 = linear(-dist1, par1.u_H1, -dist2, par2.u_H1, 0.0);
			par.v_H1 = linear(-dist1, par1.v_H1, -dist2, par2.v_H1, 0.0);
			par.ro_H2 = linear(-dist1, par1.ro_H2, -dist2, par2.ro_H2, 0.0);
			par.p_H2 = linear(-dist1, par1.p_H2, -dist2, par2.p_H2, 0.0);
			par.u_H2 = linear(-dist1, par1.u_H2, -dist2, par2.u_H2, 0.0);
			par.v_H2 = linear(-dist1, par1.v_H2, -dist2, par2.v_H2, 0.0);
			par.ro_H3 = linear(-dist1, par1.ro_H3, -dist2, par2.ro_H3, 0.0);
			par.p_H3 = linear(-dist1, par1.p_H3, -dist2, par2.p_H3, 0.0);
			par.u_H3 = linear(-dist1, par1.u_H3, -dist2, par2.u_H3, 0.0);
			par.v_H3 = linear(-dist1, par1.v_H3, -dist2, par2.v_H3, 0.0);
			par.ro_H4 = linear(-dist1, par1.ro_H4, -dist2, par2.ro_H4, 0.0);
			par.p_H4 = linear(-dist1, par1.p_H4, -dist2, par2.p_H4, 0.0);
			par.u_H4 = linear(-dist1, par1.u_H4, -dist2, par2.u_H4, 0.0);
			par.v_H4 = linear(-dist1, par1.v_H4, -dist2, par2.v_H4, 0.0);

			if (par.ro_H1 <= 0.0)
			{
				par.ro_H1 = par2.ro_H1;
			}
			if (par.p_H1 <= 0.0)
			{
				par.p_H1 = par2.p_H1;
			}

			if (par.ro_H2 <= 0.0)
			{
				par.ro_H2 = par2.ro_H2;
			}
			if (par.p_H2 <= 0.0)
			{
				par.p_H2 = par2.p_H2;
			}

			if (par.ro_H3 <= 0.0)
			{
				par.ro_H3 = par2.ro_H3;
			}
			if (par.p_H3 <= 0.0)
			{
				par.p_H3 = par2.p_H3;
			}

			if (par.ro_H4 <= 0.0)
			{
				par.ro_H4 = par2.ro_H4;
			}
			if (par.p_H4 <= 0.0)
			{
				par.p_H4 = par2.p_H4;
			}

			return;
		}
		double x, y, x2, y2, dist1, dist2, dist3;
		this->Get_Center(x, y);
		auto par2 = this->Master->par[i];
		auto par1 = this->Sosed_down->par[i];
		auto par3 = this->Sosed->par[i];
		this->Master->Get_Center(x2, y2);
		dist2 = sqrt(kv(x - x2) + kv(y - y2));
		this->Sosed->Get_Center(x2, y2);
		dist3 = sqrt(kv(x - x2) + kv(y - y2));
		this->Sosed_down->Get_Center(x2, y2);
		dist1 = sqrt(kv(x - x2) + kv(y - y2));

		par.ro = linear(-dist1, par1.ro, -dist2, par2.ro, dist3, par3.ro, 0.0);
		par.p = linear(-dist1, par1.p, -dist2, par2.p, dist3, par3.p, 0.0);
		par.u = linear(-dist1, par1.u, -dist2, par2.u, dist3, par3.u, 0.0);
		par.v = linear(-dist1, par1.v, -dist2, par2.v, dist3, par3.v, 0.0);
		par.Q = linear(-dist1, par1.Q, -dist2, par2.Q, dist3, par3.Q, 0.0);
		if (par.ro <= 0.0)
		{
			cout << "Gran.cpp    " << "Get par  354" << endl;
			cout << "Gran.cpp    " << x << " " << y << " " << x2 << " " << y2 << " " << par1.ro << " " << par2.ro << " " << par3.ro << endl;
			exit(-1);
		}
		if (par.p <= 0.0)
		{
			cout << "Gran.cpp    " << "Get par  361" << endl;
			cout << "Gran.cpp    " << x << " " << y << " " << x2 << " " << y2 << " " << par1.p << " " << par2.p << " " << par3.p << endl;
			exit(-1);
		}

		par.ro_H1 = linear(-dist1, par1.ro_H1, -dist2, par2.ro_H1, dist3, par3.ro_H1, 0.0);
		par.p_H1 = linear(-dist1, par1.p_H1, -dist2, par2.p_H1, dist3, par3.p_H1, 0.0);
		par.u_H1 = linear(-dist1, par1.u_H1, -dist2, par2.u_H1, dist3, par3.u_H1, 0.0);
		par.v_H1 = linear(-dist1, par1.v_H1, -dist2, par2.v_H1, dist3, par3.v_H1, 0.0);
		par.ro_H2 = linear(-dist1, par1.ro_H2, -dist2, par2.ro_H2, dist3, par3.ro_H2, 0.0);
		par.p_H2 = linear(-dist1, par1.p_H2, -dist2, par2.p_H2, dist3, par3.p_H2, 0.0);
		par.u_H2 = linear(-dist1, par1.u_H2, -dist2, par2.u_H2, dist3, par3.u_H2, 0.0);
		par.v_H2 = linear(-dist1, par1.v_H2, -dist2, par2.v_H2, dist3, par3.v_H2, 0.0);
		par.ro_H3 = linear(-dist1, par1.ro_H3, -dist2, par2.ro_H3, dist3, par3.ro_H3, 0.0);
		par.p_H3 = linear(-dist1, par1.p_H3, -dist2, par2.p_H3, dist3, par3.p_H3, 0.0);
		par.u_H3 = linear(-dist1, par1.u_H3, -dist2, par2.u_H3, dist3, par3.u_H3, 0.0);
		par.v_H3 = linear(-dist1, par1.v_H3, -dist2, par2.v_H3, dist3, par3.v_H3, 0.0);
		par.ro_H4 = linear(-dist1, par1.ro_H4, -dist2, par2.ro_H4, dist3, par3.ro_H4, 0.0);
		par.p_H4 = linear(-dist1, par1.p_H4, -dist2, par2.p_H4, dist3, par3.p_H4, 0.0);
		par.u_H4 = linear(-dist1, par1.u_H4, -dist2, par2.u_H4, dist3, par3.u_H4, 0.0);
		par.v_H4 = linear(-dist1, par1.v_H4, -dist2, par2.v_H4, dist3, par3.v_H4, 0.0);
	}
	else if (this->type == Extern)  // Не надо менять
	{
		par = this->Master->par[i];
		if (par.u >= Velosity_inf)
		{
			par.u = Velosity_inf;
		}
	}
	else if (this->type == Upper_wall)  // Не надо менять
	{
		par = this->Master->par[i];
	}
	else if (this->type == Axis)  
	{
		double x, y, x2, y2, dist1, dist2;
		this->Get_Center(x, y);
		//cout << "Gran.cpp    " << x << " " << y << endl;
		auto par2 = this->Master->par[i];
		//auto par3 = this->Master->par[i];
		//par3.v = -par3.v;
		auto par1 = this->Sosed_down->par[i];
		this->Master->Get_Center(x2, y2);
		dist2 = sqrt(kv(x - x2) + kv(y - y2));
		this->Sosed_down->Get_Center(x2, y2);
		//cout << "Gran.cpp    " << x2 << " " << y2 << endl;
		dist1 = sqrt(kv(x - x2) + kv(y - y2));

		par = par2;

		par.v = linear(-dist1, par1.v, -dist2, par2.v, dist2, -par2.v, 0.0);

		par.v_H1 = linear(-dist1, par1.v_H1, -dist2, par2.v_H1, dist2, -par2.v_H1, 0.0);
		par.v_H2 = linear(-dist1, par1.v_H2, -dist2, par2.v_H2, dist2, -par2.v_H2, 0.0);
		par.v_H3 = linear(-dist1, par1.v_H3, -dist2, par2.v_H3, dist2, -par2.v_H3, 0.0);
		par.v_H4 = linear(-dist1, par1.v_H4, -dist2, par2.v_H4, dist2, -par2.v_H4, 0.0);

		//par.ro = linear(-dist1, par1.ro, -dist2, par2.ro, 0.0);
		//par.p = linear(-dist1, par1.p, -dist2, par2.p,  0.0);


		/*par.v = -par.v;
		par.v_H1 = -par.v_H1;
		par.v_H2 = -par.v_H2;
		par.v_H3 = -par.v_H3;
		par.v_H4 = -par.v_H4;*/
	}
	else if (this->type == Inner_sphere)
	{
		auto par2 = this->Master->par[i];   // this->Sosed->par[i];
		double x, y;
		this->Get_Center(x, y);
		double dist = sqrt(x * x + y * y);
		double ro = 1.0 / (kv(chi_real) * dist * dist);
		double P_E = ro * chi_real * chi_real / (ggg * 10.0 * 10.0);
		double T_p = (P_E * pow(1.0 / dist, 2.0 * ggg)) / (2.0 * ro / (dist * dist));  

		par = { ro, P_E, chi_real * x / dist, chi_real * y / dist, ro,//
		0.0000001, (0.0000001 * chi_real * chi_real / (ggg * 14.1344 * 14.1344)), chi_real * x / dist, chi_real * y / dist,//
			par2.ro_H2, par2.p_H2, par2.u_H2, par2.v_H2 , //
			par2.ro_H3, par2.p_H3, par2.u_H3, par2.v_H3,//
			par2.ro_H4, par2.p_H4, par2.u_H4, par2.v_H4 };

		//par = { ro * kv(chi_real / chi_), P_E, chi_real * x / dist / (chi_real / chi_), chi_real * y / dist / (chi_real / chi_), ro * (chi_real / chi_),//
		//0.0000001, (0.0000001 * chi_real * chi_real / (ggg * 14.1344 * 14.1344)), chi_real * x / dist, chi_real * y / dist,//
		//	par2.ro_H2, par2.p_H2, par2.u_H2, par2.v_H2 , //
		//	par2.ro_H3, par2.p_H3, par2.u_H3, par2.v_H3,//
		//	par2.ro_H4, par2.p_H4, par2.u_H4, par2.v_H4 };


		//cout << "Gran.cpp    " << par.ro << " " << dist << endl;
	}
	else if (this->type == Input)
	{
		auto par2 = this->Master->par[i];
		par = { 1.0, 1.0, Velosity_inf, 0.0, 100.0, par2.ro_H1, par2.p_H1, par2.u_H1, par2.v_H1,//
		par2.ro_H2, par2.p_H2, par2.u_H2, par2.v_H2 , par2.ro_H3, par2.p_H3, par2.u_H3, par2.v_H3,//
		1.0, 0.5, Velosity_inf, 0.0 };

		//par = { 4.0, 2.5, Velosity_inf, 0.0, 100.0, par2.ro_H1, par2.p_H1, par2.u_H1, par2.v_H1,//
		//par2.ro_H2, par2.p_H2, par2.u_H2, par2.v_H2 , par2.ro_H3, par2.p_H3, par2.u_H3, par2.v_H3,//
		//1.0, 0.5, Velosity_inf, 0.0 };

	}

	return;
}


void Gran::Get_par_TVD_radial(Parametr& par, int i)  
// Здесь задаются граничные условия
// Сносит значение в текущей ячейке на текущую гарнь
// 
{
	if (this->type == Usualy)  // Не надо менять
	{
		if (this->Sosed_down == nullptr)
		{
				bool b = false;
				for (auto& i : this->Master->Grans)
				{
					if (i->type == Axis)
					{
						b = true;
						break;
					}
				}
				if (b == true)
				{
					double x, y, x2, y2, dist1, dist2, dist3;
					this->Get_Center(x, y);
					double phi = polar_angle(x, y);
					double r = sqrt(kv(x) + kv(y));
					auto par2 = this->Master->par[i];
					auto par1 = this->Master->par[i];
					par1.v = -par1.v;
					par1.v_H1 = -par1.v_H1;
					par1.v_H2 = -par1.v_H2;
					par1.v_H3 = -par1.v_H3;
					par1.v_H4 = -par1.v_H4;
					auto par3 = this->Sosed->par[i];
					this->Master->Get_Center(x2, y2);
					double r2 = sqrt(kv(x2) + kv(y2));
					double phi2 = polar_angle(x2, y2);
					dist2 = sqrt(kv(x - x2) + kv(y - y2));
					this->Sosed->Get_Center(x2, y2);
					double r3 = sqrt(kv(x2) + kv(y2));
					dist3 = sqrt(kv(x - x2) + kv(y - y2));
					double phi3 = polar_angle(x2, y2);
					this->Master->Get_Center(x2, y2);
					y2 = -y2;
					double phi1 = polar_angle(x2, y2);
					double r1 = sqrt(kv(x2) + kv(y2));
					dist1 = sqrt(kv(x - x2) + kv(y - y2));

					double fr1 = par1.u_H1 * cos(phi1) + par1.v_H1 * sin(phi1);
					double ff1 = -par1.u_H1 * sin(phi1) + par1.v_H1 * cos(phi1);

					double fr2 = par2.u_H1 * cos(phi2) + par2.v_H1 * sin(phi2);
					double ff2 = -par2.u_H1 * sin(phi2) + par2.v_H1 * cos(phi2);

					double fr3 = par3.u_H1 * cos(phi3) + par3.v_H1 * sin(phi3);
					double ff3 = -par3.u_H1 * sin(phi3) + par3.v_H1 * cos(phi3);


					double fr = linear(-dist1, fr1, -dist2, fr2, dist3, fr3, 0.0);
					par.ro = linear(-dist1, par1.ro * kv(r1), -dist2, par2.ro * kv(r2), dist3, par3.ro * kv(r3), 0.0) / kv(r);
					par.p = linear(-dist1, par1.p * pow(r1, 2.0 * ggg), -dist2, par2.p * pow(r2, 2.0 * ggg), dist3, par3.p * pow(r3, 2.0 * ggg), 0.0) / pow(r, 2.0 * ggg);
					par.p_H1 = linear(-dist1, par1.p_H1, -dist2, par2.p_H1, dist3, par3.p_H1, 0.0);
					par.ro_H1 = linear(-dist1, par1.ro_H1 * pow(r1, H_pow), -dist2, par2.ro_H1 * pow(r2, H_pow), dist3, par3.ro_H1 * pow(r3, H_pow), 0.0) / pow(r, H_pow);
					//par.ro_H1 = linear(-dist1, par1.ro_H1, -dist2, par2.ro_H1, dist3, par3.ro_H1, 0.0);
					par.u_H1 = fr * cos(phi) - ff2 * sin(phi);
					par.v_H1 = fr * sin(phi) + ff2 * cos(phi);
				}
				else
				{
					par = this->Master->par[i];
				}
				return;
			
		}
		double x, y, x2, y2, dist1, dist2, dist3;
		this->Get_Center(x, y);
		double phi = polar_angle(x, y);
		double r = sqrt(kv(x) + kv(y));
		auto par2 = this->Master->par[i];
		auto par1 = this->Sosed_down->par[i];
		auto par3 = this->Sosed->par[i];
		this->Master->Get_Center(x2, y2);
		double r2 = sqrt(kv(x2) + kv(y2));
		double phi2 = polar_angle(x2, y2);
		dist2 = sqrt(kv(x - x2) + kv(y - y2));
		this->Sosed->Get_Center(x2, y2);
		double r3 = sqrt(kv(x2) + kv(y2));
		dist3 = sqrt(kv(x - x2) + kv(y - y2));
		double phi3 = polar_angle(x2, y2);
		this->Sosed_down->Get_Center(x2, y2);
		double phi1 = polar_angle(x2, y2);
		double r1 = sqrt(kv(x2) + kv(y2));
		dist1 = sqrt(kv(x - x2) + kv(y - y2));

		double fr1 = par1.u_H1 * cos(phi1) + par1.v_H1 * sin(phi1);
		double ff1 = -par1.u_H1 * sin(phi1) + par1.v_H1 * cos(phi1);

		double fr2 = par2.u_H1 * cos(phi2) + par2.v_H1 * sin(phi2);
		double ff2 = -par2.u_H1 * sin(phi2) + par2.v_H1 * cos(phi2);

		double fr3 = par3.u_H1 * cos(phi3) + par3.v_H1 * sin(phi3);
		double ff3 = -par3.u_H1 * sin(phi3) + par3.v_H1 * cos(phi3);

		
		double fr = linear(-dist1, fr1, -dist2, fr2, dist3, fr3, 0.0);
		par.ro = linear(-dist1, par1.ro * kv(r1), -dist2, par2.ro * kv(r2), dist3, par3.ro * kv(r3), 0.0) / kv(r);
		par.p = linear(-dist1, par1.p * pow(r1, 2.0 * ggg), -dist2, par2.p * pow(r2, 2.0 * ggg), dist3, par3.p * pow(r3, 2.0 * ggg), 0.0) / pow(r, 2.0 * ggg);
		par.p_H1 = linear(-dist1, par1.p_H1, -dist2, par2.p_H1, dist3, par3.p_H1, 0.0);
		par.ro_H1 = linear(-dist1, par1.ro_H1 * pow(r1, H_pow), -dist2, par2.ro_H1 * pow(r2, H_pow), dist3, par3.ro_H1 * pow(r3, H_pow), 0.0) / pow(r, H_pow);
		//par.ro_H1 = linear(-dist1, par1.ro_H1, -dist2, par2.ro_H1, dist3, par3.ro_H1, 0.0);
		par.u_H1 = fr * cos(phi) - ff2 * sin(phi);
		par.v_H1 = fr * sin(phi) + ff2 * cos(phi);
	}

	return;
}