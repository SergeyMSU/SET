#include "Dist_func.h"
#include "Help.h"
#include <filesystem> 



Dist_func::Dist_func(const unsigned short& a1, const unsigned short& a2, const unsigned short& a3, //
	const double& c1, const double& d1, //
	const double& c2, const double& d2, //
	const double& c3, const double& d3)
{
	this->n1 = a1;
	this->n2 = a2;
	this->n3 = a3;
	this->name = "null_name";

	this->V = new double** [this->n1];
	for (int i = 0; i < this->n1; i++)
	{
		V[i] = new double* [this->n2];
	}

	for (int i = 0; i < this->n1; i++)
	{
		for (int j = 0; j < this->n2; j++)
		{
			V[i][j] = new double[this->n3];
		}
	}

	for (int i = 0; i < this->n1; i++)
	{
		for (int j = 0; j < this->n2; j++)
		{
			for (int k = 0; k < this->n3; k++)
			{
				V[i][j][k] = 0.0;
			}
		}
	}

	this->c1 = c1;
	this->d1 = d1;

	this->c2 = c2;
	this->d2 = d2;

	this->c3 = c3;
	this->d3 = d3;

	if ((a1 % 2 != 0) || (a2 % 2 != 0) || (a3 % 2 != 0))
	{
		std::cout << "ERROR nreufvyv34753247  27637424  " << std::endl;
	}

	this->a1 = pow(this->d1 + 1.0, 2.0 / this->n1);
	this->a2 = pow(this->d2 + 1.0, 2.0 / this->n2);
	this->a3 = pow(this->d3 + 1.0, 2.0 / this->n3);
}

void Dist_func::v_cyl_to_v_xyz(double v_rho, double v_phi, double x, double y, double& vx, double& vy) {
	double phi = atan2(y, x);
	vx = cos(phi) * v_rho - sin(phi) * v_phi;
	vy = sin(phi) * v_rho + cos(phi) * v_phi;

}
void Dist_func::v_cyl(double vx, double vy, double x, double y, double& w_rho, double& w_phi) {
	double phi = atan2(y, x); //polar_angle(x, y); // atan2(y, x);
	w_rho = cos(phi) * vx + sin(phi) * vy;
	w_phi = -sin(phi) * vx + cos(phi) * vy;
}


bool Dist_func::call_name(std::string name)
{
	this->name = name;
	return true;
}

bool Dist_func::Add_point(const double& V1xyz, const double& V2xyz, const double& V3xyz, const double& y, const double& z, const double& mu)
//скорости и координаты в твоей СК!
{
	double V1 = V1xyz;
	double V2 = 0;
	double V3 = 0;
	v_cyl(V2xyz, V3xyz, y, z, V2, V3);

	// Функция добавляет точку в массив
	//std::cout << "w_cyl " << V1 << " "<< V2 << " " << V3 << std::endl;

	//сначала нужно определить в какую ячейку массива записать число.
	unsigned short i = 0, j = 0, k = 0;

	// double a = pow(this->d1 + 1.0 - this->c1, 2.0 / this->n1);

	i = (int)(log(fabs(V1 - this->c1) + 1.0) / log(this->a1));
	if (i >= this->n1 / 2)
	{
		i = this->n1 / 2 - 1;
	}
	if (V1 - this->c1 < 0.0)
	{
		i = i + this->n1 / 2;
	}

	//cout << this->a1 << " " << i << endl;

	j = (int)(log(fabs(V2 - this->c2) + 1.0) / log(this->a2));
	if (j >= this->n2 / 2)
	{
		j = this->n2 / 2 - 1;
	}
	if (V2 - this->c2 < 0.0)
	{
		j = j + this->n2 / 2;
	}

	k = (int)(log(fabs(V3 - this->c3) + 1.0) / log(this->a3));
	if (k >= this->n3 / 2)
	{
		k = this->n3 / 2 - 1;
	}
	if (V3 - this->c3 < 0.0)
	{
		k = k + this->n3 / 2;
	}
	V[i][j][k] += mu;
	return true;
}

bool Dist_func::normir(const double& ccc)
{
	for (int i = 0; i < this->n1; i++)
	{
		for (int j = 0; j < this->n2; j++)
		{
			for (int k = 0; k < this->n3; k++)
			{
				this->V[i][j][k] = this->V[i][j][k] * ccc;
			}
		}
	}

	// Нормировка из-за неравномерности функции распределения

	for (int i = 0; i < this->n1; i++)
	{
		double l1 = this->c1 - 1.0 + pow(this->a1, 1.0 * i);
		double r1 = this->c1 - 1.0 + pow(this->a1, 1.0 * i + 1.0);
		if (i >= this->n1 / 2)
		{
			l1 = -(-this->c1 - 1.0 + pow(this->a1, 1.0 * i - this->n1 / 2.0));
			r1 = -(-this->c1 - 1.0 + pow(this->a1, 1.0 * i + 1.0 - this->n1 / 2.0));
		}

		for (int j = 0; j < this->n2; j++)
		{
			double l2 = this->c2 - 1.0 + pow(this->a2, 1.0 * j);
			double r2 = this->c2 - 1.0 + pow(this->a2, 1.0 * j + 1.0);
			if (j >= this->n2 / 2)
			{
				l2 = -(-this->c2 - 1.0 + pow(this->a2, 1.0 * j - this->n2 / 2.0));
				r2 = -(-this->c2 - 1.0 + pow(this->a2, 1.0 * j + 1.0 - this->n2 / 2.0));
			}

			for (int k = 0; k < this->n3; k++)
			{
				double l3 = this->c3 - 1.0 + pow(this->a3, 1.0 * k);
				double r3 = this->c3 - 1.0 + pow(this->a3, 1.0 * k + 1.0);
				if (k >= this->n3 / 2)
				{
					l3 = -(-this->c3 - 1.0 + pow(this->a3, 1.0 * k - this->n3 / 2.0));
					r3 = -(-this->c3 - 1.0 + pow(this->a3, 1.0 * k + 1.0 - this->n3 / 2.0));
				}

				//double S = fabs((r1 - l1) * (r2 * r2 - l2 * l2) * (r3 - l3) / 2.0);
				double S = fabs((r1 - l1) * (r2 - l2) * (r3 - l3));
				this->V[i][j][k] = this->V[i][j][k] / S;
			}
		}
	}


	return true;
}

bool Dist_func::print_1d(int koord)
{
	std::ofstream fout;
	std::string name_f = std::to_string(koord) + "_Dist_func_" + this->name + ".txt";
	fout.open(name_f);

	fout << "TITLE = \"HP\"  VARIABLES = \"v\", \"f\"," << "ZONE T = \"HP\"" << std::endl;

	if (koord == 1)
	{
		for (int i = 0; i < this->n1 / 2; i++)
		{
			double l = this->c1 - 1.0 + pow(this->a1, 1.0 * i);
			double r = this->c1 - 1.0 + pow(this->a1, 1.0 * i + 1);
			double S = 0.0;

			for (int j = 0; j < this->n2; j++)
			{
				for (int k = 0; k < this->n3; k++)
				{
					S = S + V[i][j][k];
				}
			}

			fout << 0.5 * (l + r) << " " << S << std::endl;
		}

		for (int i = 0 + this->n1 / 2; i < this->n1; i++)
		{
			double l = -1.0 + pow(this->a1, 1.0 * i - this->n1 / 2.0);
			double r = -1.0 + pow(this->a1, 1.0 * i + 1.0 - this->n1 / 2.0);
			double S = 0.0;

			for (int j = 0; j < this->n2; j++)
			{
				for (int k = 0; k < this->n3; k++)
				{
					S = S + V[i][j][k];
				}
			}

			fout << this->c1 - 0.5 * (l + r) << " " << S << std::endl;
		}
	}
	else if (koord == 2)
	{
		for (int j = 0; j < this->n2 / 2; j++)
		{
			double l = this->c2 - 1.0 + pow(this->a2, 1.0 * j);
			double r = this->c2 - 1.0 + pow(this->a2, 1.0 * j + 1.0);
			double S = 0.0;

			for (int i = 0; i < this->n1; i++)
			{
				for (int k = 0; k < this->n3; k++)
				{
					S = S + V[i][j][k];
				}
			}

			fout << 0.5 * (l + r) << " " << S << std::endl;
		}

		for (int j = 0 + this->n2 / 2; j < this->n2; j++)
		{
			double l = -1.0 + pow(this->a2, 1.0 * j - this->n2 / 2.0);
			double r = -1.0 + pow(this->a2, 1.0 * j + 1.0 - this->n2 / 2.0);
			double S = 0.0;

			for (int i = 0; i < this->n1; i++)
			{
				for (int k = 0; k < this->n3; k++)
				{
					S = S + V[i][j][k];
				}
			}

			fout << this->c2 - 0.5 * (l + r) << " " << S << std::endl;
		}
	}
	else if (koord == 3)
	{
		for (int k = 0; k < this->n3 / 2.0; k++)
		{
			double l = this->c3 - 1.0 + pow(this->a3, 1.0 * k);
			double r = this->c3 - 1.0 + pow(this->a3, 1.0 * k + 1.0);
			double S = 0.0;

			for (int i = 0; i < this->n1; i++)
			{
				for (int j = 0; j < this->n2; j++)
				{
					S = S + V[i][j][k];
				}
			}

			fout << 0.5 * (l + r) << " " << S << std::endl;
		}

		for (int k = 0 + this->n3 / 2; k < this->n3; k++)
		{
			double l = -1.0 + pow(this->a3, 1.0 * k - this->n3 / 2.0);
			double r = -1.0 + pow(this->a3, 1.0 * k + 1.0 - this->n3 / 2.0);
			double S = 0.0;

			for (int i = 0; i < this->n1; i++)
			{
				for (int j = 0; j < this->n2; j++)
				{
					S = S + V[i][j][k];
				}
			}

			fout << this->c3 - 0.5 * (l + r) << " " << S << std::endl;
		}
	}

	return true;
}

bool Dist_func::print_3d(void)
{
	std::ofstream fout;
	std::string name_f = "func/3D_Dist_func_" + this->name + ".txt";

	std::filesystem::path file_path = name_f;

	// Проверяем существование и доступность директории
	if (!std::filesystem::exists(file_path.parent_path())) 
	{
		cout << "ERROR iwjfheuit4h9834ut893g" << endl;
		exit(-1);
	}

	fout.open(name_f);

	if (!fout.is_open()) 
	{
		cout << "Error ugewerfwefhgoiue " << endl;
		exit(-1);
	}

	fout << this->xxx << " " << this->yyy << endl;
	fout << "TITLE = \"HP\"  VARIABLES = \"vz\", \"vr\", \"vphi\", \"f\"," << "ZONE T = \"HP\"" << std::endl;


	for (int i = 0; i < this->n1; i++)
	{

		double l1;
		double r1;
		if (i >= this->n1 / 2)
		{
			l1 = -(-this->c1 - 1.0 + pow(this->a1, 1.0 * i - this->n1 / 2.0));
			r1 = -(-this->c1 - 1.0 + pow(this->a1, 1.0 * i + 1.0 - this->n1 / 2.0));
		}
		else {
			l1 = (this->c1 - 1.0 + pow(this->a1, 1.0 * i));
			r1 = (this->c1 - 1.0 + pow(this->a1, 1.0 * i + 1.0));
		}

		for (int j = 0; j < this->n2; j++)
		{
			double l2 = this->c2 - 1.0 + pow(this->a2, 1.0 * j);
			double r2 = this->c2 - 1.0 + pow(this->a2, 1.0 * j + 1.0);
			if (j >= this->n2 / 2)
			{
				l2 = -(-this->c2 - 1.0 + pow(this->a2, 1.0 * j - this->n2 / 2.0));
				r2 = -(-this->c2 - 1.0 + pow(this->a2, 1.0 * j + 1.0 - this->n2 / 2.0));
			}

			for (int k = 0; k < this->n3; k++)
			{
				double l3 = this->c3 - 1.0 + pow(this->a3, 1.0 * k);
				double r3 = this->c3 - 1.0 + pow(this->a3, 1.0 * k + 1.0);
				if (k >= this->n3 / 2)
				{
					l3 = -(-this->c3 - 1.0 + pow(this->a3, 1.0 * k - this->n3 / 2.0));
					r3 = -(-this->c3 - 1.0 + pow(this->a3, 1.0 * k + 1.0 - this->n3 / 2.0));
				}

				fout << 0.5 * (l1 + r1) << " " << 0.5 * (l2 + r2) << " " << 0.5 * (l3 + r3) << " " << this->V[i][j][k] << std::endl;
			}
		}
	}


	return true;
}
//-----------------------------------ДЛЯ ОТЛАДКИ-------------------------------------------------------
/*int main()
{
	Dist_func test(20, 20, 20, //
		-1.5, 6, //
		0.0, 6, //
		0.0, 6);
	test.Add_point(0, 0.1, 0, 0.0, 3.14/2.0, 3536);
	//test.print_3d();
	double x = 1.0;
	double y = 0.0;
	double vx = 0.8;
	double vy = 0;
	double vz = 0.0;
	double w_rho = 0;
	double w_phi = 0;
	double v1 = 0;
	double v2 = 0;

	test.v_cyl(vx,vy,x,y,w_rho, w_phi);
	test.v_cyl_to_v_xyz(w_rho, w_phi, x, y, v1, v2);
	std::cout << vx << " " << vy << " " << vz << std::endl;
	std::cout << w_rho << " " << w_phi << std::endl;
	std::cout << v1 << " " << v2 << " " << vz << std::endl;
	return 0;
}*/