#include "Rail.h"
#include <math.h>
using namespace std;

Rail::Rail(const double& s)
{
	this->s = s;
	this->M1 = 0;
	this->M2 = 0;
	this->M3 = 0;
	this->M4 = 0;
	this->type = Rail_type::No;
}

void Rail::Init_start(Rail* T)
{
	double r, l;

	if (this->type == A)
	{
		double x = pow(R2_ / R1_, 1.0 / (this->M1 - 1));
		for (int i = 0; i < this->M1; i++)
		{
			r = R1_ * pow(x, i);  
			//r = R1_ + (R2_ - R1_) * i / (this->M1 - 1);
			auto K = new Point(r * cos(this->s), r * sin(this->s));
			this->All_point.push_back(K);
			K->type = P_U_1;
			if (i == this->M1 - 1)
			{
				this->Key_point.push_back(K);
				K->type = P_Inner_shock;
			}
			if (i == 0)
			{
				K->type = P_Inner_Boandary;
			}
		}

		for (int i = 0; i < this->M2; i++)
		{
			r = R2_ + (R3_ - R2_) * (i + 1) / (this->M2);
			auto K = new Point(r * cos(this->s), r * sin(this->s));
			this->All_point.push_back(K);
			K->type = P_U_3;
			if (i == this->M2 - 1)
			{
				this->Key_point.push_back(K);
			}
			if (i == this->M2 - 1)
			{
				K->type = P_Contact;
			}
		}

		for (int i = 0; i < this->M3; i++)
		{
			r = R3_ + (R4_ - R3_) * (i + 1) / (this->M3);
			auto K = new Point(r * cos(this->s), r * sin(this->s));
			this->All_point.push_back(K);
			K->type = P_U_4;
			if (i == this->M3 - 1)
			{
				this->Key_point.push_back(K);
				K->type = P_Outer_shock;
			}
		}

		x = pow(R5_ / R4_, 1.0 / (this->M4));
		for (int i = 0; i < this->M4; i++)
		{
			r = R4_ * pow(x, i + 1);  //R4_ + (R5_ - R4_) * (i + 1) / (this->M4);
			auto K = new Point(r * cos(this->s), r * sin(this->s));
			this->All_point.push_back(K);
			K->type = P_U_5;
			if (i == this->M4 - 1)
			{
				this->Key_point.push_back(K);
				K->type = P_Outer_Boandary;
			}
		}
	}
	else if (this->type == B)
	{
		double x = pow(R2_ / R1_, 1.0 / (this->M1 - 1));
		for (int i = 0; i < this->M1; i++)
		{
			r = R1_ * pow(x, i); //
			//r = R1_ + (R2_ - R1_) * i / (this->M1 - 1);
			auto K = new Point(r * cos(this->s), r * sin(this->s));
			this->All_point.push_back(K);
			K->type = P_U_1;
			if (i == this->M1 - 1)
			{
				this->Key_point.push_back(K);
				K->type = P_Inner_shock;
			}
			if (i == 0)
			{
				K->type = P_Inner_Boandary;
			}
		}

		for (int i = 0; i < this->M2; i++)
		{
			l = R2_ * sin(this->s) + (i + 1) * (R3_ - R2_ * sin(this->s))/(this->M2);
			auto K = new Point(R2_ * cos(this->s), l);
			this->All_point.push_back(K);
			K->type = P_U_3;
			if (i == this->M2 - 1)
			{
				this->Key_point.push_back(K);
			}
			if (i == this->M2 - 1)
			{
				K->type = P_Contact;
			}
		}

		for (int i = 0; i < this->M3; i++)
		{
			l = R3_ + (i + 1) * (R4_ - R3_) /this->M3;
			auto K = new Point(R2_ * cos(this->s), l);
			this->All_point.push_back(K);
			K->type = P_U_4;
			if (i == this->M3 - 1)
			{
				this->Key_point.push_back(K);
				K->type = P_Outer_shock;
			}
		}

		x = pow(R5_ / R4_, 1.0 / (this->M4));
		for (int i = 0; i < this->M4; i++)
		{
			l = R4_ * pow(x, i + 1);  //R4_ + (i + 1) * (R5_ - R4_) / this->M4;
			auto K = new Point(R2_ * cos(this->s), l);
			this->All_point.push_back(K);
			K->type = P_U_5;
			if (i == this->M4 - 1)
			{
				this->Key_point.push_back(K);
				K->type = P_Outer_Boandary;
			}
		}
	}
	else if (this->type == C)
	{
		double x = pow(R2_ / R1_, 1.0 / (this->M1 - 1));
		for (int i = 0; i < this->M1; i++)
		{
			r = R1_ * pow(x, i); // 
			//r = R1_ + (R2_ - R1_) * i / (this->M1 - 1);
			auto K = new Point(r * cos(this->s), r * sin(this->s));
			this->All_point.push_back(K);
			K->type = P_U_1;
			if (i == this->M1 - 1)
			{
				this->Key_point.push_back(K);
				K->type = P_Inner_shock;
			}
			if (i == 0)
			{
				K->type = P_Inner_Boandary;
			}
		}

		x = pow(Left_ / (R2_ * cos(this->s)), 1.0 / (this->M2));
		for (int i = 0; i < this->M2; i++)
		{
			l = R2_ * cos(this->s) * pow(x, i + 1); //R2_ * cos(this->s) - (i + 1) * (R2_ * cos(this->s) - Left_) / this->M2;
			auto K = new Point(l, R2_ * sin(this->s));
			this->All_point.push_back(K);
			K->type = P_U_2;
			if (i == this->M2 - 1)
			{
				this->Key_point.push_back(K);
				K->type = P_Outer_Boandary;
			}
		}

	}
}