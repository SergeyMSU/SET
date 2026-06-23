#include "Interpol_Setka.h"

using namespace std;

Interpol_Setka::Interpol_Setka(Setka* SS)
{
	double x, y;
	int n1, n2, n3, nn;
	nn = 0;
	for (int i = 0; i < SS->All_Cells.size(); i++)
	{
		SS->All_Cells[i]->Get_Center(x, y);
		auto K = new Point(x, y);
		this->All_Points.push_back(K);
		K->number = nn;
		nn++;
	}

	for (int i = 0; i < SS->M1 + SS->M2 + SS->M3 + SS->M4; i++)
	{
		SS->All_Cells[i]->Get_Center(x, y);
		auto K = new Point(x, 0.0);
		this->All_Points.push_back(K);
		K->number = nn;
		nn++;
	}

	for (int i = 0; i < SS->All_Cells.size(); i++)
	{
		for (auto& j : SS->All_Cells[i]->contour)
		{
			if (j->my_cell.size() == 4)
			{
				n1 = j->my_cell[0]->number;
				n2 = j->my_cell[1]->number;
				n3 = j->my_cell[2]->number;
				auto C = new Cell(this->All_Points[n1], this->All_Points[n2], this->All_Points[n3]);
				this->All_Cells.push_back(C);

				n1 = j->my_cell[2]->number;
				n2 = j->my_cell[1]->number;
				n3 = j->my_cell[3]->number;
				C = new Cell(this->All_Points[n1], this->All_Points[n2], this->All_Points[n3]);
				this->All_Cells.push_back(C);
			}
		}
	}


}

void Interpol_Setka::Print_Cell(string name)
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
