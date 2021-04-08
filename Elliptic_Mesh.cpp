#include <iostream>
#include <fstream>
#include <vector>
#include<iomanip>
using namespace std;
#include "Elliptic_Mesh.h"

int main()
{
	initialize();

	read_boundary_points();

	relaxation_method();

	output_grid();
}

void output_grid()
{
	fstream file;
	file.open("naca0012/grid.p3d", ios_base::out);
	
	//输出plot3d格式的网格
	int numberOfGridZones = 1;
	file << numberOfGridZones << "\n";
	file << NI << "\t";
	file << NJ << "\n";
	file << setiosflags(ios::right);
	file << setiosflags(ios::scientific);
	file << setprecision(15);
	//for (int i = 0; i < NI; ++i)
	//{
	//	for (int j = 0; j < NJ; ++j)
	//	{
	//		file << globalCoordX[i][j] << "\t";
	//	}		
	//}

	for (int j = 0; j < NJ; ++j)
	{
		for (int i = 0; i < NI; ++i)
		{
			file << globalCoordX[i][j] << "\t";
		}
	}
	file << "\n";

	//for (int i = 0; i < NI; ++i)
	//{
	//	for (int j = 0; j < NJ; ++j)
	//	{
	//		file << globalCoordY[i][j] << "\t";
	//	}
	//}
	for (int j = 0; j < NJ; ++j)
	{
		for (int i = 0; i < NI; ++i)
		{
			file << globalCoordY[i][j] << "\t";
		}
	}
	file << "\n";
	file.close();
	cout << "grid generation completed, program ends!" << endl;
}

void relaxation_method()
{
	//solve x	
	double error = 1e30; 
	int iter = 0;
	double omega = 1.1;
	do
	{
		double errorL1 = 0;
		for (int i = 1; i < NI - 1; ++i)
		{
			for (int j = 1; j < NJ - 1; ++j)
			{
				double oldX = globalCoordX[i][j];
				double oldY = globalCoordY[i][j];

				double xdiff_j = globalCoordX[i][j + 1] - globalCoordX[i][j - 1];
				double ydiff_j = globalCoordY[i][j + 1] - globalCoordY[i][j - 1];

				double xdiff_i = globalCoordX[i + 1][j] - globalCoordX[i - 1][j];
				double ydiff_i = globalCoordY[i + 1][j] - globalCoordY[i - 1][j];

				double hybrid_diff_x = globalCoordX[i + 1][j + 1] - globalCoordX[i - 1][j + 1]
								     - globalCoordX[i + 1][j - 1] + globalCoordX[i - 1][j - 1];
				double hybrid_diff_y = globalCoordY[i + 1][j + 1] - globalCoordY[i - 1][j + 1]
								     - globalCoordY[i + 1][j - 1] + globalCoordY[i - 1][j - 1];

				double Ae = (xdiff_j * xdiff_j + ydiff_j * ydiff_j) / 4.0;
				double An = (xdiff_i * xdiff_i + ydiff_i * ydiff_i) / 4.0;
				double Aw = Ae;
				double As = An;
				double Ac = -(Ae + Aw + An + As) + 1e-40;
				double Sij_x = (xdiff_i * xdiff_j + ydiff_i * ydiff_j) * hybrid_diff_x / 8.0;
				double Sij_y = (xdiff_i * xdiff_j + ydiff_i * ydiff_j) * hybrid_diff_y / 8.0;

				//globalCoordX[i][j] = 1.0 / Ac * (Sij_x - Ae * globalCoordX[i + 1][j] - Aw * globalCoordX[i - 1][j]
				//									 - An * globalCoordX[i][j + 1] - As * globalCoordX[i][j - 1]);
				//globalCoordY[i][j] = 1.0 / Ac * (Sij_y - Ae * globalCoordY[i + 1][j] - Aw * globalCoordY[i - 1][j]
				//									   - An * globalCoordY[i][j + 1] - As * globalCoordY[i][j - 1]);

				globalCoordX[i][j] = ( 1.0 - omega ) * oldX + omega * (
									 1.0 / Ac * (Sij_x - Ae * globalCoordX[i + 1][j] - Aw * globalCoordX[i - 1][j]
													   - An * globalCoordX[i][j + 1] - As * globalCoordX[i][j - 1]) );
				globalCoordY[i][j] = ( 1.0 - omega ) * oldY + omega * (
									 1.0 / Ac * (Sij_y - Ae * globalCoordY[i + 1][j] - Aw * globalCoordY[i - 1][j]
													   - An * globalCoordY[i][j + 1] - As * globalCoordY[i][j - 1]) );
				
				double errorL1_x = abs(globalCoordX[i][j] - oldX);
				double errorL1_y = abs(globalCoordY[i][j] - oldY);

				errorL1 += errorL1_x;
				errorL1 += errorL1_y;

				if (isnan(errorL1))
				{
					int kkk = 1;
				}				
			}
		}

		errorL1 /=  (2.0 * NI * NJ);

		iter++;
		error = errorL1;
		if (iter % 200 == 0)
		{
			cout << "iter = " << iter << "\terror = " << error << endl;
		}		
	} while (error > 1e-13);

	cout << "Coordinates solved. " << endl;
}

void initialize()
{
	NI = 201;
	NJ = 71;

	//数据是按列存储的,同一ksi存在同一数组中
	globalCoordX.resize(NI);
	globalCoordY.resize(NI);
	for (int i = 0; i < NI; ++i)
	{
		globalCoordX[i].resize(NJ);
		globalCoordY[i].resize(NJ);
	}
}

void read_boundary_points()
{
	//数据是按列存储的
	read_BC("naca0012/BC1.x", globalCoordX, globalCoordY, "bottom");
	read_BC("naca0012/BC2.x", globalCoordX, globalCoordY, "up");
	read_BC("naca0012/BC3.x", globalCoordX, globalCoordY, "right");
	read_BC("naca0012/BC3.x", globalCoordX, globalCoordY, "left");
}

void read_BC(string fileName, vector < vector<double> >& globalCoordX, vector < vector<double> >& globalCoordY, string pos)
{
	fstream file;
	file.open(fileName, ios_base::in);
	int numberOfZones;
	file >> numberOfZones;

	vector<int> ni(numberOfZones);
	vector<int> nj(numberOfZones);
	vector<int> numberOfPoints(numberOfZones);
	int totalNumberOfPoints = 0;
	for (int iZone = 0; iZone < numberOfZones; ++iZone)
	{
		file >> ni[iZone];
		file >> nj[iZone];
		numberOfPoints[iZone] = ni[iZone] * nj[iZone];
		totalNumberOfPoints += numberOfPoints[iZone];
	}

	vector<double> xCoord(totalNumberOfPoints);
	vector<double> yCoord(totalNumberOfPoints);
	int pointCount = 0;
	for (int iZone = 0; iZone < numberOfZones; ++iZone)
	{
		for (int iPoint = 0; iPoint < numberOfPoints[iZone]; ++iPoint)
		{
			file >> xCoord[iPoint + pointCount];
		}
		for (int iPoint = 0; iPoint < numberOfPoints[iZone]; ++iPoint)
		{
			file >> yCoord[iPoint + pointCount];
		}
		pointCount += numberOfPoints[iZone];
	}

	if (pos == "bottom")
	{
		for (int i = 0; i < NI; ++i)
		{
			globalCoordX[i][0] = xCoord[i];
			globalCoordY[i][0] = yCoord[i];
		}		
	}
	else if (pos == "up")
	{
		for (int i = 0; i < NI; ++i)
		{
			globalCoordX[i][NJ - 1] = xCoord[i];
			globalCoordY[i][NJ - 1] = yCoord[i];
		}
	}
	else if (pos == "left")
	{
		for (int j = 0; j < NJ; ++j)
		{
			globalCoordX[0][j] = xCoord[j];
			globalCoordY[0][j] = yCoord[j];
		}
	}
	else if (pos == "right")
	{
		for (int j = 0; j < NJ; ++j)
		{
			globalCoordX[NI-1][j] = xCoord[j];
			globalCoordY[NI-1][j] = yCoord[j];
		}
	}

	int kkk = 1;
}