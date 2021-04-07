#include <iostream>
#include <fstream>
#include <vector>
using namespace std;
#include "Elliptic_Mesh.h"

int main()
{
	initialize();
	read_boundary_points();
}

void initialize()
{
	NI = 202;
	NJ = 71;

	globalCoordX.resize(NI);
	for (int i = 0; i < NI; ++i)
	{
		globalCoordX[i].resize(NJ);
	}

	globalCoordY.resize(NI);
	for (int i = 0; i < NI; ++i)
	{
		globalCoordY[i].resize(NJ);
	}
}

void read_boundary_points()
{
	read_BC("naca0012/BC1.x", globalCoordX, globalCoordY, "bottom");
	read_BC("naca0012/BC2.x", globalCoordX, globalCoordY, "up");
	read_BC("naca0012/BC3.x", globalCoordX, globalCoordY, "left");
	read_BC("naca0012/BC3.x", globalCoordX, globalCoordY, "right");
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