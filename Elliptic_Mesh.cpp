#include <iostream>
#include<fstream>
#include <vector>
using namespace std;
#include "Elliptic_Mesh.h"

int main()
{
	read_boundary_points();
}

void read_boundary_points()
{
	//read_BC("naca0012/BC1.x");
	//read_BC("naca0012/BC2.x");
	read_BC("naca0012/BC3.x");
}

void read_BC(string fileName)
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
	int kkk = 1;
}