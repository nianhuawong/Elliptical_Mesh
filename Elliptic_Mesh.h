#pragma once
int NI, NJ;
vector < vector<double> > globalCoordX;
vector < vector<double> > globalCoordY;

void read_boundary_points();
void read_BC(string fileName, vector < vector<double> > &globalCoordX, vector < vector<double> >& globalCoordY, string pos);
void initialize();
