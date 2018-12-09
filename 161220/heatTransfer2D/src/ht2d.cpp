/*
 * ht2d.cpp
 *
 *  Created on: Dec 19, 2016
 *      Author: Syed Ahmad Raza
 */

#include <iostream>
#include <fstream>
#include <sstream>  // for int to string conversion
#include <cmath>

using namespace std;

// value of thermal diffusivity assumed for pure silver (99.9%)
const double alpha = 1.6563e-4;

const double timeStep = 0.001;

// boundary conditions
const double T_xi = 25.0;
const double T_xf = 50.0;
const double T_yi = 25.0;
const double T_yf = 50.0;

// initial conditions
const double T_i = 25.0;

// length and width of the specimen
const double length = 1.0;
const double width = 1.0;

// number of intervals on the x and y axes and their respective values
const int nx = 300;
const int ny = 300;
const double x = 0.0; const double y = 0.0;
const double deltaX = length / (nx * 1.0);
const double deltaY = width / (ny * 1.0);

double T[nx + 1][ny + 1] = { 0.0 };

void gridFiler()
{
    ofstream fileX; ofstream fileY;
    fileX.open("./data/coordinatesX.dat");
    fileY.open("./data/coordinatesY.dat");
    for (int i = 0; i <= nx; ++i)
    {
        fileX << x + (i * deltaX) << "\t";
    }
    for (int j = 0; j <= nx; ++j)
    {
        fileY << y + (j * deltaY) << "\t";
    }
}

void solFiler(string fileName)
{
    ofstream file;
        file.open(fileName);
        for (int i = 0; i <= nx; ++i)
        {
            for (int j = 0; j <= ny; ++j)
            {
                file << T[i][j] << "\t";
            }
            file << endl;
        }
}

void solver(double deltaT, int time, double T[][ny + 1])
{
    const double totalTimeSteps = time / deltaT;

    // Applying the boundary conditions
    // Two vertical boundaries
    for (int i = 0; i <= ny; ++i)
    {
        T[i][0] = T_xi;
        T[i][ny] = T_xf;
    }
    // Two horizontal boundaries
    for (int j = 0; j <= nx; ++j)
    {
        T[0][j] = T_yi;
        T[ny][j] = T_yf;
    }
    // Applying the initial condition
    for (int i = 1; i < nx; ++i)
    {
        for (int j = 1; j < ny; ++j)
            T[i][j] = T_i;
    }
    // Determining the numerical solution
    for (int t = 1; t <= totalTimeSteps; ++t)
    {
        for (int i = 1; i < nx; ++i)
        {
            for (int j = 1; j < ny; ++j)
            {
                T[i][j] = T[i][j] + alpha * deltaT *
                        (
                         ((T[i + 1][j] - 2 * T[i][j] + T[i - 1][j])
                                 / (deltaX * deltaX))
                         + ((T[i][j + 1] - 2 * T[i][j] + T[i][j - 1])
                                 / (deltaY * deltaY))
                        );
            }
        }
        // Filing the result at various values of time
        if ((t == 100) || ((t % 1000) == 0))
        {
            ostringstream convertT; convertT << (t * deltaT);
            solFiler(("./data/ht2dCase02T" + convertT.str() + ".dat").c_str());
        }

    }
}

int main()
{
    gridFiler();
    solver(timeStep, 3000, T);

}
