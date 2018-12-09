/*
* laplaceUnequalGridNumerical.cpp
*
*  Created on: 2017-04-09
*      Author: Syed Ahmad Raza
*/

#include <cmath>        // math functions
#include <iostream>     // functions for input and output to console
#include <fstream>      // functions for file input and output
#include <sstream>      // functions for string conversion
#include <algorithm>    // functions for ranges of elements (arrays)

using namespace std;

const double pi(3.14159);

// Domain and range
const double L = 1.0;
const double W = 1.0;

// Number of nodes on the x and y axes
const int nx = 100;
const int ny = 100;

// Boundary Conditions
double Q1 = 1.0;
double Q2 = 2.0;

// Relaxation parameter
double w = 1.8;

// Indexers and other variables
int i = 0;
int j = 0;
int m = 0;
double diffMax = 0.1;
double tolerance = 1e-9;

// Arrays
double Q[nx][ny] = { {0.0} };       // Phi
double Qo[nx][ny] = { {0.0} };      // Phi old
double X[nx] = {0.0};               // x coordinates
double Xd[nx-1] = {0.0};            // x interval sizes
double Y[ny] = {0.0};               // y coordinates
double Yd[ny-1] = {0.0};            // y interval sizes

// Intervals on x and y axes
const double dx = L/1.0*(nx-1);
const double dy = W/1.0*(ny-1);

void gridderFiler()
{
    ofstream fileX;
    fileX.open("../data/coordinatesX.dat");
    for (i = 0; i < nx; ++i)
    {
        X[i] = (0.5 * L) * (1 + sin( pi * ((i * 1.0 / (nx-1)) - 0.5)));
        fileX << X[i] << '\n';
    }
    ofstream fileXd;
    fileXd.open("../data/coordinatesXd.dat");
    for (i = 0; i < nx-1; ++i)
    {
        Xd[i] = X[i+1] - X[i];
        fileXd << Xd[i] << '\n';
    }

    ofstream fileY;
    fileY.open("../data/coordinatesY.dat");
    for (j = 0; j < ny; ++j)
    {
        Y[j] = W * (sin( pi * (j * 1.0 / (ny-1)) * 0.5));
        fileY << Y[j] << '\n';
    }
    ofstream fileYd;
    fileYd.open("../data/coordinatesYd.dat");
    for (j = 0; j < ny-1; ++j)
    {
        Yd[j] = Y[j+1] - Y[j];
        fileYd << Yd[j] << '\n';
    }
}

void solFiler(string fileName)
{
    ofstream file;
    file.open(fileName);
    for (j = 0; j < ny; ++j)
    {
        for (i = 0; i < nx; ++i)
        {   file << Q[i][j] << '\t';  }
        file << '\n';
    }
}

void errorFiler(string fileName)
{
    // Writing data for error calculation
    ofstream fileE;
    fileE.open(fileName);
    for (j = 0; j < ny; ++j)
    {
        fileE << Q[nx/2][j] << '\t';
    }
}

void solver()
{
    // Generating and filing the grid;
    gridderFiler();

    // Applying the boundary conditions
    // Vertical boundaries
    for (j = 0; j < ny; ++j)
    {
        Q[0][j] = Q1;
        Q[nx-1][j] = Q1;
    }
    // Horizontal boundaries
    for (i = 0; i < nx; ++i)
    {
        Q[i][0] = Q1;
        Q[i][ny-1] = Q2;
    }

    // Applying iterations using SOR method
    while (m <= 10 || ((diffMax >= tolerance) && m <= 10000))
    {
        cout << ++m << endl;
        for (i = 1; i < nx-1; ++i)
        {
            for (j = 1; j < ny-1; ++j)
            {
                Q[i][j]
                = (1.0 - w)*Q[i][j] +
                w*         // term 2
                (
                 (         // term 2a
                  2*(Q[i+1][j] + Q[i-1][j])/
                  (Xd[i]*Xd[i] + Xd[i-1]*Xd[i-1])
                 ) -
                 (         // term 2b
                  2*(Q[i+1][j] - Q[i-1][j])*(Xd[i] - Xd[i-1])/
                  ( (Xd[i] + Xd[i-1])*(Xd[i]*Xd[i] + Xd[i-1]*Xd[i-1]) )
                 ) +
                 (         // term 2c
                  2*(Q[i][j+1] + Q[i][j-1])/
                  (Yd[j]*Yd[j] + Yd[j-1]*Yd[j-1])
                 ) -
                 (         // term 2d
                  2*(Q[i][j+1] - Q[i][j-1])*(Yd[j] - Yd[j-1])/
                  ( (Yd[j] + Yd[j-1])*(Yd[j]*Yd[j] + Yd[j-1]*Yd[j-1]) )
                 )
                )*
                1.0/4.0*       // term 3
                (
                 (         // term 3 numerator
                  (Xd[i]*Xd[i] + Xd[i-1]*Xd[i-1])*
                  (Yd[j]*Yd[j] + Yd[j-1]*Yd[j-1])
                 )/
                 (         // term 3 denominator
                  Xd[i]*Xd[i] + Xd[i-1]*Xd[i-1] +
                  Yd[j]*Yd[j] + Yd[j-1]*Yd[j-1]
                 )
                );
            }
        }
        diffMax = 0.0;
        for (i = 1; i < nx-1; ++i)
        {
            for (j = 1; j < ny-1; ++j)
            {
                diffMax = max(abs(Q[i][j] - Qo[i][j]), diffMax);
                Qo[i][j] = Q[i][j];
            }
        }

        // Writing iteration data to an output file
        if (m % 100 == 0)
        {
            ostringstream mStr;
            mStr << m;
            solFiler("../data/laplaceUnequalGrid" + mStr.str() + ".dat");
        }
    }
    solFiler("../data/laplaceUnequalGridFinal.dat");
    errorFiler("../data/laplaceNumericalError.dat");
}

int main()
{
    solver();
    return 0;
}
