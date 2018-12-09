/*
* laplaceUnequalGridAnalytical.cpp
*
*  Created on: 2017-04-11
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

// Boundary conditions
double Q1 = 1.0;
double Q2 = 2.0;

// Indexers
int i = 0;
int j = 0;
int m = 0;

double Q[nx][ny] = { {0.0} };       // Phi
double X[nx] = {0.0};               // x coordinates
double Y[ny] = {0.0};               // y coordinates
double Q3 = Q2 - Q1;
double lambda = 0.0;
double summation = 0.0;

void gridder()
{
    for (i = 0; i < nx; ++i)
    {
        X[i] = (0.5 * L) * (1 + sin( pi * ((i * 1.0 / (nx-1)) - 0.5)));
    }
    for (j = 0; j < ny; ++j)
    {
        Y[j] = W * (sin( pi * (j * 1.0 / (ny-1)) * 0.5));
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

void analytical()
{
    // Generating and filing the grid;
    gridder();

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

    for (i = 1; i < nx-1; ++i)
    {
        for (j = 1; j < ny-1; ++j)
        {
            summation = 0;
            for (m = 1; m <= 200; ++m)
            {
                lambda = (1.0 * m * pi) / L;
                summation += ( (pow(-1, m+1) + 1) / m )
                             * (
                                sin(lambda * X[i]) * sinh(lambda * Y[j])
                                / sinh(lambda * W)
                               );
            }
            Q[i][j] = (2.0 * Q3 / pi) * summation + Q1;
        }
    }
    solFiler("../data/laplaceAnalytical.dat");
    errorFiler("../data/laplaceAnalyticalError.dat");
}

int main()
{
    analytical();
    cout << X[nx/2];    // x = 0.507933
    return 0 ;
}
