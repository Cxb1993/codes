/*
 * analytical.cpp
 *
 *  Created on: Mar 23, 2017
 *      Author: Syed Ahmad Raza
 */

#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

const double pi(3.14159);

// Domain and range
const double L = 1.0;
const double W = 1.0;

// Number of nodes on the x and y axes
const int nx = 50;
const int ny = 50;

// Boundary conditions
double Q1 = 1.0;
double Q2 = 2.0;

// Indexers
int i = 0;
int j = 0;
int n = 0;

double Q[nx][ny] = { {0.0} };

void analytical()
{
    // Phi array and variable initialization

    double Q3 = Q2 - Q1;
    double lambda = 0.0;
    double summation = 0.0;
    double x = 0.0;     // using equal grid intervals
    double y = 0.0;

    // Two vertical boundary conditions
    for (j = 0; j < ny; ++j)
    {
        Q[0][j] = Q1;
        Q[nx-1][j] = Q1;
    }

    // Two horizontal boundary conditions
    for (i = 0; i < nx; ++i)
    {
        Q[i][0] = Q1;
        Q[i][ny-1] = Q2;
    }

    for (i = 1; i < nx - 1; ++i)
    {
        x += L / (nx - 1);
        cout << "x = " << x << '\n';
        y = 0.0;
        for (j = 1; j < ny - 1; ++j)
        {
            y += W / (ny - 1);
            cout << "y = " << y << '\n';
            summation = 0;
            for (n = 1; n <= 200; ++n)
            {
                lambda = (1.0 * n * pi) / L;
                summation += ( (pow(-1, n+1) + 1) / n )
                             * (
                                sin(lambda * x) * sinh(lambda * y)
                                / sinh(lambda * W)
                               );
            }
            Q[i][j] = (2.0 * Q3 / pi) * summation + Q1;
        }
    }
    // Writing data to an output file
    ofstream file;
    file.open("../data/laplaceAnalytical.dat");
    for (j = 0; j < ny; ++j)
    {
        for (i = 0; i < nx; ++i)
        {
            file << Q[i][j] << "\t";
        }
        file << "\n";
    }
    // Writing data for error calculation
    ofstream fileE;
    fileE.open("../data/laplaceAnalyticalError.dat");
    for (j = 0; j < ny; ++j)
    {
        fileE << Q[25][j] << '\t';
    }
    fileE << '\n';
    for (j = 0; j < ny; ++j)
    {
        fileE << j * (W / (ny - 1)) << '\t';
    }
}

int main()
{
    analytical();
    return 0 ;
}
