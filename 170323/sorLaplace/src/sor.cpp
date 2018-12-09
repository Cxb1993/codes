/*
 * sor.cpp
 *
 *  Created on: Mar 27, 2017
 *      Author: Syed Ahmad Raza
 */

#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

using namespace std;

const double pi(3.14159);

// Domain and range
const double L = 1.0;
const double W = 1.0;

// Number of nodes on the x and y axes
const int nx = 50;
const int ny = 50;

// Intervals on x and y axes
const double dx = L/(1.0*nx - 1);
const double dy = W/(1.0*ny - 1);

// Boundary conditions
double Q1 = 1.0;
double Q2 = 2.0;

// Indexers
int i = 0;
int j = 0;
int m = 0;

// Relaxation parameter
double w = 1.8;

double Qmean = 1.0;
double Qmeano = 0.0;

double E = 1.0;

double Q[nx][ny] = { {0.0} };

string fileString = "";

void fileWriter(string fileName)
{
    ofstream file;
    file.open("../data/laplaceSOR" + fileName + ".dat");
    for (j = 0; j < ny; ++j)
    {
        for (i = 0; i < nx; ++i)
            file << Q[i][j] << '\t';
        file << '\n';
    }
}

void sor(double error)
{
    m = 0;
    // Applying boundary condtions
    for (j = 0; j < ny; ++j)
    {
        Q[0][j] = Q1;
        Q[nx-1][j] = Q1;
    }
    for (i = 0; i < nx; ++i)
    {
        Q[i][0] = Q1;
        Q[i][ny-1] = Q2;
    }

    // Applying the SOR iterations
    while (E >= error && m <= 10000)
    {
        ++m;
        for (i = 1; i < nx - 1; ++i)
        {
            for (j = 1; j < ny - 1; ++j)
            {
                Q[i][j] = (1.0 - w)*Q[i][j] + w*
                         (
                          (Q[i+1][j] + Q[i-1][j])/dx/dx +
                          (Q[i][j+1] + Q[i][j-1])/dy/dy
                         )*
                         (
                          dx*dx*dy*dy/2.0/(dx*dx + dy*dy)
                         );
            }
        }
        Qmeano = Qmean;
        Qmean = 0.0;
        for (i = 0; i < nx; ++i)
            for (j = 0; j < ny; ++j)
                Qmean += Q[i][j];
        Qmean = Qmean/(1.0*nx*ny);
        E = fabs(Qmean - Qmeano);

//        // Writing iteration data to an output file
//        if (m % 100 == 0)
//        {
//            ostringstream mStr;
//            mStr << m;
//            fileString = mStr.str();
//            fileWriter(fileString);
//        }
    }
    ostringstream eStr;
    eStr << error;
    fileWriter("Final" + eStr.str());

    // Writing data for error calculation
    ofstream fileE;
    fileE.open("../data/laplaceSORError" + eStr.str() + ".dat");
    for (j = 0; j < ny; ++j)
    {
        fileE << Q[25][j] << '\t';
    }
    fileE << '\n';
    for (j = 0; j < ny; ++j)
    {
        fileE << j*dy << '\t';
    }
    cout << error << '\t' << m << '\n';
}

int main()
{
    sor(1e-2);
    sor(1e-3);
    sor(1e-4);
    sor(1e-5);
    return 0 ;
}
