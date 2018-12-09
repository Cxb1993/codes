/*
 * navierAnalytical.cpp
 *
 *  Created on: 2017-09-27
 *      Author: Syed Ahmad Raza
 */

#include <cmath>        // math functions
#include <fstream>      // functions for file input
#include <iostream>     // functions for input and output to console
#include <sstream>      // string to number conversion
// #include <algorithm>    // functions for ranges of elements (arrays)

#include "constants.h"
#include "gridder.h"
#include "filers.h"
#include "navierAnalytical.h"

using namespace std;

double Ua[ny] = {0.0};
double Va[ny] = {0.0};
double Pavg[ny] = {0.0};
double Ulast[ny] = {0.0};
double pavg_inlet = 0.0;
double l2Norm = 0.0;

void averagePressure()
{
    // Use following if running the code as standalone
    ifstream pressureData;
    string row;
    pressureData.open("../data/navierFVD_P_Final.dat");
    int j = 0;
    while (getline(pressureData, row))
    {
        stringstream ss(row);
        ss >> Pavg[j];
        pavg_inlet += Pavg[j];
        ++j;
    }

    // // Use following if running along navierFVD()
    // for (int j = 0; j < ny; ++j)
    // {
    //     pavg_inlet += P[nx-1][j];
    // }

    pavg_inlet = pavg_inlet / static_cast<double>(ny);
}

void navierAnalytical()
{
    averagePressure();
    double deltaP = pconst - pavg_inlet;
    double Uavg = D * D * (- deltaP) / (12 * mu * L);

    // // Use following if running the code as standalone
    // ifstream coordinateYData;
    // string row;
    // coordinateYData.open("../data/coordinateY.dat");
    // int j = 0;
    // while(getline(coordinateYData, row))
    // {
    //     stringstream ss(row);
    //     ss >> Y[j];
    //     ++j;
    // }

    for (int j = 0; j < (ny / 2); ++j)
    {
        Ua[j] = Uavg * 1.5 * (1.0 - ((2.0 * Y[j] / D) * (2.0 * Y[j] / D)));
        // Note that the above formula assumes y values starting from central axis of the 2D pipe
    }
    filer1(Ua, ny, "navierAnalytical_U");
    cout << '\n' << "deltaP is "<< deltaP << '\n';
    cout << '\n' << "Uavg is "<< Uavg << '\n';
}

void navierConvergence()
{
    l2Norm = 0.0;

    // // Use following if running the code as standalone
    // ifstream uVelocityData;
    // string row;
    // uVelocityData.open("../data/navierFVD_U_Final.dat");
    // int j = 0;
    // while (getline(uVelocityData, row))
    // {
    //     stringstream ss(row);
    //     double temp = 0.0;
    //     for (int i = 0; i < (nx-1); ++i)
    //     {
    //         ss >> temp;
    //     }
    //     ss >> Ulast[j];
    //     ++j;
    // }
    // for (int j = 0; j < (ny / 2); ++j)
    // {
    //     l2Norm += (Ua[j] - Ulast[j+(ny/2)]) * (Ua[j] - Ulast[j+(ny/2)]);    // Currently, this will only work for equal grids because of the formulation of naveirAnalytical()
    // }

    // Use following if running along navierFVD()
    for (int j = 0; j < (ny / 2); ++j)
    {
        l2Norm += (Ua[j] - U[nx-1][j+(ny/2)]) * (Ua[j] - U[nx-1][j+(ny/2)]);    // Currently, this will only work for equal grids because of the formulation of naveirAnalytical()
    }

    l2Norm = pow(l2Norm / (ny / 2), 0.5);
}
