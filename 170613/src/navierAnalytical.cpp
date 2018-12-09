/*
 * navierAnalytical.cpp
 *
 *  Created on: 2017-06-13
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

double Ua[nx][ny] = { {0.0} };
double Va[nx][ny] = { {0.0} };;
double Pavg[ny] = {0.0};

void averagePressure()
{
    ifstream pressureData;
    string row;
    string num;
    double temp2;
    pressureData.open("../data/navierFVD_P_Final.dat");

    int j = 0;
    while (getline(pressureData, row))
    {
        stringstream ss(row);
        ss >> Pavg[j];
        cout << Pavg[j] << '\n';
        // row >> num;
        // row >> temp1;
        // temp1 >> Pavg[j];
        // cout << Pavg[j];
        ++j;
    }
}

// void analytical()
// {
//     double deltaP = 368.0;
//     double Uavg = ly*ly*deltaP/(12*mu*lx);
//     for (int i = 0; i < nx; ++i)
//     {
//         for (int j = 0; j < ny; ++j)
//         {
//             Ua[i][j] = (3.0/2.0)*(1.0 - (Y[j]/(ly/2.0))*(Y[j]/(ly/2.0)))*Uavg;
//         }
//     }
//     filer2(Ua, nx, ny, "navierAnalytical_U");
// }

int main()
{
    averagePressure();
    return 0;
}
