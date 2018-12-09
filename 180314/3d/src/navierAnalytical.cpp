/*
 * navierAnalytical.cpp
 *
 *  Created on: 2017-09-27
 *      Author: Syed Ahmad Raza
 */

#include <cmath>        // math functions
#include <fstream>      // functions for file operations
#include <sstream>      // string to number conversion
#include <iostream>      // string to number conversion
#include <limits>
// #include <algorithm>    // functions for ranges of elements (arrays)

#include "constants.h"
#include "gridder.h"
#include "filers.h"
#include "navierFVD.h"
#include "navierAnalytical.h"

using namespace std;

double l2Norm       = 0.0;
double pavg_inlet   = 0.0;
double pavg_outlet  = 0.0;
double deltaP       = 0.0;
double uAvg         = 0.0;

double Ua[ny]     = {0.0};
double Va[ny]     = {0.0};

/*
averagePressure function is to be used to calculate the average pressure at
the inlet for the case when a Neumann boundary condition is used at inlet in
the numerical solution.
*/
void averagePressure()
{
    // Following loop can be omitted if running along navierFVD()
    ifstream pressureData;
    string row;
    pressureData.open(filePath + velScheme + "_ny-" + to_string(ny)
                      + "_numerical_P-Final" + fileExt);
    int j = 0;
    while (getline(pressureData, row))
    {
        stringstream ss(row);
        for (int i = 0; i < nx-1; ++i)
        {
            ss >> P[i][j];
        }
        ++j;
    }
    // For LHS pressure:
    for (int j = 0; j < ny; ++j)
    {
        pavg_inlet += P[0][j];
    }
    pavg_inlet = pavg_inlet / static_cast<double>(ny);
    // For RHS pressure:
    for (int j = 0; j < ny; ++j)
    {
        pavg_outlet += P[nx-1][j];
    }
    pavg_outlet = pavg_outlet / static_cast<double>(ny);
}

void navierAnalytical()
{
    averagePressure();   // for Neumann boundary condition at inlet
    // pavg_inlet = pin;   // for Dirichlet condition at inlet (optional)
    // pavg_outlet = pout;   // for Dirichlet condition at outlet (optional)
    // pavg_inlet = 1226800;
    // pavg_outlet = 1226640;
    deltaP = pavg_outlet - pavg_inlet;
    uAvg = D * D * (-deltaP) / (12 * mu * L);

    /*
    The below formula for analytical solution assumes y values starting from
    central axis of the 2D pipe. Therefore, the solution is coded in two steps.
    */
    int q = static_cast<int>(ceil(static_cast<double>(2*ny-3) / 2)) - 2;
    double y = 0.0;
    // cout.precision(numeric_limits<double>::digits10 + 2);
    cout << "j\t" << "q\t" << "y\t\t" << "y/D" << "\n";
    for (int j = 0; j < static_cast<int>(ceil(static_cast<double>(ny-2) / 2));
        ++j)
    {
        y = -Y[q];
        Ua[j] = uAvg * 1.5 * (1.0 - ((2.0 * y / D) * (2.0 * y / D)));
        cout<< j << "\t" << q << "\t" << y << "\t\t" << y / D << "\n";
        ----q;
    }
    cout << "\nj\t" << "q\t" << "y\t\t" << "y/D" << "\n";
    q = 1;
    for (int j = static_cast<int>(ceil(static_cast<double>(ny-2) / 2));
        j < (ny-2); ++j)
    {
        y = Y[q];
        Ua[j] = uAvg * 1.5 * (1.0 - ((2.0 * y / D) * (2.0 * y / D)));
        cout<< j << "\t" << q << "\t" << y << "\t\t" << y / D << "\n";
        ++++q;
    }
    filer1(Ua, (2*ny-3), velScheme + "_ny-" + to_string(ny) + "_exact_U");
}

void navierComparison()
{
    l2Norm = 0.0;

    // Following loop can be omitted if running along navierFVD()
    ifstream uVelocityData;
    string row;
    uVelocityData.open(filePath + velScheme + "_ny-" + to_string(ny)
                       + "_numerical_U-Final" + fileExt);
    int j = 0;
    while (getline(uVelocityData, row))
    {
        stringstream ss(row);
        for (int i = 0; i < nx-1; ++i)
        {
            ss >> U[i][j];
        }
        ++j;
    }
    for (int j = 0; j < ny-1; ++j)
    {
        l2Norm += ((Ua[j] - U[3*nx/4][j]) * (Ua[j] - U[3*nx/4][j]));
    }
    l2Norm = pow(l2Norm / static_cast<double>(ny), 0.5);
}
