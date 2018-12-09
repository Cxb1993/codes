/*
 * main.cpp
 *
 *  Created on: 2017-09-27
 *      Author: Syed Ahmad Raza
 */

#include <cmath>        // math functions
#include <iostream>     // functions for input and output to console
#include <fstream>      // functions for file output
// #include <algorithm>    // functions for ranges of elements (arrays)

#include "constants.h"
#include "gridder.h"
#include "filers.h"
#include "navierFVD.h"
#include "navierAnalytical.h"

using namespace std;

// g++ -std=c++14 -Wall gridder.cpp filers.cpp navierFVD.cpp navierAnalytical.cpp main.cpp -o main & ./main

void infoWriter()
{
    ofstream fileInfo;
    fileInfo.open("../data/characteristics.txt");
    fileInfo << "Finite Volume Method case for 2D flow inside a pipe\n";
    fileInfo << '\n';
    fileInfo << "Variable\t\tValues\n";
    fileInfo << "x_y_ratio\t\t" << x_y_ratio << '\n';
    fileInfo << "D\t\t\t\t" << D << " m\n";
    fileInfo << "L\t\t\t\t" << L << " m\n";
    fileInfo << "ny\t\t\t\t" << ny << '\n';
    fileInfo << "nx\t\t\t\t" << nx << '\n';
    fileInfo << '\n';
    fileInfo << "uin\t\t\t\t" << uin << " m/s\n";
    fileInfo << "vin\t\t\t\t" << vin << " m/s\n";
    fileInfo << "pconst\t\t\t" << pconst << " Pa\n";
    fileInfo << '\n';
    fileInfo << "nu\t\t\t\t" << nu << '\n';
    fileInfo << "rho\t\t\t\t" << rho << '\n';
    fileInfo << '\n';
    fileInfo << "omega\t\t\t" << omega << '\n';
    fileInfo << "maxTimesteps\t" << maxTimesteps << '\n';
    fileInfo << "maxPressIters\t" << maxPressIters << '\n';
    fileInfo << "uResidual\t\t" << uResidual << '\n';
    fileInfo << "vResidual\t\t" << vResidual << '\n';
    fileInfo << "pResidual\t\t" << pResidual << '\n';
    fileInfo << "dt\t\t\t\t" << dt << " s\n";
    fileInfo << '\n';
    fileInfo << "Re\t\t\t\t" << Re << '\n';
    fileInfo << "Entry length\t" << 0.05 * Re * D << " m\n";
    fileInfo << "CFL\t\t\t\t" << (uin * dt / dx) + (vin * dt / dy) << '\n';
}

int main()
{
    infoWriter();
    // navierFVD();
    navierAnalytical();
    navierConvergence();

    // Write the final time step and L2 norm to characteristics file
    ofstream fileInfo;
    fileInfo.open("../data/characteristics.txt", ios_base::app); // open for editing
    fileInfo << '\n' << "Final time step is " << t << '\n'; // will not work when code is run without navierFVD();
    fileInfo << '\n' << "L2 norm is " << l2Norm << '\n';

    return 0;
}
