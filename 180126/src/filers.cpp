/*
 * filers.cpp
 *
 *  Created on: 2017-09-27
 *      Author: Syed Ahmad Raza
 */

#include <iostream>     // functions for input and output to console
#include <iomanip>      // functions for setting width
#include <fstream>      // functions for file output
#include <sstream>      // for file name conversion
#include <string>       // for fileName
#include <limits>       // for "numeric_limits<double>::digits10"
                        // it is used to set the maximum precision possible

#include "constants.h"
#include "filers.h"

using namespace std;

ofstream progressFile;
const int mx = nx;
const int my = ny;

void filer1(double Z[], int nk, string fileName)
/*
Function to file a one-dimensional array.
*/
{
    ofstream fileZ;
    fileZ.precision(numeric_limits<double>::digits10 + 2);
    fileZ.open(filePath + fileName + fileExt);
    for (int k = 0; k < nk; ++k)
    {
        fileZ << Z[k] << "\n";
    }
}

void filer2(double Q[][my], int mx, int my, string fileName)
/*
Function to file a two-dimensional array.
*/
{
    ofstream fileQ;
    fileQ.precision(numeric_limits<double>::digits10 + 2);
    fileQ.open(filePath + fileName + fileExt);
    for (int j = 0; j < my; ++j)
    {
        for (int i = 0; i < mx; ++i)
        {
            fileQ << Q[i][j] << "\t";
        }
        fileQ << "\n";
    }
}

void filerAllSol(int timeStep)
{
    // default argument for timeStep = -1; defined in filers.h
    if (timeStep == -1)
    {
        filer2(U, nx, ny,
               velScheme + "_ny=" + to_string(ny) + "_numerical_U-Final");
        filer2(V, nx, ny,
               velScheme + "_ny=" + to_string(ny) + "_numerical_V-Final");
        filer2(P, nx, ny,
               velScheme + "_ny=" + to_string(ny) + "_numerical_P-Final");
    }
    else
    {
        string timeStr = to_string(timeStep);
        // // If to_string does not work (possible in some compilers), use:
        // ostringstream timeStrTemp;
        // timeStrTemp << timeStep;
        // string timeStr = timeStrTemp.str();
        filer2(U, nx, ny,
               velScheme + "_ny=" + to_string(ny) + "_numerical_U-" + timeStr);
        filer2(V, nx, ny,
               velScheme + "_ny=" + to_string(ny) + "_numerical_V-" + timeStr);
        filer2(P, nx, ny,
               velScheme + "_ny=" + to_string(ny) + "_numerical_P-" + timeStr);
    }
}

void progressFileCreator()
{
    progressFile.open(filePath + velScheme + "_ny=" + to_string(ny)
           + "_progressLog" + fileExt);
}

void progressFiler()
{
    // progressFile.open(filePath + velScheme + "_ny=" + to_string(ny)
           // + "_progressLog" + fileExt, ios_base::app);    // open for editing
    // Write column headers
    if (t == 1)
    {
        progressFile << setw(12) << "Timestep"
                     << setw(12) << "P iters"
                     << setw(17) << "Max. mChange"
                     << setw(17) << "Max. pChange"
                     << setw(17) << "Max. uChange"
                     << endl << endl;
    }
    // Write timestep data
    progressFile << setw(12) << t
                 << setw(12) << pIter
                 << setw(17) << mChangeMax
                 << setw(17) << pChangeMax
                 << setw(17) << uChangeMax
                 << endl;
}

void infoFiler()
{
    ofstream f;
    f.open(filePath + velScheme + "_ny=" + to_string(ny)
           + "_characteristics" + fileExt);

    // Write heading and column headers
    f << "Finite Volume Method case for 2D flow inside a pipe\n";
    f << endl;
    f << setw(15) << "Variable" << "Values" << "\n";
    f << endl;

    f << setw(15) << "velScheme" << velScheme << "\n";
    f << setw(15) << "x_y_ratio" << x_y_ratio << "\n";
    f << setw(15) << "D" << D << " m\n";
    f << setw(15) << "L" << L << " m\n";
    f << setw(15) << "ny" << ny << "\n";
    f << setw(15) << "nx" << nx << "\n";
    f << endl;

    f << setw(15) << "nu" << nu << "\n";
    f << setw(15) << "rho" << rho << "\n";
    f << endl;

    f << setw(15) << "omega" << omega << "\n";
    f << setw(15) << "maxTimesteps" << maxTimesteps << "\n";
    f << setw(15) << "maxPressIters" << maxPressIters << "\n";
    f << setw(15) << "uResidual" << uResidual << "\n";
    f << setw(15) << "vResidual" << vResidual << "\n";
    f << setw(15) << "pResidual" << pResidual << "\n";
    f << setw(15) << "dt" << dt << " s\n";
    f << endl;

    f << setw(15) << "uin" << uin << " m/s\n";
    f << setw(15) << "vin" << vin << " m/s\n";
    f << setw(15) << "pbegin" << pbegin << " Pa\n";
    // f << setw(15) << "pin" << pin << " Pa (unused for Neumann)\n";
    // f << setw(15) << "pout" << pout << " Pa (unused for Neumann)\n";
    f << endl;

    // Write data from analytical function
    f << "Data from analytical function\n\n";
    f << setw(15) << "uAvg" << uAvg << " m/s\n";
    streamsize sf = f.precision();   // store current precision value
    f.precision(numeric_limits<double>::digits10 + 2);
    f << setw(15) << "pavg_inlet" << pavg_inlet << " Pa\n";
    f << setw(15) << "pavg_outlet" << pavg_outlet << " Pa\n";
    f.precision(sf);                 // revert to previous precision value
    f << setw(15) << "deltaP" << deltaP << " Pa\n";
    f << endl;

    // Write calculated variables
    f << "Data from calculated variables\n\n";
    f << setw(15) << "Re" << Re << "\n";
    f << setw(15) << "entryLength" << entryLength << " m\n";
    f << setw(15) << "CFL" << CFL << "\n";
    f << endl;

    // Final timestep and l2Norm
    f << setw(15) << "Final timestep" << t << "\n\n";
        // will not work when code is run without navierFVD();
    f << setprecision(numeric_limits<double>::digits10 + 2)
      << setw(15) << "l2Norm" << l2Norm << "\n";
    f << endl;
}
