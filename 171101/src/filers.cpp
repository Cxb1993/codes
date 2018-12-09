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

#include "constants.h"
#include "filers.h"

using namespace std;

void filer1(double Z[], int nk, string fileName)
/*
Function to output one-dimensional array to a folder named 'data' in the
parent directory.
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

void filer2(double Q[][ny], int nx, int ny, string fileName)
/*
Function to output two-dimensional array to a folder named 'data' in the
parent directory.
*/
{
    ofstream fileQ;
    fileQ.precision(numeric_limits<double>::digits10 + 2);
    fileQ.open(filePath + fileName + fileExt);
    for (int j = 0; j < ny; ++j)
    {
        for (int i = 0; i < nx; ++i)
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
               velScheme + "_ny-" + to_string(ny) + "_numerical_U-Final");
        filer2(V, nx, ny,
               velScheme + "_ny-" + to_string(ny) + "_numerical_V-Final");
        filer2(P, nx, ny,
               velScheme + "_ny-" + to_string(ny) + "_numerical_P-Final");
    }
    else
    {
        string timeStr = to_string(timeStep);
        // // If to_string does not work (possible in some compilers), use:
        // ostringstream timeStrTemp;
        // timeStrTemp << timeStep;
        // string timeStr = timeStrTemp.str();
        filer2(U, nx, ny,
               velScheme + "_ny-" + to_string(ny) + "_numerical_U-" + timeStr);
        filer2(V, nx, ny,
               velScheme + "_ny-" + to_string(ny) + "_numerical_V-" + timeStr);
        filer2(P, nx, ny,
               velScheme + "_ny-" + to_string(ny) + "_numerical_P-" + timeStr);
    }
}

void progressPrinter()
{
    {
        // Print column headers
        if (t == 1)
        {
            printf("\n%9s %7s %16s \t %16s\n", "Time step", "P iters", "deltaP_max", "deltaU_max");
        }
        // Print timestep data
        printf("%9d %7d %16f \t %16e\n", t, pIter, pChange, uChange);
    }
}

void infoWriter()
{
    // FILE WRITER

    // Create and open file
    ofstream f;
    // f.precision(numeric_limits<double>::digits10 + 2);
    f.open(filePath + velScheme + "_ny-" + to_string(ny) + "_characteristics" + fileExt);

    // Write heading and column headers
    f << "Finite Volume Method case for 2D flow inside a pipe\n";
    f << endl;
    f << setw(15) << left << "Variable" << "Values" << "\n";
    f << endl;

    f << setw(15) << left << "velScheme" << velScheme << "\n";
    f << setw(15) << left << "x_y_ratio" << x_y_ratio << "\n";
    f << setw(15) << left << "D" << D << " m\n";
    f << setw(15) << left << "L" << L << " m\n";
    f << setw(15) << left << "ny" << ny << "\n";
    f << setw(15) << left << "nx" << nx << "\n";
    f << endl;

    f << setw(15) << left << "nu" << nu << "\n";
    f << setw(15) << left << "rho" << rho << "\n";
    f << endl;

    f << setw(15) << left << "omega" << omega << "\n";
    f << setw(15) << left << "maxTimesteps" << maxTimesteps << "\n";
    f << setw(15) << left << "maxPressIters" << maxPressIters << "\n";
    f << setw(15) << left << "uResidual" << uResidual << "\n";
    f << setw(15) << left << "vResidual" << vResidual << "\n";
    f << setw(15) << left << "pResidual" << pResidual << "\n";
    f << setw(15) << left << "dt" << dt << " s\n";
    f << endl;

    f << setw(15) << left << "uin" << uin << " m/s\n";
    f << setw(15) << left << "vin" << vin << " m/s\n";
    f << setw(15) << left << "pbegin" << pbegin << " Pa\n";
    f << setw(15) << left << "pin" << pin << " Pa (unused for Neumann)\n";
    f << setw(15) << left << "pout" << pout << " Pa (unused for Neumann)\n";
    f << endl;

    // Write data from analytical function
    f << "Data from analytical function\n\n";
    f << setw(15) << left << "uAvg" << uAvg << " m/s\n";
    streamsize sf = f.precision();   // store current precision value
    f.precision(numeric_limits<double>::digits10 + 2);
    f << setw(15) << left << "pavg_inlet" << pavg_inlet << " Pa\n";
    f << setw(15) << left << "pavg_outlet" << pavg_outlet << " Pa\n";
    f.precision(sf);                 // revert to previous precision value
    f << setw(15) << left << "deltaP" << deltaP << " Pa\n";
    f << endl;

    // Write calculated variables
    f << "Data from calculated variables\n\n";
    f << setw(15) << left << "Re" << Re << "\n";
    f << setw(15) << left << "entryLength" << entryLength << " m\n";
    f << setw(15) << left << "CFL" << CFL << "\n";
    f << endl;

    // f.open(filePath + velScheme + "_characteristics" + fileExt,
    //               ios_base::app);    // open for editing

    f << setw(15) << left << "Final timestep" << t << "\n\n";
        // will not work when code is run without navierFVD();
    f << setprecision(numeric_limits<double>::digits10 + 2)
             << setw(15) << left << "l2Norm" << l2Norm << "\n";
    f << endl;


    // SCREEN PRINTER

    // Print final timestep
    cout << endl << setw(15) << left << "Final timestep" << t << "\n\n";

    // Print analytical data
    cout << "Analytical data\n";
    cout << setw(15) << left << "uAvg" << uAvg << "\n";
    cout << setw(15) << left << "pavg_inlet" << pavg_inlet << "\n";
    cout << setw(15) << left << "pavg_outlet" << pavg_outlet << "\n";
    cout << setw(15) << left << "deltaP" << deltaP << "\n\n";

    // Print case characteristics
    cout << "Case characteristics\n";
    cout << setw(15) << left << "velScheme" << velScheme << "\n";
    cout << setw(15) << left << "Re" << Re << "\n";
    cout << setw(15) << left << "entryLength" << entryLength << "\n";
    cout << setw(15) << left << "CFL" << CFL << "\n\n";

    streamsize sc = cout.precision();   // store current precision value
    cout.precision(numeric_limits<double>::digits10 + 2);
    cout << setw(15) << left << "l2Norm" << l2Norm << "\n\n";
    cout.precision(sc);                 // revert to previous precision value

    if (L <= entryLength)
    {
        cout << "WARNING: Boundary Layer is not fully developed.\n\n";
    }
}
