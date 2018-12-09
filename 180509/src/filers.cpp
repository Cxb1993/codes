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
using namespace Numeric_lib;

ofstream progressFile;


void filer2DArray(const Matrix<double,2>& Q,
                  const Matrix<double,1>& Xp,
                  const Matrix<double,1>& Yp,
                  const int xIndexLimit, const int yIndexLimit,
                  const int xOffset, const int yOffset,
                  const string colHeader, const string fileName)
{
    /*
    Function to file a 2D array.
    */
    ofstream fileQ;
    fileQ.precision(numeric_limits<double>::digits10 + 2);
    fileQ.open(filePath + fileName + ".csv");

    fileQ << "x, y, z, " << colHeader << "\n";

    for (int j = 0; j <= yIndexLimit; ++j)
    {
        for (int i = 0; i <= xIndexLimit; ++i)
        {
            fileQ << Xp(i) << ", " << Yp(j) << ", 0.0, "
                  << Q(i+xOffset,j+yOffset) << "\n";
        }
    }
}

void filer2DSolution()
{
    /*
    Function to file the complete solution of 2D problem.
    */
    ofstream fileQ;
    fileQ.precision(numeric_limits<double>::digits10 + 2);
    fileQ.open(filePath + fileUniqueName + "_" + to_string(nxy) + "_sol"
               + ".csv");

    fileQ << "x, y, z, pressure, Vx, Vy"
          << "\n";

    int offSet = ghostCells / 2;
    for (int j = 0; j <= ny-5; ++j)
    {
        for (int i = 0; i <= nx-5; ++i)
        {
            fileQ << Xs(i) << ", " << Ys(j) << ", 0.0, "
                  << P(i+offSet,j+offSet) << ", "
                  << U(i+offSet,j+offSet) << ", "
                  << V(i+offSet,j+offSet) << "\n";
        }
    }
}

void filer1(const Matrix<double,1>& R, const int nk, const string fileName)
/*
Function to file a one-dimensional array.
*/
{
    ofstream fileR;
    fileR.precision(numeric_limits<double>::digits10 + 2);
    fileR.open(filePath + fileName + fileExt);

    for (int k = 0; k < nk; ++k)
    {
        fileR << R(k) << "\n";
    }
}

void filer2(const Matrix<double,2>& T,
            const int mx, const int my, const string fileName)
/*
Function to file a two-dimensional array.
*/
{
    ofstream fileT;
    fileT.precision(numeric_limits<double>::digits10 + 2);
    fileT.open(filePath + fileName + fileExt);

    for (int j = 0; j < my; ++j)
    {
        for (int i = 0; i < mx; ++i)
        {
            fileT << T(i,j) << "\t";
        }
        fileT << "\n";
    }
}

void filerAllSol(const int timeStep)
{
    // default argument for timeStep = -1; defined in filers.h
    if (timeStep == -1)
    {
        // filer2(U, nx, ny,
        //         fileUniqueName + "_" + to_string(nxy) + "_U-final");
        // filer2(V, nx, ny,
        //         fileUniqueName + "_" + to_string(nxy) + "_V-final");
        // filer2(P, nx, ny,
        //         fileUniqueName + "_" + to_string(nxy) + "_P-final");
        // filer2(MC, nx, ny,
        //         fileUniqueName + "_" + to_string(nxy) + "_MC-final");
        // filer2(PC, nx, ny,
        //         fileUniqueName + "_" + to_string(nxy) + "_PC-final");

        filer2DArray(U, Xa, Ys, nx-4, ny-5, ghostCells/2-1, ghostCells/2,
                     "Vx",
                     fileUniqueName + "_" + to_string(nxy) + "_U-last");
        filer2DArray(V, Xs, Ya, nx-5, ny-4, ghostCells/2, ghostCells/2-1,
                     "Vy",
                     fileUniqueName + "_" + to_string(nxy) + "_V-last");
        filer2DArray(P, Xs, Ys, nx-5, ny-5, ghostCells/2, ghostCells/2,
                     "P",
                     fileUniqueName + "_" + to_string(nxy) + "_P-last");
        filer2DSolution();
    }
    else
    {
        string timeStr = to_string(timeStep);
        /*
        If to_string does not work (possible in some compilers), use:
        ostringstream timeStrTemp;
        timeStrTemp << timeStep;
        string timeStr = timeStrTemp.str();
        */
        // filer2(U, nx, ny,
        //         fileUniqueName + "_" + to_string(nxy) + "_U-" + timeStr);
        // filer2(V, nx, ny,
        //         fileUniqueName + "_" + to_string(nxy) + "_V-" + timeStr);
        // filer2(P, nx, ny,
        //         fileUniqueName + "_" + to_string(nxy) + "_P-" + timeStr);
        // filer2(MC, nx, ny,
        //         fileUniqueName + "_" + to_string(nxy) + velScheme
        //         + "_MC-" + timeStr);
        // filer2(PC, nx, ny,
        //         fileUniqueName + "_" + to_string(nxy) + velScheme
        //         + "_PC-" + timeStr);

        filer2DArray(U, Xa, Ys, nx-4, ny-5, ghostCells/2-1, ghostCells/2,
                     "Vx",
                     fileUniqueName + "_" + to_string(nxy) + "_U-" + timeStr);
        filer2DArray(V, Xs, Ya, nx-5, ny-4, ghostCells/2, ghostCells/2-1,
                     "Vy",
                     fileUniqueName + "_" + to_string(nxy) + "_V-" + timeStr);
        filer2DArray(P, Xs, Ys, nx-5, ny-5, ghostCells/2, ghostCells/2,
                     "P",
                     fileUniqueName + "_" + to_string(nxy) + "_P-" + timeStr);
        filer2DSolution();
    }
}

void filerCoordinates()
{
    /*
    Filing the coordinate files
    */

    filer1(Xa, (nx-3), filePath + fileUniqueName + "_" + to_string(nxy)
            + "_coordinateX");
    filer1(Ya, (ny-3), filePath + fileUniqueName + "_" + to_string(nxy)
            + "_coordinateY");
    filer1(Xs, (nx-4), filePath + fileUniqueName + "_" + to_string(nxy)
            + "_coordinateXs");
    filer1(Ys, (ny-4), filePath + fileUniqueName + "_" + to_string(nxy)
            + "_coordinateYs");

    filer1(Dx, (nx), filePath + fileUniqueName + "_" + to_string(nxy)
            + "_coordinateDx");
    filer1(Dxs, (nx), filePath + fileUniqueName + "_" + to_string(nxy)
            + "_coordinateDxs");

    filer1(Dy, (ny), filePath + fileUniqueName + "_" + to_string(nxy)
            + "_coordinateDy");
    filer1(Dys, (ny), filePath + fileUniqueName + "_" + to_string(nxy)
            + "_coordinateDys");
}

void filerCreateProgress()
{
    progressFile.open(filePath + fileUniqueName + "_" + to_string(nxy)
                      + "_logProgress.txt");
}

void filerProgress()
{
    // Write column headers
    if (t == 1)
    {
        progressFile
            << setw(12) << "Timestep"
            << setw(12) << "P iters"
            << setw(17) << "Max. mChange"
            << setw(17) << "Max. pChange"
            << setw(17) << "Max. uChange"
            << setw(17) << "Max. vChange"
            << endl
            << setw(12) << ""
            << setw(12) << ""
            << setw(11) << "i" << setw(2) << ", "
            << setw(4)  << "j"
            << setw(11) << "i" << setw(2) << ", "
            << setw(4)  << "j"
            << setw(11) << "i" << setw(2) << ", "
            << setw(4)  << "j"
            << setw(11) << "i" << setw(2) << ", "
            << setw(4)  << "j"
            << endl << endl;
    }
    // Write timestep data
    progressFile
        << setw(12) << t
        << setw(12) << pIter
        << setw(17) << mChangeMax
        << setw(17) << pChangeMax
        << setw(17) << uChangeMax
        << setw(17) << vChangeMax
        << endl
        << setw(12) << ""
        << setw(12) << ""
        << setw(11) << mChangeMax_i << setw(2) << ", "
        << setw(4)  << mChangeMax_j
        << setw(11) << pChangeMax_i << setw(2) << ", "
        << setw(4)  << pChangeMax_j
        << setw(11) << uChangeMax_i << setw(2) << ", "
        << setw(4)  << uChangeMax_j
        << setw(11) << vChangeMax_i << setw(2) << ", "
        << setw(4)  << vChangeMax_j
        << endl;
}

void filerInfoStart()
{
    ofstream f;
    f.open(filePath + fileUniqueName + "_" + to_string(nxy)
           + "_logCharacteristics.txt");

    // Write heading and column headers for case characteristics
    f   << left << "Finite Volume Method for 2D flow"
        << endl << endl
        << setw(20) << "Variable" << "Values" << endl
        << endl
        << setw(20) << "velScheme" << velScheme << endl
        << setw(20) << "ny" << ny << endl
        << setw(20) << "nx" << nx << endl
        << setw(20) << "ghostCells" << ghostCells << endl
        // << setw(20) << "CFL" << CFL << endl
        << endl

        << setw(20) << "Re" << Re << endl
        << setw(20) << "nu" << nu << endl
        << setw(20) << "rho" << rho << endl
        << endl

        << setw(20) << "omega" << omega << endl
        << setw(20) << "maxTimesteps" << maxTimesteps << endl
        << setw(20) << "maxPressIters" << maxPressIters << endl
        << setw(20) << "uResidual" << uResidual << endl
        << setw(20) << "vResidual" << vResidual << endl
        << setw(20) << "pResidual" << pResidual << endl
        << setw(20) << "dt" << dt << " s" << endl
        << endl

        << setw(20) << "x_y_ratio" << x_y_ratio << endl
        << setw(20) << "W" << W << " m" << endl
        << setw(20) << "L" << L << " m" << endl
        << setw(20) << "uIn" << uIn << " m/s" << endl
        << setw(20) << "vIn" << vIn << " m/s" << endl
        << setw(20) << "uInitial" << uInitial << " m/s" << endl
        << setw(20) << "vInitial" << vInitial << " m/s" << endl
        << setw(20) << "pInitial" << pInitial << " Pa" << endl
        // << setw(20) << "pIn" << pIn << " Pa (unused for Neumann)" << endl
        // << setw(20) << "pout" << pout << " Pa (unused for Neumann)" << endl
        << endl;

    // Write comments, if any, for the current case
    f   << setw(20) << "fileUniqueName" << fileUniqueName << endl
        << setw(20) << "Comments" << comments << endl
        << endl;

    f.close();
}

void filerInfoEnd()
{
    ofstream f;
    f.open(filePath + fileUniqueName + "_" + to_string(nxy)
           + "_logCharacteristics.txt", ios_base::app);
            // open for editing

    // Write case completion data
    f   << endl << left << "Case completion data" << endl << endl
        << setw(20) << "Max. mass residual" << mChangeMax << endl
        << setw(20) << "Max. pr. residual" << pChangeMax << endl
        << setw(20) << "Max. change in u" << uChangeMax << endl
        << setw(20) << "Max. change in v" << vChangeMax << endl
        << setw(20) << "Final timestep" << t << endl
        << setw(20) << "Running time" << scriptRunningTime << " min" << endl;
}
