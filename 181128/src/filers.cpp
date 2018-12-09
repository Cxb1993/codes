/*
 * filers.cpp
 *
 *     Project: Finite Volume Navier-Stokes Solver
 *      Author: Syed Ahmad Raza
 */

#include <iomanip>      // functions for setting width
#include <fstream>      // functions for file output
#include <sstream>      // for file name conversion
#include <limits>       // for "numeric_limits<double>::digits10"
                        // it is used to set the maximum precision possible

#include "constants.h"
#include "filers.h"

using namespace std;
using namespace Numeric_lib;

/*
Function to file a 3D array. It can be used to map x, y and z direction
velocities at the correct coordinates as per staggered grid.
*/
void filer3DArray(const Matrix<double,3>& Q,
                  const Matrix<double,1>& Xp,
                  const Matrix<double,1>& Yp,
                  const Matrix<double,1>& Zp,
                  const int xIndexLimit,
                  const int yIndexLimit,
                  const int zIndexLimit,
                  const int xOffset, const int yOffset, const int zOffset,
                  const string colHeader, const string fileName)
{
    ofstream fileQ;
    fileQ.precision(numeric_limits<double>::digits10 + 2);
    fileQ.open(filePath + fileName + ".dat");

    fileQ << "TITLE= \"3D Navier-Stokes Array\"" << "\n"
          << "VARIABLES= x, y, z, " << colHeader << "\n"
          << "ZONE T= \"Single Zone\""
          << ", I=" << to_string(xIndexLimit+1)
          << ", J=" << to_string(yIndexLimit+1)
          << ", K=" << to_string(zIndexLimit+1)
          << ", F=POINT" << "\n";

    for (int k = 0; k <= zIndexLimit; ++k)
    {
        for (int j = 0; j <= yIndexLimit; ++j)
        {
            for (int i = 0; i <= xIndexLimit; ++i)
            {
                fileQ << Xp(i) << "\t" << Yp(j) << "\t" << Zp(k) << "\t"
                      << Q(i+xOffset,j+yOffset,k+zOffset) << "\n";
            }
        }
    }
}

/*
Function to file the complete solution of 3D problem, with values centered at
the cell center.
*/
void filer3DSolution(const string fileName)
{
    ofstream fileQ;
    fileQ.precision(numeric_limits<double>::digits10 + 2);
    fileQ.open(filePath + fileName + ".dat");

    fileQ << "TITLE= \"3D Navier-Stokes Solution\"" << "\n"
          << "VARIABLES= x, y, z, P, Vx, Vy, Vz, ETA, FX, FY, FZ" << "\n"
          << "ZONE T= \"Single Zone\""
          << ", I=" << to_string(nx-4)
          << ", J=" << to_string(ny-4)
          << ", K=" << to_string(nz-4)
          << ", F=POINT" << "\n";

    int offSet = ghostCells / 2;
    for (int k = 0; k <= nz-5; ++k)
    {
        for (int j = 0; j <= ny-5; ++j)
        {
            for (int i = 0; i <= nx-5; ++i)
            {
                fileQ << Xs(i) << "\t" << Ys(j) << "\t" << Zs(k) << "\t"
                      << P  (i+offSet,j+offSet,k+offSet) << "\t"
                      << U  (i+offSet,j+offSet,k+offSet) << "\t"
                      << V  (i+offSet,j+offSet,k+offSet) << "\t"
                      << W  (i+offSet,j+offSet,k+offSet) << "\t"
                      << ETA(i+offSet,j+offSet,k+offSet) << "\t"
                      << FX (i+offSet,j+offSet,k+offSet) << "\t"
                      << FY (i+offSet,j+offSet,k+offSet) << "\t"
                      << FZ (i+offSet,j+offSet,k+offSet) << "\n";
            }
        }
    }
}

/*
Function to file a one-dimensional array.
*/
void filer1(const Matrix<double,1>& R, const int nk, const string fileName)
{
    ofstream fileR;
    fileR.precision(numeric_limits<double>::digits10 + 2);
    fileR.open(filePath + fileName + ".dat");

    for (int k = 0; k < nk; ++k)
    {
        fileR << R(k) << "\n";
    }
}

/*
Function to file a two-dimensional array.
*/
void filer2(const Matrix<double,2>& T,
            const int mx, const int my, const string fileName)
{
    ofstream fileT;
    fileT.precision(numeric_limits<double>::digits10 + 2);
    fileT.open(filePath + fileName + ".dat");

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
        filer3DArray(U, Xa, Ys, Zs, nx-4, ny-5, nz-5,
                     ghostCells/2-1, ghostCells/2, ghostCells/2, "Vx",
                     fileUniqueName + "_U-final");
        filer3DArray(V, Xs, Ya, Zs, nx-5, ny-4, nz-5,
                     ghostCells/2, ghostCells/2-1, ghostCells/2, "Vy",
                     fileUniqueName + "_V-final");
        filer3DArray(W, Xs, Ys, Za, nx-5, ny-5, nz-4,
                     ghostCells/2, ghostCells/2, ghostCells/2-1, "Vz",
                     fileUniqueName + "_W-final");
        filer3DSolution(fileUniqueName + "_sol-final");
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
        filer3DArray(U, Xa, Ys, Zs, nx-4, ny-5, nz-5,
                     ghostCells/2-1, ghostCells/2, ghostCells/2, "Vx",
                     fileUniqueName + "_U-" + timeStr);
        filer3DArray(V, Xs, Ya, Zs, nx-5, ny-4, nz-5,
                     ghostCells/2, ghostCells/2-1, ghostCells/2, "Vy",
                     fileUniqueName + "_V-" + timeStr);
        filer3DArray(W, Xs, Ys, Za, nx-5, ny-5, nz-4,
                     ghostCells/2, ghostCells/2, ghostCells/2-1, "Vz",
                     fileUniqueName + "_W-" + timeStr);
        filer3DSolution(fileUniqueName + "_sol-" + timeStr);
    }
}

/*
Filing the coordinate files
*/
void filerCoordinates()
{
    filer1(Xa, (nx-3), filePath + fileUniqueName + "_coordinateX");
    filer1(Xs, (nx-4), filePath + fileUniqueName + "_coordinateXs");

    filer1(Ya, (ny-3), filePath + fileUniqueName + "_coordinateY");
    filer1(Ys, (ny-4), filePath + fileUniqueName + "_coordinateYs");

    filer1(Za, (nz-3), filePath + fileUniqueName + "_coordinateZ");
    filer1(Zs, (nz-4), filePath + fileUniqueName + "_coordinateZs");

    filer1(Dx, nx, filePath + fileUniqueName + "_coordinateDx");
    filer1(Dxs, nx, filePath + fileUniqueName + "_coordinateDxs");

    filer1(Dy, ny, filePath + fileUniqueName + "_coordinateDy");
    filer1(Dys, ny, filePath + fileUniqueName + "_coordinateDys");

    filer1(Dz, nz, filePath + fileUniqueName + "_coordinateDz");
    filer1(Dzs, nz, filePath + fileUniqueName + "_coordinateDzs");
}

void filerInfoStart()
{
    ofstream f;
    f.open(filePath + fileUniqueName + "_logCharacteristics.txt");

    // Write heading and column headers for case characteristics
    f   << left << "Finite Volume Method for 3D flow"
        << endl << endl
        << setw(20) << "Variable" << "Values" << endl
        << endl
        << setw(20) << "velScheme" << velScheme << endl
        << setw(20) << "ny" << ny << endl
        << setw(20) << "nx" << nx << endl
        << setw(20) << "nz" << nz << endl
        << setw(20) << "ghostCells" << ghostCells << endl
        << setw(20) << "nSubGrids" << nSubGrids << endl
        // << setw(20) << "CFL" << CFL << endl
        << endl

        << setw(20) << "Re" << Re << endl
        << endl

        << setw(20) << "omega" << omega << endl
        << setw(20) << "maxTimesteps" << maxTimesteps << endl
        << setw(20) << "maxPrIters" << maxPrIters << endl
        << setw(20) << "uResidual" << uResidual << endl
        << setw(20) << "vResidual" << vResidual << endl
        << setw(20) << "wResidual" << wResidual << endl
        << setw(20) << "pResidual" << pResidual << endl
        << setw(20) << "dt" << dt << endl
        << setw(20) << "numOfThreads" << numOfThreads << endl
        << endl

        << setw(20) << "L" << L << endl
        << setw(20) << "B" << B << endl
        << setw(20) << "H" << H << endl
        << setw(20) << "uIn" << uIn << endl
        << setw(20) << "vIn" << vIn << endl
        << setw(20) << "wIn" << wIn << endl
        << setw(20) << "uInitial" << uInitial << endl
        << setw(20) << "vInitial" << vInitial << endl
        << setw(20) << "wInitial" << wInitial << endl
        << setw(20) << "pInitial" << pInitial << endl
        // << setw(20) << "pIn" << pIn << " (unused for Neumann)" << endl
        // << setw(20) << "pout" << pout << " (unused for Neumann)" << endl
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
    f.open(filePath + fileUniqueName
        + "_logCharacteristics.txt", ios_base::app);
            // open for editing

    // Write case completion data
    f   << endl << left << "Case completion data" << endl << endl
        << setw(20) << "Max. mass residual" << mChangeMax << endl
        << setw(20) << "Max. pr. residual" << pChangeMax << endl
        << setw(20) << "Max. change in u" << uChangeMax << endl
        << setw(20) << "Max. change in v" << vChangeMax << endl
        << setw(20) << "Max. change in w" << wChangeMax << endl
        << setw(20) << "tPrConvrg" << tPrConvrg
                    <<  "\t\t\t// 1st timestep when pIter < maxPrIters" << endl
        << setw(20) << "Final timestep" << t << endl
        << setw(20) << "Running time" << scriptRunningTime << " min" << endl;
}
