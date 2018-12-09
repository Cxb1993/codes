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

ofstream CDFile;
ofstream CLFile;

/*
Function to file a 3D array. It can be used to map x, y and z direction
velocities at the correct coordinates as per staggered grid.
*/
void filer3DArray(const Matrix<double,2>& Q,
                  const Matrix<double,1>& Xp,
                  const Matrix<double,1>& Yp,

                  const int xIndexLimit,
                  const int yIndexLimit,

                  const int xOffset, const int yOffset,
                  const string colHeader, const string fileName)
{
    ofstream fileQ;
    fileQ.precision(numeric_limits<double>::digits10 + 2);
    fileQ.open(filePath + fileName + ".dat");

    fileQ << "TITLE= \"3D Navier-Stokes Array\"" << "\n"
          << "VARIABLES= x, y, " << colHeader << "\n"
          << "ZONE T= \"Single Zone\""
          << ", I=" << to_string(xIndexLimit+1)
          << ", J=" << to_string(yIndexLimit+1)
          << ", F=POINT" << "\n";


        for (int j = 0; j <= yIndexLimit; ++j)
        {
            for (int i = 0; i <= xIndexLimit; ++i)
            {
                fileQ << Xp(i) << "\t" << Yp(j) << "\t"
                      << Q(i+xOffset,j+yOffset) << "\n";

            }
        }

}

/*
Function to file the complete solution of 3D problem.
*/
void filer3DSolution(const string fileName)
{
    ofstream fileQ;
    fileQ.precision(numeric_limits<double>::digits10 + 2);
    fileQ.open(filePath + fileName + ".dat");

    fileQ << "TITLE= \"2D Navier-Stokes Solution\"" << "\n"
          << "VARIABLES= x, y, P, Vx, Vy, VF_X, VF_Y,eta, STREAM, VORT" << "\n"
          << "ZONE T= \"Single Zone\""
          << ", I=" << to_string(nx-4)
          << ", J=" << to_string(ny-4)
          << ", F=POINT" << "\n";

    int offSet = ghostCells / 2;
    for (int j = 0; j <= ny-5; ++j)
        {
            for (int i = 0; i <= nx-5; ++i)
            {
                fileQ << Xs(i) << "\t" << Ys(j) << "\t"
                      << P(i+offSet,j+offSet) << "\t"
                      << U(i+offSet,j+offSet) << "\t"
                      << V(i+offSet,j+offSet) << "\t"
                      << VF_X(i+offSet,j+offSet) << "\t"
                      << VF_Y(i+offSet,j+offSet) << "\t"
                      << Eta(i+offSet,j+offSet) << "\t"
                      << STR (i+offSet,j+offSet) << "\t"
                      << VORT(i+offSet,j+offSet) << "\n";
            }
        }
        fileQ //<< "CD ="<< CD << "\n"
          //<< "CL ="<< CL << "\n"
          << "CD ="<< CD_Temp << "\n"
          << "CL ="<< CL_Temp << "\n";
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
        filer3DArray(U, Xa, Ys, nx-4, ny-5,
                     ghostCells/2-1, ghostCells/2, "Vx",
                     fileUniqueName + "_" + to_string(nxyz) + "_U-final");
        filer3DArray(V, Xs, Ya, nx-5, ny-4,
                     ghostCells/2, ghostCells/2-1, "Vy",
                     fileUniqueName + "_" + to_string(nxyz) + "_V-final");
        filer3DSolution(fileUniqueName + "_" + to_string(nxyz)
                        + "_sol-final");
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
        filer3DArray(U, Xa, Ys, nx-4, ny-5,
                     ghostCells/2-1, ghostCells/2, "Vx",
                     fileUniqueName + "_" + to_string(nxyz) + "_U-" + timeStr);
        filer3DArray(V, Xs, Ya, nx-5, ny-4,
                     ghostCells/2, ghostCells/2-1,"Vy",
                     fileUniqueName + "_" + to_string(nxyz) + "_V-" + timeStr);
        filer3DSolution(fileUniqueName + "_" + to_string(nxyz) + "_sol-"
                        + timeStr);
    }
}

/*
Filing the coordinate files
*/
void filerCoordinates()
{
    filer1(Xa, (nx-3), filePath + fileUniqueName + "_" + to_string(nxyz)
            + "_coordinateX");
    filer1(Xs, (nx-4), filePath + fileUniqueName + "_" + to_string(nxyz)
            + "_coordinateXs");
    filer1(Ya, (ny-3), filePath + fileUniqueName + "_" + to_string(nxyz)
            + "_coordinateY");
    filer1(Ys, (ny-4), filePath + fileUniqueName + "_" + to_string(nxyz)
            + "_coordinateYs");
    filer1(Dx, nx, filePath + fileUniqueName + "_" + to_string(nxyz)
            + "_coordinateDx");
    filer1(Dxs, nx, filePath + fileUniqueName + "_" + to_string(nxyz)
            + "_coordinateDxs");

    filer1(Dy, ny, filePath + fileUniqueName + "_" + to_string(nxyz)
            + "_coordinateDy");
    filer1(Dys, ny, filePath + fileUniqueName + "_" + to_string(nxyz)
            + "_coordinateDys");
}

void filerInfoStart()
{
    ofstream f;
    f.open(filePath + fileUniqueName + "_" + to_string(nxyz)
           + "_logCharacteristics.txt");

    // Write heading and column headers for case characteristics
    f   << left << "Finite Volume Method for 3D flow"
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
        << setw(20) << "maxPrIters" << maxPrIters << endl
        << setw(20) << "uResidual" << uResidual << endl
        << setw(20) << "vResidual" << vResidual << endl

        << setw(20) << "pResidual" << pResidual << endl
        << setw(20) << "dt" << dt << " s" << endl
        << setw(20) << "numOfThreads" << numOfThreads << endl
        << endl

        << setw(20) << "x_y_ratio" << x_y_ratio << endl
        << setw(20) << "x_z_ratio" << x_z_ratio << endl
        << setw(20) << "L" << L << " m" << endl
        << setw(20) << "B" << B << " m" << endl
        //<< setw(20) << "H" << H << " m" << endl
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

void filerCreateCD()
{
    CDFile.open(filePath + fileUniqueName + "_" + to_string(nxyz)
                      + "_CD.dat");
}

void filerCD()
{
    // Write timestep data
    CDFile << (t*dt)
        << "\t"<< CD_Temp
        << endl;
}

void filerCreateCL()
{
    CLFile.open(filePath + fileUniqueName + "_" + to_string(nxyz)
                      + "_CL.dat");
}

void filerCL()
{
    // Write timestep data
    CLFile << (t*dt)
        << "\t" << CL_Temp
        << endl;
}

void filerInfoEnd()
{
    ofstream f;
    f.open(filePath + fileUniqueName + "_" + to_string(nxyz)
           + "_logCharacteristics.txt", ios_base::app);
            // open for editing

    // Write case completion data
    f   << endl << left << "Case completion data" << endl << endl
        << setw(20) << "Max. mass residual" << mChangeMax << endl
        << setw(20) << "Max. pr. residual" << pChangeMax << endl
        << setw(20) << "Max. change in u" << uChangeMax << endl
        << setw(20) << "Max. change in v" << vChangeMax << endl
        << setw(20) << "CD" << CD_Temp << " " << endl
        << setw(20) << "CL" << CL_Temp << " " << endl
        //<< setw(20) << "Max. change in w" << wChangeMax << endl
        << setw(20) << "tPrConvrg" << tPrConvrg
                    <<  "\t\t\t// 1st timestep when pIter < maxPrIters" << endl
        << setw(20) << "Final timestep" << t << endl
        << setw(20) << "Running time" << scriptRunningTime << " min" << endl;
}
