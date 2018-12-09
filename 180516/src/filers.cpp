/*
 * filers.cpp
 *
 *     Project: Finite Volume Navier-Stokes Solver
 *      Author: Syed Ahmad Raza
 */

#include <iomanip>      // functions for setting width
#include <fstream>      // functions for file output
#include <sstream>      // for file name conversion
// #include <string>       // for fileName - already in constants.h
#include <limits>       // for "numeric_limits<double>::digits10"
                        // it is used to set the maximum precision possible

#include "constants.h"
#include "filers.h"

using namespace std;
using namespace Numeric_lib;

ofstream progressFile;

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
    /*
    Function to file a 3D array. It can be used to map x, y and z direction
    velocities at the correct coordinates as per staggered grid.
    */
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

void filer3DSolution(const string fileName)
{
    /*
    Function to file the complete solution of 3D problem.
    */
    ofstream fileQ;
    fileQ.precision(numeric_limits<double>::digits10 + 2);
    fileQ.open(filePath + fileName + ".dat");

    fileQ << "TITLE= \"2D Navier-Stokes Solution\"" << "\n"
          << "VARIABLES= x, y, z, P, Vx, Vy, Vz" << "\n"
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
                      << P(i+offSet,j+offSet,k+offSet) << "\t"
                      << U(i+offSet,j+offSet,k+offSet) << "\t"
                      << V(i+offSet,j+offSet,k+offSet) << "\t"
                      << W(i+offSet,j+offSet,k+offSet) << "\n";
            }
        }
    }
}

// void filer2DArray(const Matrix<double,2>& Q,
//                   const Matrix<double,1>& Xp,
//                   const Matrix<double,1>& Yp,
//                   const int xIndexLimit, const int yIndexLimit,
//                   const int xOffset, const int yOffset,
//                   const string colHeader, const string fileName)
// {
//     /*
//     Function to file a 2D array. It can be used to map x and y direction
//     velocities at the correct coordinates as per staggered grid.
//     */
//     ofstream fileQ;
//     fileQ.precision(numeric_limits<double>::digits10 + 2);
//     fileQ.open(filePath + fileName + ".dat");
//
//     fileQ << "TITLE= \"2D Navier-Stokes Array\"" << "\n"
//           << "VARIABLES= x, y, z, " << colHeader << "\n"
//           << "ZONE T= \"Single Zone\""
//           << ", I=" << to_string(xIndexLimit+1)
//           << ", J=" << to_string(yIndexLimit+1)
//           << ", K=1"
//           << ", F=POINT" << "\n";
//
//     for (int j = 0; j <= yIndexLimit; ++j)
//     {
//         for (int i = 0; i <= xIndexLimit; ++i)
//         {
//             fileQ << Xp(i) << "\t" << Yp(j) << "\t" << "0.0" << "\t"
//                   << Q(i+xOffset,j+yOffset) << "\n";
//         }
//     }
// }

// void filer2DSolution(const string fileName)
// {
//     /*
//     Function to file the complete solution of 2D problem.
//     */
//     ofstream fileQ;
//     fileQ.precision(numeric_limits<double>::digits10 + 2);
//     fileQ.open(filePath + fileName + ".dat");
//
//     fileQ << "TITLE= \"2D Navier-Stokes Solution\"" << "\n"
//           << "VARIABLES= x, y, z, P, Vx, Vy" << "\n"
//           << "ZONE T= \"Single Zone\""
//           << ", I=" << to_string(nx-4)
//           << ", J=" << to_string(ny-4)
//           << ", K=1"
//           << ", F=POINT" << "\n";
//
//     int offSet = ghostCells / 2;
//     for (int j = 0; j <= ny-5; ++j)
//     {
//         for (int i = 0; i <= nx-5; ++i)
//         {
//             fileQ << Xs(i) << "\t" << Ys(j) << "\t" << 0.0 << "\t"
//                   << P(i+offSet,j+offSet) << "\t"
//                   << U(i+offSet,j+offSet) << "\t"
//                   << V(i+offSet,j+offSet) << "\n";
//         }
//     }
// }

void filer1(const Matrix<double,1>& R, const int nk, const string fileName)
/*
Function to file a one-dimensional array.
*/
{
    ofstream fileR;
    fileR.precision(numeric_limits<double>::digits10 + 2);
    fileR.open(filePath + fileName + ".dat");

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
        // filer2(U, nx, ny,
        //         fileUniqueName + "_" + to_string(nxyz) + "_U-final");
        // filer2(V, nx, ny,
        //         fileUniqueName + "_" + to_string(nxyz) + "_V-final");
        // filer2(P, nx, ny,
        //         fileUniqueName + "_" + to_string(nxyz) + "_P-final");
        // filer2(MC, nx, ny,
        //         fileUniqueName + "_" + to_string(nxyz) + "_MC-final");
        // filer2(PC, nx, ny,
        //         fileUniqueName + "_" + to_string(nxyz) + "_PC-final");

        // filer2DArray(U, Xa, Ys, nx-4, ny-5,
        //              ghostCells/2-1, ghostCells/2, "Vx",
        //              fileUniqueName + "_" + to_string(nxyz) + "_U-final");
        // filer2DArray(V, Xs, Ya, nx-5, ny-4,
        //              ghostCells/2, ghostCells/2-1, "Vy",
        //              fileUniqueName + "_" + to_string(nxyz) + "_V-final");
        // filer2DSolution(fileUniqueName + "_" + to_string(nxyz)
        //                 + "_sol-final");

        filer3DArray(U, Xa, Ys, Zs, nx-4, ny-5, nz-5,
                     ghostCells/2-1, ghostCells/2, ghostCells/2, "Vx",
                     fileUniqueName + "_" + to_string(nxyz) + "_U-final");
        filer3DArray(V, Xs, Ya, Zs, nx-5, ny-4, nz-5,
                     ghostCells/2, ghostCells/2-1, ghostCells/2, "Vy",
                     fileUniqueName + "_" + to_string(nxyz) + "_V-final");
        filer3DArray(W, Xs, Ys, Za, nx-5, ny-5, nz-4,
                     ghostCells/2, ghostCells/2, ghostCells/2-1, "Vz",
                     fileUniqueName + "_" + to_string(nxyz) + "_W-final");
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
        // filer2(U, nx, ny,
        //         fileUniqueName + "_" + to_string(nxyz) + "_U-" + timeStr);
        // filer2(V, nx, ny,
        //         fileUniqueName + "_" + to_string(nxyz) + "_V-" + timeStr);
        // filer2(P, nx, ny,
        //         fileUniqueName + "_" + to_string(nxyz) + "_P-" + timeStr);
        // filer2(MC, nx, ny,
        //         fileUniqueName + "_" + to_string(nxyz) + velScheme
        //         + "_MC-" + timeStr);
        // filer2(PC, nx, ny,
        //         fileUniqueName + "_" + to_string(nxyz) + velScheme
        //         + "_PC-" + timeStr);

        // filer2DArray(U, Xa, Ys, nx-4, ny-5,
        //              ghostCells/2-1, ghostCells/2, "Vx",
        //              fileUniqueName + "_" + to_string(nxyz) + "_U-"
        //              + timeStr);
        // filer2DArray(V, Xs, Ya, nx-5, ny-4,
        //              ghostCells/2, ghostCells/2-1, "Vy",
        //              fileUniqueName + "_" + to_string(nxyz) + "_V-"
        //              + timeStr);
        // filer2DSolution(fileUniqueName + "_" + to_string(nxyz) + "_sol-"
        //                 + timeStr);

        filer3DArray(U, Xa, Ys, Zs, nx-4, ny-5, nz-5,
                     ghostCells/2-1, ghostCells/2, ghostCells/2, "Vx",
                     fileUniqueName + "_" + to_string(nxyz) + "_U-"
                     + timeStr);
        filer3DArray(V, Xs, Ya, Zs, nx-5, ny-4, nz-5,
                     ghostCells/2, ghostCells/2-1, ghostCells/2, "Vy",
                     fileUniqueName + "_" + to_string(nxyz) + "_V-"
                     + timeStr);
        filer3DArray(W, Xs, Ys, Za, nx-5, ny-5, nz-4,
                     ghostCells/2, ghostCells/2, ghostCells/2-1, "Vz",
                     fileUniqueName + "_" + to_string(nxyz) + "_W-"
                     + timeStr);
        filer3DSolution(fileUniqueName + "_" + to_string(nxyz)
                        + "_sol-final");
    }
}

void filerCoordinates()
{
    /*
    Filing the coordinate files
    */

    filer1(Xa, (nx-3), filePath + fileUniqueName + "_" + to_string(nxyz)
            + "_coordinateX");
    filer1(Xs, (nx-4), filePath + fileUniqueName + "_" + to_string(nxyz)
            + "_coordinateXs");

    filer1(Ya, (ny-3), filePath + fileUniqueName + "_" + to_string(nxyz)
            + "_coordinateY");
    filer1(Ys, (ny-4), filePath + fileUniqueName + "_" + to_string(nxyz)
            + "_coordinateYs");

    filer1(Za, (nz-3), filePath + fileUniqueName + "_" + to_string(nxyz)
            + "_coordinateZ");
    filer1(Zs, (nz-4), filePath + fileUniqueName + "_" + to_string(nxyz)
            + "_coordinateZs");

    filer1(Dx, nx, filePath + fileUniqueName + "_" + to_string(nxyz)
            + "_coordinateDx");
    filer1(Dxs, nx, filePath + fileUniqueName + "_" + to_string(nxyz)
            + "_coordinateDxs");

    filer1(Dy, ny, filePath + fileUniqueName + "_" + to_string(nxyz)
            + "_coordinateDy");
    filer1(Dys, ny, filePath + fileUniqueName + "_" + to_string(nxyz)
            + "_coordinateDys");

    filer1(Dz, nz, filePath + fileUniqueName + "_" + to_string(nxyz)
            + "_coordinateDz");
    filer1(Dzs, nz, filePath + fileUniqueName + "_" + to_string(nxyz)
            + "_coordinateDzs");
}

void filerCreateProgress()
{
    progressFile.open(filePath + fileUniqueName + "_" + to_string(nxyz)
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
            << setw(5)  << "i" << setw(2) << ", "
            << setw(4)  << "j" << setw(2) << ", "
            << setw(4)  << "k"
            << setw(5)  << "i" << setw(2) << ", "
            << setw(4)  << "j" << setw(2) << ", "
            << setw(4)  << "k"
            << setw(5)  << "i" << setw(2) << ", "
            << setw(4)  << "j" << setw(2) << ", "
            << setw(4)  << "k"
            << setw(5)  << "i" << setw(2) << ", "
            << setw(4)  << "j" << setw(2) << ", "
            << setw(4)  << "k"
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
        << setw(17) << wChangeMax
        << endl
        << setw(12) << ""
        << setw(12) << ""
        << setw(5)  << mChangeMax_i << setw(2) << ", "
        << setw(4)  << mChangeMax_j << setw(2) << ", "
        << setw(4)  << mChangeMax_k
        << setw(5)  << pChangeMax_i << setw(2) << ", "
        << setw(4)  << pChangeMax_j << setw(2) << ", "
        << setw(4)  << pChangeMax_k
        << setw(5)  << uChangeMax_i << setw(2) << ", "
        << setw(4)  << uChangeMax_j << setw(2) << ", "
        << setw(4)  << uChangeMax_k
        << setw(5)  << vChangeMax_i << setw(2) << ", "
        << setw(4)  << vChangeMax_j << setw(2) << ", "
        << setw(4)  << vChangeMax_k
        << setw(5)  << wChangeMax_i << setw(2) << ", "
        << setw(4)  << wChangeMax_j << setw(2) << ", "
        << setw(4)  << wChangeMax_k
        << endl;
}

void filerInfoStart()
{
    ofstream f;
    f.open(filePath + fileUniqueName + "_" + to_string(nxyz)
           + "_logCharacteristics.txt");

    // Write heading and column headers for case characteristics
    f   << left << "Finite Volume Method for 2D flow"
        << endl << endl
        << setw(20) << "Variable" << "Values" << endl
        << endl
        << setw(20) << "velScheme" << velScheme << endl
        << setw(20) << "ny" << ny << endl
        << setw(20) << "nx" << nx << endl
        << setw(20) << "nz" << nz << endl
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
        << setw(20) << "wResidual" << wResidual << endl
        << setw(20) << "pResidual" << pResidual << endl
        << setw(20) << "dt" << dt << " s" << endl
        << endl

        << setw(20) << "x_y_ratio" << x_y_ratio << endl
        << setw(20) << "x_z_ratio" << x_z_ratio << endl
        << setw(20) << "L" << L << " m" << endl
        << setw(20) << "B" << B << " m" << endl
        << setw(20) << "H" << H << " m" << endl
        << setw(20) << "uIn" << uIn << " m/s" << endl
        << setw(20) << "vIn" << vIn << " m/s" << endl
        << setw(20) << "uInitial" << uInitial << " m/s" << endl
        << setw(20) << "vInitial" << vInitial << " m/s" << endl
        << setw(20) << "wInitial" << wInitial << " m/s" << endl
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
    f.open(filePath + fileUniqueName + "_" + to_string(nxyz)
           + "_logCharacteristics.txt", ios_base::app);
            // open for editing

    // Write case completion data
    f   << endl << left << "Case completion data" << endl << endl
        << setw(20) << "Max. mass residual" << mChangeMax << endl
        << setw(20) << "Max. pr. residual" << pChangeMax << endl
        << setw(20) << "Max. change in u" << uChangeMax << endl
        << setw(20) << "Max. change in v" << vChangeMax << endl
        << setw(20) << "Max. change in w" << wChangeMax << endl
        << setw(20) << "Final timestep" << t << endl
        << setw(20) << "Running time" << scriptRunningTime << " min" << endl;
}
