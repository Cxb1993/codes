/*
 * printers.cpp
 *
 *     Project: Finite Volume Navier-Stokes Solver
 *      Author: Syed Ahmad Raza
 */

 #include <iostream>     // functions for input and output to console
 #include <iomanip>      // functions for setting width
 #include <sstream>      // for file name conversion
 #include <limits>       // for "numeric_limits<double>::digits10"
                         // it is used to set the maximum precision possible

 #include "constants.h"
 #include "printers.h"

using namespace std;

void printerProgress()
{
    // Print column headers
    if (t == 1 || t % 50 == 0)
    {
        cout
            << endl
            << setw(12) << "Timestep"
            << setw(12) << "P iters"
            << setw(17) << "Max. mChange"
            << setw(17) << "Max. pChange"
            << setw(17) << "Max. uChange"
            << setw(17) << "Max. vChange"
           // << setw(17) << "Max. wChange"
            << endl << endl;
    }
    // Print timestep data
    cout
        << setw(12) << t
        << setw(12) << pIter
        << setw(17) << mChangeMax
        << setw(17) << pChangeMax
        << setw(17) << uChangeMax
        << setw(17) << vChangeMax
        //<< setw(17) << wChangeMax
        << endl;
}

void printerInfoEnd()
{
    // Write heading and column headers for case characteristics
    cout << left << "Finite Volume Method for 3D flow"
         << endl << endl
         << setw(20) << "Variable" << "Values" << endl
         << endl
         << setw(20) << "velScheme" << velScheme << endl
         << setw(20) << "nx" << nx << endl
         << setw(20) << "ny" << ny << endl
         //<< setw(20) << "nz" << nz << endl
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
        // << setw(20) << "wResidual" << wResidual << endl
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
         << setw(20) << "wIn" << wIn << " m/s" << endl
         << setw(20) << "uInitial" << uInitial << " m/s" << endl
         << setw(20) << "vInitial" << vInitial << " m/s" << endl
         //<< setw(20) << "wInitial" << wInitial << " m/s" << endl
         << setw(20) << "pInitial" << pInitial << " Pa" << endl
         // << setw(20) << "pIn" << pIn << " Pa (unused for Neumann)" << endl
         // << setw(20) << "pout" << pout << " Pa (unused for Neumann)" << endl
        << endl;

    // Write comments, if any, for the current case
    cout << setw(20) << "fileUniqueName" << fileUniqueName << endl
         << setw(20) << "Comments" << comments << endl
         << endl;

    // Print case completion data
    cout << endl << "Case completion data" << endl << endl
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
         << setw(20) << "Running time" << scriptRunningTime << " minutes"
                     << endl;
}
