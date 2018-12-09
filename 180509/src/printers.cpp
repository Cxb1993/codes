/*
 * printers.cpp
 *
 *  Created on: 2017-12-26
 *      Author: Syed Ahmad Raza
 */

 #include <iostream>     // functions for input and output to console
 #include <iomanip>      // functions for setting width
 #include <fstream>      // functions for file output
 #include <sstream>      // for file name conversion
 #include <limits>       // for "numeric_limits<double>::digits10"
                         // it is used to set the maximum precision possible

 #include "constants.h"
 #include "printers.h"

using namespace std;

void printerProgress()
{
    // Print column headers
    if (t == 1)
    {
        cout
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
    // Print timestep data
    cout
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

void printerInfoEnd()
{
    // Write heading and column headers for case characteristics
    cout << left << "Finite Volume Method for 2D flow"
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
    cout << setw(20) << "fileUniqueName" << fileUniqueName << endl
         << setw(20) << "Comments" << comments << endl
         << endl;

    // Print case completion data
    cout << endl << "Case completion data" << endl << endl
         << setw(20) << "Max. mass residual" << mChangeMax << endl
         << setw(20) << "Max. pr. residual" << pChangeMax << endl
         << setw(20) << "Max. change in u" << uChangeMax << endl
         << setw(20) << "Max. change in v" << vChangeMax << endl
         << setw(20) << "Final timestep" << t << endl
         << setw(20) << "Running time" << scriptRunningTime << " min" << endl;
}
