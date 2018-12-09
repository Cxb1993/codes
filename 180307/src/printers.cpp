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
 #include <string>       // for fileName
 #include <limits>       // for "numeric_limits<double>::digits10"
                         // it is used to set the maximum precision possible

 #include "constants.h"
 #include "printers.h"

using namespace std;

void progressPrinter()
{
    // Print column headers
    if (t == 1)
    {
        cout << setw(12) << "Timestep"
             << setw(12) << "P iters"
             << setw(17) << "Max. mChange"
             << setw(17) << "Max. pChange"
             << setw(17) << "Max. uChange"
             << setw(17) << "Max. vChange"
             << endl;
    }
    // Print timestep data
    cout << setw(12) << t
         << setw(12) << pIter
         << setw(17) << mChangeMax
         << setw(17) << pChangeMax
         << setw(17) << uChangeMax
         << setw(17) << vChangeMax
         << endl;
}

void infoPrinter()
{
    // Print final timestep
    cout << endl << left << setw(15) << "Final timestep" << t << endl << endl;

    // // Print analytical data
    // cout << left << "Analytical data\n";
    // cout << setw(15) << "uAvg" << uAvg << " m/s\n";
    // cout << setw(15) << "pavg_inlet" << pavg_inlet << " Pa\n";
    // cout << setw(15) << "pavg_outlet" << pavg_outlet << " Pa\n";
    // cout << setw(15) << "deltaP" << deltaP << " Pa\n\n";

    // Print case characteristics
    cout << "Case characteristics\n";
    cout << setw(15) << "velScheme" << velScheme << "\n";
    cout << setw(15) << "Re" << Re << "\n";
    cout << setw(15) << "entryLength" << entryLength << " m\n";
    cout << setw(15) << "CFL" << CFL << "\n\n";

    streamsize sc = cout.precision();   // store current precision value
    cout.precision(numeric_limits<double>::digits10 + 2);
    cout << setw(15) << "l2Norm" << l2Norm << "\n\n";
    cout.precision(sc);                 // revert to previous precision value
    cout << setw(15) << "Running time" << scriptRunningTime << " minutes\n\n";

    if (L <= entryLength)
    {
        cout << "WARNING: For 2D channel flow boundary layer is not fully\
        developed.\n\n";
    }
}
