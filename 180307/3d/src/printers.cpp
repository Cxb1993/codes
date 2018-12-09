/*
 * printers.cpp
 *
 *  Created on: 2018-01-26
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
        printf("\n%9s %7s %16s \t %16s\n", "Timestep", "P iters",
               "Max. pChange", "Max. uChange");
    }
    // Print timestep data
    printf("%9d %7d %16f \t %16e\n", t, pIter, pChange, uChange);
}

void infoPrinter()
{
    // Print final timestep
    cout << endl << setw(15) << left << "Final timestep" << t << "\n\n";

    // Print analytical data
    // cout << "Analytical data\n";
    // cout << setw(15) << left << "uAvg" << uAvg << "\n";
    // cout << setw(15) << left << "pavg_inlet" << pavg_inlet << "\n";
    // cout << setw(15) << left << "pavg_outlet" << pavg_outlet << "\n";
    // cout << setw(15) << left << "deltaP" << deltaP << "\n\n";

    // Print case characteristics
    cout << "Case characteristics\n";
    cout << setw(15) << left << "velScheme" << velScheme << "\n";
    cout << setw(15) << left << "Re" << Re << "\n";
    cout << setw(15) << left << "entryLength" << entryLength << "\n";
    cout << setw(15) << left << "CFL" << CFL << "\n\n";

    // streamsize sc = cout.precision();   // store current precision value
    // cout.precision(numeric_limits<double>::digits10 + 2);
    // cout << setw(15) << left << "l2Norm" << l2Norm << "\n\n";
    // cout.precision(sc);                 // revert to previous precision value

    if (L <= entryLength)
    {
        cout << "WARNING: Boundary Layer is not fully developed.\n\n";
    }
}
