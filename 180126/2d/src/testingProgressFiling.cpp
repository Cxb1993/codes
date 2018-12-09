/*
 * testingProgressFiling.cpp
 *
 *  Created on: 2018-03-06
 *      Author: Syed Ahmad Raza
 *
 * Testing how to store console output to a file
 */

 #include <iostream>     // functions for input and output to console
 #include <iomanip>      // functions for setting width
 #include <fstream>      // functions for file output
 #include <sstream>      // for file name conversion
 #include <string>       // for fileName
 #include <limits>       // for "numeric_limits<double>::digits10"
                         // it is used to set the maximum precision possible

using namespace std;

int t = 1;
ofstream progressFile;

void progressFileCreator()
{
    progressFile.open("progressLog.dat");
    progressFile << "Initialized";
    // progressFile.close();
}

void progressFiler()
{
    // progressFile.open("progressLog.dat", ios_base::app);    // open for editing
    // Write column headers
    if (t == 1)
    {
        progressFile << setw(20) << "Timestep"
                     << setw(20) << "pIter"
                     << setw(20) << "Max. pChange"
                     << setw(20) << "Max. uChange" << endl << endl;
    }
    // Write timestep data
    progressFile << setw(20) << t << endl;
}

int main()
{
    progressFileCreator();
    // progressFile.open("progressLog.dat");
    for (t = 1; t < 10; ++t)
    {
        progressFiler();
    }
    return 0;
}