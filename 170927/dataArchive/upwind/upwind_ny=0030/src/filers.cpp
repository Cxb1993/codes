/*
 * filers.cpp
 *
 *  Created on: 2017-09-27
 *      Author: Syed Ahmad Raza
 */

#include <fstream>      // functions for file output
#include <string>       // for file name conversion (to_string)
#include <sstream>      // for file name conversion

#include "constants.h"
#include "filers.h"

using namespace std;

// Function to output one-dimensional array to a folder named 'data' in the parent directory
void filer1(double Z[], int nk, string fileName)
{
    ofstream fileZ;
    fileZ.open("../data/" + fileName + ".dat");
    for (int k = 0; k < nk; ++k)
    {
        fileZ << Z[k] << '\n';
    }
}

// Function to output two-dimensional array to a folder named 'data' in the parent directory
void filer2(double Q[][ny], int nx, int ny, string fileName)
{
    ofstream fileQ;
    fileQ.open("../data/" + fileName + ".dat");
    for (int j = 0; j < ny; ++j)
    {
        for (int i = 0; i < nx; ++i)
        {
            fileQ << Q[i][j] << '\t';
        }
        fileQ << '\n';
    }
}

void filerAllSol(string fileNamePrefix, int timeStep)
{
    if (timeStep == -1)
    {
        filer2(U, nx, ny, fileNamePrefix + "_U_Final");
        filer2(V, nx, ny, fileNamePrefix + "_V_Final");
        filer2(P, nx, ny, fileNamePrefix + "_P_Final");
    }
    else
    {
        // string timeStr = to_string(timeStep);    // does not work in some compilers
        ostringstream timeStrTemp;
        timeStrTemp << timeStep;
        string timeStr = timeStrTemp.str();
        filer2(U, nx, ny, fileNamePrefix + "_U_" + timeStr);
        filer2(V, nx, ny, fileNamePrefix + "_V_" + timeStr);
        filer2(P, nx, ny, fileNamePrefix + "_P_" + timeStr);
    }
}
