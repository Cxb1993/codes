/*
 * filer.cpp
 *
 *  Created on: 2017-04-20
 *      Author: Syed Ahmad Raza
 */

#include <fstream>      // functions for file input and output
#include <string>       // for file name

#include "constants.h"
#include "gridder.h"
#include "filers.h"
#include "navierFVD.h"

using namespace std;

void gridFiler(double Z[], int nk, string fileName)
{
    ofstream fileZ;
    fileZ.open(fileName);
    for (int k = 0; k < nk; ++k)
    {
        fileZ << Z[k] << '\n';
    }
}

void solFiler(double Q[][ny], int nx, int ny, string fileName)
{
    ofstream file;
    file.open(fileName);
    for (int j = 0; j < ny; ++j)
    {
        for (int i = 0; i < nx; ++i)
        {
            file << Q[i][j] << '\t';
        }
        file << '\n';
    }
}
