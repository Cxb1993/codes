/*
 * filename.cpp
 *
 *  Created on: 2017-09-27
 *      Author: Syed Ahmad Raza
 */

#include <cmath>        // math functions

#include "constants.h"
#include "filers.h"
#include "gridder.h"

using namespace std;

double X[2*nx-1] = {0.0};        // x coordinates
double Y[2*ny-1] = {0.0};        // y coordinates
// double XS[nx] = {0.0};        // staggered x coordinates
// double YS[ny] = {0.0};        // staggered y coordinates

void gridder()
{
    for (int i = 0; i < (2*nx-1); ++i)
    {
        X[i] = (0.5 * lx) * (1 + sin( pi * ((i * 1.0 / (nx-1)) - 0.5)));;
        // X[i] = static_cast<double>(i) * dx / 2;
    }
    for (int j = 0; j < (2*ny-1); ++j)
    {
        Y[j] = ly * (sin( pi * (j * 1.0 / (ny-1)) * 0.5));
        // Y[j] = static_cast<double>(j) * dy / 2;
    }
    filer1(X, (2*nx-1), filePath + "coordinateX_ny-" + to_string(ny));
    filer1(Y, (2*ny-1), filePath + "coordinateY_ny-" + to_string(ny));
}
