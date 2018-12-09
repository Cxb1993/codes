/*
 * filename.cpp
 *
 *  Created on: 2017-04-17
 *      Author: Syed Ahmad Raza
 */

#include <cmath>        // math functions

#include "constants.h"
#include "gridder.h"
#include "filers.h"
#include "navierFVD.h"

using namespace std;

double X[nx] = {0.0};        // x coordinates
double Y[ny] = {0.0};        // y coordinates

void gridder()
{
    for (int i = 0; i < nx; ++i)
    {
        // X[i] = (0.5 * lx) * (1 + sin( pi * ((i * 1.0 / (nx-1)) - 0.5)));
        X[i] = static_cast<double>(i)*dx;
    }
    for (int j = 0; j < ny; ++j)
    {
        // Y[j] = ly * (sin( pi * (j * 1.0 / (ny-1)) * 0.5));
        Y[j] = static_cast<double>(j)*dy;
    }
}
