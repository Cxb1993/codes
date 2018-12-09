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

double X[nx-1]      = {0.0};        // x coordinates
double Y[ny-1]      = {0.0};        // y coordinates
double Dx[nx-2]     = {0.0};
double Dx1[nx-2]    = {0.0};
double Dx2[nx-2]    = {0.0};
double Dy[ny-2]     = {0.0};
double Dy1[ny-2]    = {0.0};
double Dy2[ny-2]    = {0.0};

void gridder()
{
    for (int i = 0; i < (nx-1); ++i)
    {
        // X[i] = (0.5 * lx) * (1 + sin( pi *
        //     ((static_cast<double>(i) / static_cast<double>(nx-1)) - 0.5)));
        // X[i] = lx *
        // (sin( pi * (static_cast<double>(i) / static_cast<double>(ny-1)) * 0.5));
        X[i] = static_cast<double>(i) * L / static_cast<double>(nx-1);
    }
    for (int j = 0; j < (ny-1); ++j)
    {
        // Y[j] = ly *
        // (sin( pi * (static_cast<double>(j) / static_cast<double>(ny-1)) * 0.5));
        Y[j] = static_cast<double>(j) * D / static_cast<double>(ny-1);
    }
    filer1(X, (nx-1), filePath + "coordinateX_ny-" + to_string(ny));
    filer1(Y, (ny-1), filePath + "coordinateY_ny-" + to_string(ny));

    for (int i = 0; i < (nx-2); ++i)
    {
        if (i == 0)
        {
            Dx[i] = X[i+1] - X[i];
            Dx1[i] = (X[i+2] - X[i]) / 2.0;
            Dx2[i] = Dx[i] - (Dx1[i] - Dx[i]);
        }
        else if (i == (nx-2)-1)
        {
            Dx[i] = X[i+1] - X[i];
            Dx2[i] = (X[i+1] - X[i-1]) / 2.0;
            Dx1[i] = Dx[i] + (Dx[i] - Dx2[i]);
        }
        else
        {
            Dx[i] = X[i+1] - X[i];
            Dx1[i] = (X[i+2] - X[i]) / 2.0;
            Dx2[i] = (X[i+1] - X[i-1]) / 2.0;
        }
    }
    filer1(Dx, (nx-2), filePath + "coordinateDx_ny-" + to_string(ny));
    filer1(Dx1, (nx-2), filePath + "coordinateDx1_ny-" + to_string(ny));
    filer1(Dx2, (nx-2), filePath + "coordinateDx2_ny-" + to_string(ny));

    for (int j = 0; j < (ny-2); ++j)
    {
        if (j == 0)
        {
            Dy[j] = Y[j+1] - Y[j];
            Dy1[j] = (Y[j+2] - Y[j]) / 2.0;
            Dy2[j] = Dy[j] - (Dy1[j] - Dy[j]);
        }
        else if (j == (ny-2)-1)
        {
            Dy[j] = Y[j+1] - Y[j];
            Dy2[j] = (Y[j+1] - Y[j-1]) / 2.0;
            Dy1[j] = Dy[j] + (Dy[j] - Dy2[j]);
        }
        else
        {
            Dy[j] = Y[j+1] - Y[j];
            Dy1[j] = (Y[j+2] - Y[j]) / 2.0;
            Dy2[j] = (Y[j+1] - Y[j-1]) / 2.0;
        }
    }
    filer1(Dy, (ny-2), filePath + "coordinateDy_ny-" + to_string(ny));
    filer1(Dy1, (ny-2), filePath + "coordinateDy1_ny-" + to_string(ny));
    filer1(Dy2, (ny-2), filePath + "coordinateDy2_ny-" + to_string(ny));
}
