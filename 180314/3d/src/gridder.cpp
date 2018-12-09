/*
 * gridder.cpp
 *
 *  Created on: 2018-01-26
 *      Author: Syed Ahmad Raza
 */

#include <cmath>        // math functions

#include "constants.h"
#include "filers.h"
#include "gridder.h"

using namespace std;

double X[nx+1]      = {0.0};        // x coordinates
double Y[ny+1]      = {0.0};        // y coordinates
double Z[nz+1]      = {0.0};        // y coordinates
double Dx[nx+1]     = {0.0};
double Dxs[nx+1]    = {0.0};
double Dy[ny+1]     = {0.0};
double Dys[ny+1]    = {0.0};
double Dz[nz+1]     = {0.0};
double Dzs[nz+1]    = {0.0};

void gridder()
{
    // Defining x coordinates
    for (int i = 1; i < nx; ++i)
    {
        X[i] = static_cast<double>(i-1) * L / static_cast<double>(nx-2);
    }
    X[0]    = X[1];
    X[nx]   = X[nx-1];

    // Defining y coordinates
    for (int j = 1; j < ny; ++j)
    {
        Y[j] = 0.5 * H *
                (
                 tan(
                     pi * (
                           ((static_cast<double>(j-1)/static_cast<double>(ny-2))
                           / 2.0) - 0.25
                          )
                    ) + 1
                );
        // Y[j] = static_cast<double>(j-1) * H / static_cast<double>(ny-2);
    }
    Y[0]    = Y[1];
    Y[ny]   = Y[ny-1];

    // Defining z coordinates
    for (int k = 1; k < nz; ++k)
    {
        Z[k] = 0.5 * B *
                (
                 tan(
                     pi * (
                           ((static_cast<double>(k-1)/static_cast<double>(nz-2))
                           / 2.0) - 0.25
                          )
                    ) + 1
                );
        // Z[k] = static_cast<double>(k-1) * B / static_cast<double>(nz-2);
    }
    Z[0]    = Z[1];
    Z[nz]   = Z[nz-1];

    filer1(X, (nx+1), filePath + "coordinateX_ny-" + to_string(ny));
    filer1(Y, (ny+1), filePath + "coordinateY_ny-" + to_string(ny));
    filer1(Z, (nz+1), filePath + "coordinateZ_ny-" + to_string(ny));

    // Defining x-axis intervals
    for (int i = 0; i < (nx+1); ++i)
    {
        if ((i > 1) && (i < (nx-1)))
        {
            Dx[i]   = X[i+1] - X[i];
            Dxs[i]  = (X[i+1] - X[i-1]) / 2.0;
        }
        else if (i >= (nx-1))
        {
            Dx[i]   = Dx[nx-2];
            Dxs[i]  = Dxs[nx-2];
        }
    }
    Dx[1]   = X[1+1] - X[1];
    Dxs[1]  = Dxs[2];
    Dx[0]   = Dx[1];
    Dxs[0]  = Dxs[2];

    filer1(Dx, (nx+1), filePath + "coordinateDx_ny-" + to_string(ny));
    filer1(Dxs, (nx+1), filePath + "coordinateDxs_ny-" + to_string(ny));

    // Defining y-axis intervals
    for (int j = 0; j < (ny+1); ++j)
    {
        if ((j > 1) && (j < (ny-1)))
        {
            Dy[j]   = Y[j+1] - Y[j];
            Dys[j]  = (Y[j+1] - Y[j-1]) / 2.0;
        }
        else if (j >= (ny-1))
        {
            Dy[j]   = Dy[ny-2];
            Dys[j]  = Dys[ny-2];
        }
    }
    Dy[1]   = Y[1+1] - Y[1];
    Dys[1]  = Dys[2];
    Dy[0]   = Dy[1];
    Dys[0]  = Dys[2];

    filer1(Dy, (ny+1), filePath + "coordinateDy_ny-" + to_string(ny));
    filer1(Dys, (ny+1), filePath + "coordinateDys_ny-" + to_string(ny));

    // Defining z-axis intervals
    for (int k = 0; k < (nz+1); ++k)
    {
        if ((k > 1) && (k < (nz-1)))
        {
            Dz[k]   = Z[k+1] - Z[k];
            Dzs[k]  = (Z[k+1] - Z[k-1]) / 2.0;
        }
        else if (k >= (nz-1))
        {
            Dz[k]   = Dz[nz-2];
            Dzs[k]  = Dzs[nz-2];
        }
    }
    Dz[1]   = Z[1+1] - Z[1];
    Dzs[1]  = Dzs[2];
    Dz[0]   = Dz[1];
    Dzs[0]  = Dzs[2];

    filer1(Dz, (nz+1), filePath + "coordinateDz_ny-" + to_string(ny));
    filer1(Dzs, (nz+1), filePath + "coordinateDzs_ny-" + to_string(ny));
}
