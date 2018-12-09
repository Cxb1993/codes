/*
 * gridder.cpp
 *
 *  Created on: 2017-09-27
 *      Author: Syed Ahmad Raza
 */

#include <cmath>        // math functions

#include "constants.h"
#include "filers.h"
#include "gridder.h"

using namespace std;

double X[nx+1]      = {0.0};        // x coordinates
double Y[ny+1]      = {0.0};        // y coordinates
double Dx[nx+1]     = {0.0};
double Dxs[nx+1]    = {0.0};
double Dy[ny+1]     = {0.0};
double Dys[ny+1]    = {0.0};

void gridder()
{
    for (int i = 1; i < nx; ++i)
    {
        // X[i] = (0.5 * lx) * (1 + sin( pi *
        //     ((static_cast<double>(i) / static_cast<double>(nx-1)) - 0.5)));
        // X[i] = lx *
        // (sin( pi * (static_cast<double>(i) / static_cast<double>(ny-1)) * 0.5));
        X[i] = static_cast<double>(i-1) * L / static_cast<double>(nx-2);
    }
    X[0]    = X[1];
    X[nx]   = X[nx-1];
    for (int j = 1; j < ny; ++j)
    {
        // Y[j] = 0.5 * D *
        //         (
        //          tan(
        //              pi * (
        //                    ((static_cast<double>(j-1)/static_cast<double>(ny-2))
        //                    / 2.0) - 0.25
        //                   )
        //             ) + 1
        //         );
        Y[j] = static_cast<double>(j-1) * D / static_cast<double>(ny-2);
    }
    Y[0]    = Y[1];
    Y[ny]   = Y[ny-1];
    filer1(X, (nx+1), filePath + "ny=" + to_string(ny) + fileUniqueName
            + "_coordinateX");
    filer1(Y, (ny+1), filePath + "ny=" + to_string(ny) + fileUniqueName
            + "_coordinateY");

    for (int i = 0; i < (nx+1); ++i)
    {
        if ((i > 1) && (i < (nx-1)))
        {
            Dx[i]   = ( X[i+1] - X[i]   );
            Dxs[i]  = ( X[i+1] - X[i-1] ) / 2.0;
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

    filer1(Dx, (nx+1), filePath + "ny=" + to_string(ny) + fileUniqueName
            + "_coordinateDx");
    filer1(Dxs, (nx+1), filePath + "ny=" + to_string(ny) + fileUniqueName
            + "_coordinateDxs");

    for (int j = 0; j < (ny+1); ++j)
    {
        if ((j > 1) && (j < (ny-1)))
        {
            Dy[j]   = ( Y[j+1] - Y[j]   );
            Dys[j]  = ( Y[j+1] - Y[j-1] ) / 2.0;
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

    filer1(Dy, (ny+1), filePath + "ny=" + to_string(ny) + fileUniqueName
            + "_coordinateDy");
    filer1(Dys, (ny+1), filePath + "ny=" + to_string(ny) + fileUniqueName
            + "_coordinateDys");
}
