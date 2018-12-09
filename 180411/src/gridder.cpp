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
using namespace Numeric_lib;

Matrix<double,1> X(nx+1);
Matrix<double,1> Y(ny+1);
Matrix<double,1> Xa(nx-3);
Matrix<double,1> Ya(ny-3);
Matrix<double,1> Dx(nx);
Matrix<double,1> Dxs(nx);
Matrix<double,1> Dy(ny);
Matrix<double,1> Dys(ny);

void gridder()
{
    // Define the grid coordinates

    for (int i = 2; i <= nx-2; ++i)
    {
        // X[i] = (0.5 * lx)
        //        * (1 + sin( pi * ((static_cast<double>(i)
        //        / static_cast<double>(nx-1)) - 0.5)));
        // X[i] = lx *
        //        (sin( pi * (static_cast<double>(i)
        //        / static_cast<double>(ny-1)) * 0.5));
        X(i) = static_cast<double>(i-2) * L / static_cast<double>(nx-4);
    }

    for (int j = 2; j <= ny-2; ++j)
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
        Y(j) = static_cast<double>(j-2) * D / static_cast<double>(ny-4);
    }

    // Define the x and y direction grid lengths

    for (int i = 3; i <= nx-3; ++i)
    {
        Dx (i)  = ( X(i+1) - X(i)   );
        Dxs(i)  = ( X(i+1) - X(i-1) ) / 2.0;
    }

    for (int j = 3; j <= ny-3; ++j)
    {
        Dy (j)  = ( Y(j+1) - Y(j)   );
        Dys(j)  = ( Y(j+1) - Y(j-1) ) / 2.0;
    }

    // Ghost boundary grid lengths

    Dx (2)    = X(3) - X(2);
    Dx (1)    = Dx(2);
    Dx (0)    = Dx(2);
    Dx (nx-2) = Dx(nx-3);
    Dx (nx-1) = Dx(nx-3);

    Dxs(2)    = Dxs(3);
    Dxs(1)    = Dxs(2);
    Dxs(0)    = Dxs(2);
    Dxs(nx-2) = Dxs(nx-3);
    Dxs(nx-1) = Dxs(nx-3);

    Dy (2)    = Y(3) - Y(2);
    Dy (1)    = Dy(2);
    Dy (0)    = Dy(2);
    Dy (ny-2) = Dy(ny-3);
    Dy (ny-1) = Dy(ny-3);

    Dys(2)    = Dys(3);
    Dys(1)    = Dys(2);
    Dys(0)    = Dys(2);
    Dys(ny-2) = Dys(ny-3);
    Dys(ny-1) = Dys(ny-3);

    // Modifying the index of X and Y arrays to represent the actual grid

    for (int i = 0; i <= nx-4; ++i)
    {
        Xa(i) = X(i+2);
    }
    for (int j = 0; j <= ny-4; ++j)
    {
        Ya(j) = Y(j+2);
    }

    // Filing the coordinate files

    filer1(Xa, (nx-3), filePath + "nxy=" + to_string(nxy) + fileUniqueName
            + "_coordinateX");
    filer1(Ya, (ny-3), filePath + "nxy=" + to_string(nxy) + fileUniqueName
            + "_coordinateY");

    filer1(Dx, (nx), filePath + "nxy=" + to_string(nxy) + fileUniqueName
            + "_coordinateDx");
    filer1(Dxs, (nx), filePath + "nxy=" + to_string(nxy) + fileUniqueName
            + "_coordinateDxs");

    filer1(Dy, (ny), filePath + "nxy=" + to_string(nxy) + fileUniqueName
            + "_coordinateDy");
    filer1(Dys, (ny), filePath + "nxy=" + to_string(nxy) + fileUniqueName
            + "_coordinateDys");
}
