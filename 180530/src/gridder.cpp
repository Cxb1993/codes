/*
 * gridder.cpp
 *
 *     Project: Finite Volume Navier-Stokes Solver
 *      Author: Syed Ahmad Raza
 */

#include <cmath>        // math functions
// #include <iostream>

#include "constants.h"
#include "filers.h"
#include "gridder.h"

// using namespace std;
using namespace Numeric_lib;

// Initial grid coordinates for evaluating grid lengths
Matrix<double,1> X(nx+1);
Matrix<double,1> Y(ny+1);
Matrix<double,1> Z(nz+1);

// Actual grid cooridinates (with adjusted index)
Matrix<double,1> Xa(nx-3);
Matrix<double,1> Ya(ny-3);
Matrix<double,1> Za(nz-3);

// Midpoints of grid coordinates
Matrix<double,1> Xs(nx-4);
Matrix<double,1> Ys(ny-4);
Matrix<double,1> Zs(ny-4);

// Grid lengths
Matrix<double,1> Dx(nx);
Matrix<double,1> Dxs(nx);
Matrix<double,1> Dy(ny);
Matrix<double,1> Dys(ny);
Matrix<double,1> Dz(nz);
Matrix<double,1> Dzs(nz);

void gridder()
{
    // Define the grid coordinates

    for (int i = 2; i <= nx-2; ++i)
    {
        // // For finer grid intervals in increasing x direction:
        // X(i) = L * sin( 0.5 * pi
        //                * (static_cast<double>(i-2)
        //                   / static_cast<double>(nx-4))
        //               );

        // // For finer grid intervals in the start and end of x axis:
        // X(i) = (0.5 * L)
        //        * ( 1 + sin(pi*(- 0.5
        //                        + (static_cast<double>(i-2) /
        //                           static_cast<double>(nx-4) )
        //                       )
        //                   )
        //          );

        // For equal grid intervals:
        X(i) = L * (static_cast<double>(i-2) / static_cast<double>(nx-4));
        // In other words, L * (i/nx), where (i/nx) is simply the x position
    }

    for (int j = 2; j <= ny-2; ++j)
    {
        // // For finer grid intervals in increasing y direction:
        // Y(j) = B * sin( 0.5 * pi
        //                * (static_cast<double>(j-2)
        //                   / static_cast<double>(ny-4))
        //               );

        // // For finer grid intervals in the start and end of y axis:
        // Y(j) = (0.5 * B)
        //        * ( 1 + sin(pi*(- 0.5
        //                        + (static_cast<double>(j-2) /
        //                           static_cast<double>(ny-4) )
        //                       )
        //                   )
        //          );

        // For equal grid intervals:
        Y(j) = B * (static_cast<double>(j-2) / static_cast<double>(ny-4));
        // In other words, B * (j/ny), where (j/ny) is simply the y position
    }

    for (int k = 2; k <= nz-2; ++k)
    {
        // // For finer grid intervals in increasing y direction:
        // Z(k) = H * sin( 0.5 * pi
        //                * (static_cast<double>(k-2)
        //                   / static_cast<double>(nz-4))
        //               );

        // // For finer grid intervals in the start and end of z axis:
        // Z(k) = (0.5 * H)
        //        * ( 1 + sin(pi*(- 0.5
        //                        + (static_cast<double>(k-2) /
        //                           static_cast<double>(nz-4) )
        //                       )
        //                   )
        //          );

        // For equal grid intervals:
        Z(k) = H * (static_cast<double>(k-2) / static_cast<double>(nz-4));
        // In other words, H * (k/nz), where (i/nz) is simply the y position
    }

    // Define each of the directional grid lengths

    for (int i = 2; i <= nx-4; ++i)
    {
        Dx (i)  = ( X(i+1) - X(i) )      ;
        Dxs(i)  = ( X(i+2) - X(i) ) / 2.0;
    }

    for (int j = 2; j <= ny-4; ++j)
    {
        Dy (j)  = ( Y(j+1) - Y(j) )      ;
        Dys(j)  = ( Y(j+2) - Y(j) ) / 2.0;
    }

    for (int k = 2; k <= nz-4; ++k)
    {
        Dz (k)  = ( Z(k+1) - Z(k) )      ;
        Dzs(k)  = ( Z(k+2) - Z(k) ) / 2.0;
    }

    // Ghost boundary grid lengths

    Dx (1)    = Dx(2);
    Dx (0)    = Dx(2);
    Dx (nx-3) = X(nx-2) - X(nx-3);
    Dx (nx-2) = Dx(nx-3);
    Dx (nx-1) = Dx(nx-3);

    Dxs(1)    = Dxs(2);
    Dxs(0)    = Dxs(2);
    Dxs(nx-3) = Dxs(nx-4);
    Dxs(nx-2) = Dxs(nx-4);
    Dxs(nx-1) = Dxs(nx-4);

    Dy (1)    = Dy(2);
    Dy (0)    = Dy(2);
    Dy (ny-3) = Y(ny-2) - Y(ny-3);
    Dy (ny-2) = Dy(ny-3);
    Dy (ny-1) = Dy(ny-3);

    Dys(1)    = Dys(2);
    Dys(0)    = Dys(2);
    Dys(ny-3) = Dys(ny-4);
    Dys(ny-2) = Dys(ny-4);
    Dys(ny-1) = Dys(ny-4);

    Dz (1)    = Dz(2);
    Dz (0)    = Dz(2);
    Dz (nz-3) = Z(nz-2) - Z(nz-3);
    Dz (nz-2) = Dz(nz-3);
    Dz (nz-1) = Dz(nz-3);

    Dzs(1)    = Dzs(2);
    Dzs(0)    = Dzs(2);
    Dzs(nz-3) = Dzs(nz-4);
    Dzs(nz-2) = Dzs(nz-4);
    Dzs(nz-1) = Dzs(nz-4);

    // Modifying the index of X, Y and Z arrays to represent the actual grid

    for (int i = 0; i <= nx-4; ++i)
    {
        Xa(i) = X(i+2);
    }

    for (int j = 0; j <= ny-4; ++j)
    {
        Ya(j) = Y(j+2);
    }

    for (int k = 0; k <= nz-4; ++k)
    {
        Za(k) = Z(k+2);
    }

    // Defining the midpoint values of the grids

    for (int i = 0; i <= nx-5; ++i)
    {
        Xs(i) = ( Xa(i+1) + Xa(i) ) / 2.0;
    }

    for (int j = 0; j <= ny-5; ++j)
    {
        Ys(j) = ( Ya(j+1) + Ya(j) ) / 2.0;
    }

    for (int k = 0; k <= nz-5; ++k)
    {
        Zs(k) = ( Za(k+1) + Za(k) ) / 2.0;
    }
}
