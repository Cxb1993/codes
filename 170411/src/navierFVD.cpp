/*
 * navierFVD.cpp
 *
 *  Created on: 2017-04-24
 *      Author: Syed Ahmad Raza
 */

#include <cmath>        // math functions
#include <iostream>     // functions for input and output to console
#include <sstream>      // functions for string conversion
// #include <algorithm>    // functions for ranges of elements (arrays)

#include "constants.h"
#include "gridder.h"
#include "filers.h"
#include "navierFVD.h"

using namespace std;

double U[nx][ny] = { {0.0} };
double V[nx][ny] = { {0.0} };
double Ur[nx][ny] = { {0.0} };
double Vr[nx][ny] = { {0.0} };
double FU[nx][ny] = { {0.0} };
double FV[nx][ny] = { {0.0} };
double FUo[nx][ny] = { {0.0} };
double FVo[nx][ny] = { {0.0} };
double P[nx][ny] = { {0.0} };
double Po[nx][ny] = { {0.0} };

void initial()
{
    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < ny; ++j)
        {
            U[i][j] = 0.0;
            V[i][j] = 0.0;
            P[i][j] = pin;
        }
    }
}

void boundary()
{
    // Applying the boundary conditions
    // Left boundary
    for (int j = 0; j < ny; ++j)
    {
        U[0][j] = uin;
        V[0][j] = vin;
        P[0][j] = P[1][j];
    }
    // Right boundary
    for (int j = 0; j < ny; ++j)
    {
        U[nx-1][j] = U[nx-2][j];
        V[nx-1][j] = V[nx-2][j];
        P[nx-1][j] = pin;
    }
    // Bottom boundary
    for (int i = 0; i < nx; ++i)
    {
        U[i][0] = 0.0;
        V[i][0] = 0.0;
        P[i][0] = P[i][1];
    }
    // Top boundary
    for (int i = 0; i < ny; ++i)
    {
        U[i][ny-1] = 0.0;
        V[i][ny-1] = 0.0;
        P[i][ny-1] = P[i][ny-2];
    }
}

void momentum()
{
    for (int i = 1; i < nx-1; ++i)
    {
        for (int j = 1; j < ny-1; ++j)
        {
            // Upwind scheme: flow is left to right, bottom to top
            double Ue = U[i][j];
            double Uw = U[i-1][j];
            double Vn = V[i][j];
            double Vs = V[i][j-1];

            FU[i][j] = -(Ue*Ue - Uw*Uw)/dx + (nu/dx)
                    *( (U[i+1][j] - U[i][j])/dx - (U[i][j] - U[i-1][j])/dx );
            FV[i][j] = -(Vn*Vn - Vs*Vs)/dy + (nu/dy)
                    *( (V[i][j+1] - V[i][j])/dy - (V[i][j] - V[i][j-1])/dy );

            Ur[i][j] = 3.0*FU[i][j]/2.0 - FUo[i][j]/2.0;
            Vr[i][j] = 3.0*FV[i][j]/2.0 - FVo[i][j]/2.0;

            FUo[i][j] = FU[i][j];
            FVo[i][j] = FV[i][j];
        }
    }
}

void momentumTimeStep()
{
    for (int i = 1; i < nx-1; ++i)
    {
        for (int j = 1; j < ny-1; ++j)
        {
            U[i][j] += dt*Ur[i][j];
            V[i][j] += dt*Vr[i][j];
        }
    }
}

void pressure()
{
    int m = 0;
    double rmsValue = 1.0;
    while (m <= maxPressIterations && rmsValue > residual)
    {
        double summation = 0.0;
        for (int i = 1; i < nx - 1; ++i)
        {
            for (int j = 1; j < ny - 1; ++j)
            {
                P[i][j] = (1.0 - omega)*P[i][j] + omega*
                (
                    (
                    - P[i+1][j]*dy/dx1 - P[i-1][j]*dy/dx2
                    - P[i][j+1]*dx/dy1 - P[i][j-1]*dx/dy2
                    + ( (U[i][j] - U[i-1][j])*dy + (V[i][j] - V[i][j-1])*dx ) /dt
                    )
                    /
                    (-dy/dx1 - dy/dx2 - dx/dy1 - dx/dy2)
                );
                summation += (P[i][j] - Po[i][j])*(P[i][j] - Po[i][j]);
                Po[i][j] = P[i][j];
            }
        }
        rmsValue = pow(summation / (nx*ny), 0.5);
        cout << m++ << ' ';
    }
}

void velocityCorrector()
{
    for (int i = 1; i < nx - 1; ++i)
    {
        for (int j = 1; j < ny - 1; ++j)
        {
            U[i][j] -= dt/rho * (P[i+1][j] - P[i-1][j])/2.0/dx;
            V[i][j] -= dt/rho * (P[i][j+1] - P[i][j-1])/2.0/dy;
        }
    }
}
