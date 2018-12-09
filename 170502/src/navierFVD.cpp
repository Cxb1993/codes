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
double data[nx][ny] = { {0.0} };

void initial()
{
    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < ny; ++j)
        {
            U[i][j] = 0.0;
            V[i][j] = 0.0;
            P[i][j] = 0.0;
        }
    }
}

void boundary()
{
    // Applying the boundary conditions
    // Top boundary
    for (int i = 0; i < ny; ++i)
    {
        U[i][ny-1] = -U[i][ny-2];       // Dirichlet
        V[i][ny-1] = 0.0;               // Dirichlet
        // P[i][ny-1] = P[i][ny-2];     // Neumann
    }
    // Bottom boundary
    for (int i = 0; i < nx; ++i)
    {
        U[i][0] = -U[i][1];             // Dirichlet
        V[i][0] = 0.0;                  // Dirichlet
        // P[i][0] = P[i][1];           // Neumann
    }
    // Left boundary
    for (int j = 0; j < ny; ++j)
    {
        U[0][j] = uin;                  // Dirichlet
        V[0][j] = -V[0][j] + 2*vin;     // Dirichlet
        // P[0][j] = pin;               // Dirichlet
    }
    // Right boundary
    for (int j = 0; j < ny; ++j)
    {
        U[nx-1][j] = U[nx-2][j];        // Neumann
        V[nx-1][j] = V[nx-2][j];        // Neumann
        // P[nx-1][j] = -P[nx-2][j];    // Dirichlet
    }
}

void momentum()
{
    for (int i = 1; i < nx-1; ++i)
    {
        for (int j = 1; j < ny-1; ++j)
        {
            // Upwind scheme: flow is left to right, bottom to top
            double ue, uw, un, us, ve, vw, vn, vs;
            if (U[i][j] >= 0.0)
            {
                ue = U[i][j];
                uw = U[i-1][j];
            }
            else
            {
                ue = U[i+1][j];
                uw = U[i][j];
            }
            if (V[i][j] >= 0.0)
            {
                vn = V[i][j];
                vs = V[i][j-1];
            }
            else
            {
                vn = V[i][j+1];
                vs = V[i][j];
            }
            un = U[i][j] + (U[i][j+1] - U[i][j])/dy;
            us = U[i][j] - (U[i][j] - U[i][j-1])/dy;
            ve = V[i][j] + (V[i+1][j] - V[i][j])/dx;
            vw = V[i][j] - (V[i][j] - V[i-1][j])/dx;

            FU[i][j] = -(ue*ue - uw*uw)/dx -(un*vn - us*vs)/dy
                    + nu
                    *(
                      (
                       (U[i+1][j] - U[i][j])/dx - (U[i][j] - U[i-1][j])/dx
                      )/dx
                      +
                      (
                       (U[i][j+1] - U[i][j])/dy - (U[i][j] - U[i][j-1])/dy
                      )/dy
                     );
            FV[i][j] = -(ue*ve - uw*vw)/dx -(vn*vn - vs*vs)/dy
                     + nu
                     *(
                       (
                        (V[i+1][j] - V[i][j])/dx - (V[i][j] - V[i-1][j])/dx
                       )/dx
                       +
                       (
                        (V[i][j+1] - V[i][j])/dy - (V[i][j] - V[i][j-1])/dy
                       )/dy
                      );
            if (t==0)
            {
                Ur[i][j] = FU[i][j];
                Vr[i][j] = FV[i][j];
            }
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

void temp()
{
    for (int i = 1; i < nx-1; ++i)
    {
        for (int j = 1; j < ny-1; ++j)
        {
            data[i][j] =    (
                            (U[i][j] - U[i-1][j])*dy +
                            (V[i][j] - V[i][j-1])*dx
                            )/dt;
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
