/*
 * navierFVD.cpp
 *
 *  Created on: 2017-09-27
 *      Author: Syed Ahmad Raza
 */

#include <cmath>        // math functions
// #include <algorithm>    // functions for ranges of elements (arrays)

#include "constants.h"
#include "gridder.h"
#include "filers.h"
#include "navierFVD.h"

using namespace std;

int t               = 0;
double uChange      = 1.0;
double pChange      = 1.0;
int pIter           = -1;

double U[nx][ny]    = { {0.0} };
double V[nx][ny]    = { {0.0} };
double Ur[nx][ny]   = { {0.0} };
double Vr[nx][ny]   = { {0.0} };
double Uo[nx][ny]   = { {0.0} };
double Vo[nx][ny]   = { {0.0} };
double FU[nx][ny]   = { {0.0} };
double FV[nx][ny]   = { {0.0} };
double FUo[nx][ny]  = { {0.0} };
double FVo[nx][ny]  = { {0.0} };
double P[nx][ny]    = { {0.0} };

void initial()
{
    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < ny; ++j)
        {
            U[i][j] = ubegin;
            V[i][j] = vbegin;
            P[i][j] = pbegin;
        }
    }
}

void velBoundary()
{
    // Applying the velocity boundary conditions
    // Top boundary
    for (int i = 0; i < nx; ++i)
    {
        U[i][ny-1] = -U[i][ny-2];       // Dirichlet
        // U[i][ny-1] = -U[i][ny-2] + 2*uin;       // Dirichlet
        V[i][ny-1] = 0.0;               // Dirichlet
    }
    // Bottom boundary
    for (int i = 0; i < nx; ++i)
    {
        U[i][0] = -U[i][1];             // Dirichlet
        // U[i][0] = -U[i][1] -2*uin;      // Dirichlet
        V[i][0] = 0.0;                  // Dirichlet
    }
    // Left boundary
    for (int j = 0; j < ny; ++j)
    {
        U[0][j] = uin;                  // Dirichlet
        V[0][j] = -V[1][j];             // Dirichlet
        // V[0][j] = -V[1][j] + 2*vin;     // Dirichlet - attempt to model nonzero v velocity vector at inlet
    }
    // Right boundary
    for (int j = 0; j < ny; ++j)
    {
        U[nx-1][j] = U[nx-2][j];        // Neumann
        V[nx-1][j] = V[nx-2][j];        // Neumann
    }
}

void pressBoundary()
{
    // Applying the pressure boundary conditions

    // Left boundary
    for (int j = 0; j < ny; ++j)
    {
        P[0][j] = pin;                  // Dirichlet
        // P[0][j] = P[1][j];              // Neumann
    }
    // Right boundary
    for (int j = 0; j < ny; ++j)
    {
        P[nx-1][j] = pout;              // Dirichlet
        // P[nx-1][j] = P[nx-2][j];        // Neumann
    }
    // Top boundary
    for (int i = 0; i < nx; ++i)
    {
        P[i][ny-1] = P[i][ny-2];        // Neumann
    }
    // Bottom boundary
    for (int i = 0; i < nx; ++i)
    {
        P[i][0] = P[i][1];              // Neumann
    }
}

void momentum(string velocityScheme)
{
    for (int i = 1; i < nx-1; ++i)
    {
        for (int j = 1; j < ny-1; ++j)
        {
            double ue, uw, un, us, ve, vw, vn, vs;

            // Central scheme
            if (velocityScheme == "ct")
            {
                ue = (U[i][j] + U[i+1][j]) / 2.0;
                uw = (U[i-1][j] + U[i][j]) / 2.0;
                vn = (V[i][j] + V[i][j+1]) / 2.0;
                vs = (V[i][j-1] + V[i][j]) / 2.0;
            }
            // Upwind scheme
            else if (velocityScheme == "up")
            {
                if (U[i][j] > 0.0)
                {
                    ue = U[i][j];
                    uw = U[i-1][j];
                }
                else
                {
                    ue = U[i+1][j];
                    uw = U[i][j];
                }
                if (V[i][j] > 0.0)
                {
                    vn = V[i][j];
                    vs = V[i][j-1];
                }
                else
                {
                    vn = V[i][j+1];
                    vs = V[i][j];
                }
            }
            // QUICK scheme
            else if (velocityScheme == "qk")
            {
                if (U[i][j] > 0.0)  // positive
                {
                    if (i == 1 || i == nx-2)
                    {
                        ue = U[i][j];
                        uw = U[i-1][j];
                    }
                    else
                    {
                        ue = U[i][j]*6/8 + U[i+1][j]*3/8 - U[i-1][j]/8;
                        uw = U[i-1][j]*6/8 + U[i][j]*3/8 - U[i-2][j]/8;
                    }
                }
                else                // negative
                {
                    if (i == 1 || i == nx-2)
                    {
                        ue = U[i+1][j];
                        uw = U[i][j];
                    }
                    else
                    {
                        ue = U[i+1][j]*6/8 + U[i][j]*3/8 - U[i+2][j]/8;
                        uw = U[i][j]*6/8 + U[i-1][j]*3/8 - U[i+1][j]/8;
                    }
                }
                if (V[i][j] > 0.0)  // positive
                {
                    if (j == 1 || j == ny-2)
                    {
                        vn = V[i][j];
                        vs = V[i][j-1];
                    }
                    else
                    {
                        vn = V[i][j]*6/8 + V[i][j+1]*3/8 - V[i][j-1]/8;
                        vs = V[i][j-1]*6/8 + V[i][j]*3/8 - V[i][j-2]/8;
                    }
                }
                else                // negative
                {
                    if (j == 1 || j == ny-2)
                    {
                        vn = V[i][j+1];
                        vs = V[i][j];
                    }
                    else
                    {
                        vn = V[i][j+1]*6/8 + V[i][j]*3/8 - V[i][j+2]/8;
                        vs = V[i][j]*6/8 + V[i][j-1]*3/8 - V[i][j+1]/8;
                    }
                }
            }
            else
            {
                printf("\nUndefined velocity scheme selected\n");
                // initialize ue, uw, vn, vs to avoid warning
                // "=Wmaybe-uninitialized"
                ue = (U[i][j] + U[i+1][j]) / 2.0;
                uw = (U[i-1][j] + U[i][j]) / 2.0;
                vn = (V[i][j] + V[i][j+1]) / 2.0;
                vs = (V[i][j-1] + V[i][j]) / 2.0;
            }

            un = (U[i][j+1] + U[i][j]) / 2;
            us = (U[i][j] + U[i][j-1]) / 2;
            ve = (V[i+1][j] + V[i][j]) / 2;
            vw = (V[i][j] + V[i-1][j]) / 2;

            FU[i][j] = - (ue*ue - uw*uw)/dx - (un*vn - us*vs)/dy
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
            FV[i][j] = - (ue*ve - uw*vw)/dx - (vn*vn - vs*vs)/dy
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
            if (t == 1)     // Euler scheme for the first time-step
            {
                Ur[i][j] = FU[i][j];
                Vr[i][j] = FV[i][j];
            }
            else
            {
                Ur[i][j] = 3.0*FU[i][j]/2.0 - FUo[i][j]/2.0;
                Vr[i][j] = 3.0*FV[i][j]/2.0 - FVo[i][j]/2.0;
            }

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
            U[i][j] += dt * Ur[i][j];
            V[i][j] += dt * Vr[i][j];
        }
    }
}

void pressure()
{
    int m = 0;
    double pTemp = 0.0;
    pChange = 1.0;
    while (m <= maxPressIters && pChange > pResidual)
    {
        pChange = 0.0;
        for (int i = 1; i < nx-1; ++i)
        {
            for (int j = 1; j < ny-1; ++j)
            {
                pTemp =
                    (
                    - P[i+1][j]*dy/dx1 - P[i-1][j]*dy/dx2
                    - P[i][j+1]*dx/dy1 - P[i][j-1]*dx/dy2
                    + ( (U[i][j] - U[i-1][j])*dy + (V[i][j] - V[i][j-1])*dx ) *(rho / dt)
                    )
                    /
                    (- dy/dx1 - dy/dx2 - dx/dy1 - dx/dy2);

                pChange += abs(pTemp - P[i][j]);

                // Applying the SOR method through following formula:
                P[i][j] = P[i][j] + omega * (pTemp - P[i][j]);
                // P[i][j] = (1.0 - omega)*P[i][j] + omega * pTemp; //same
            }
        }
        pChange = pChange / static_cast<double>(nx * ny);
        pressBoundary();
        ++m;
    }
    pIter = --m;
}

void velUpdater()
{
    for (int i = 1; i < nx-1; ++i)
    {
        for (int j = 1; j < ny-1; ++j)
        {
            U[i][j] = U[i][j] - (dt / rho) * (P[i+1][j] - P[i][j])/dx;
            V[i][j] = V[i][j] - (dt / rho) * (P[i][j+1] - P[i][j])/dy;

            uChange += abs(Uo[i][j] - U[i][j]);

            Uo[i][j] = U[i][j];
            Vo[i][j] = V[i][j];
        }
    }
    uChange = uChange / static_cast<double>(nx * ny);
    // vChange = vChange / static_cast<double>(nx * ny);
}

void navierFVD()
{
    initial();
    velBoundary();
    pressBoundary();
    filerAllSol(t);
    // End of timestep
    t = 1;
    while   (t <= 10 || (t <= maxTimesteps && uChange > uResidual))
    {
        uChange = 0.0;

        momentum(velScheme);
        momentumTimeStep();
        velBoundary();
        pressure();
        velUpdater();
        velBoundary();
        progressPrinter();
        // if ((t <= 2) || (t % 1000 == 0))
        // {
        //     filerAllSol(t);
        // }
        // End of timestep
        ++t;
    }
    --t;
    filerAllSol();
}