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
#include "printers.h"
#include "navierFVD.h"

using namespace std;

int t               = 0;
int pIter           = -1;
double uChange      = 1.0;
double pChange      = 1.0;
double uRun         = 0.0;  // used in momentum and momentumTimeStep
double vRun         = 0.0;  // used in momentum and momentumTimeStep

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
        // U[i][0] = -U[i][1] -2*uin;   // Dirichlet - attempt to model
                                        // nonzero v velocity vector at inlet
        V[i][0] = 0.0;                  // Dirichlet
    }
    // Left boundary
    for (int j = 0; j < ny; ++j)
    {
        U[0][j] = uin;                  // Dirichlet
        V[0][j] = -V[1][j];             // Dirichlet
        // V[0][j] = -V[1][j] + 2*vin;  // Dirichlet - attempt to model
                                        // nonzero v velocity vector at inlet
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
        // P[0][j] = pin;                  // Dirichlet
        P[0][j] = P[1][j];              // Neumann
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
    double ue, uw, un, us, ve, vw, vn, vs;
    for (int i = 1; i < nx-1; ++i)
    {
        for (int j = 1; j < ny-1; ++j)
        {
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

            un = (U[i][j+1] + U[i][j]) / 2.0;
            us = (U[i][j] + U[i][j-1]) / 2.0;
            ve = (V[i+1][j] + V[i][j]) / 2.0;
            vw = (V[i][j] + V[i-1][j]) / 2.0;

            FU[i][j] = - (ue*ue - uw*uw)/Dx1[i] - (un*vn - us*vs)/Dy[j]
                       + nu
                       * (
                          (
                           (U[i+1][j] - U[i][j])/Dx[i+1]
                           - (U[i][j] - U[i-1][j])/Dx[i]
                          )/Dx1[i]
                          +
                          (
                           (U[i][j+1] - U[i][j])/Dy1[j]
                           - (U[i][j] - U[i][j-1])/Dy2[j]
                          )/Dy[j]
                         );
            FV[i][j] = - (ue*ve - uw*vw)/Dx[i] - (vn*vn - vs*vs)/Dy1[j]
                       + nu
                       * (
                          (
                           (V[i+1][j] - V[i][j])/Dx1[i]
                           - (V[i][j] - V[i-1][j])/Dx2[i]
                          )/Dx[i]
                          +
                          (
                           (V[i][j+1] - V[i][j])/Dy[j+1]
                           - (V[i][j] - V[i][j-1])/Dy[j]
                          )/Dy1[j]
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
        pChange = 0.0;      // to find maximum change in p during timestep
        double change = 0.0;
        for (int i = 1; i < nx-1; ++i)
        {
            for (int j = 1; j < ny-1; ++j)
            {
                pTemp =
                    (
                     - P[i+1][j]*Dy[j]/Dx1[i] - P[i-1][j]*Dy[j]/Dx2[i]
                     - P[i][j+1]*Dx[i]/Dy1[j] - P[i][j-1]*Dx[i]/Dy2[j]
                     + ( (U[i][j] - U[i-1][j])*Dy[j]
                        + (V[i][j] - V[i][j-1])*Dx[i] )
                     *(rho / dt)
                    )
                    /
                    (- Dy[j]/Dx1[i] - Dy[j]/Dx2[i]
                     - Dx[i]/Dy1[j] - Dx[i]/Dy2[j]);

                change = abs(pTemp - P[i][j]);
                if (change > pChange)
                {
                    pChange = change;
                }
                // Applying the SOR method through following formula:
                P[i][j] = P[i][j] + omega * (pTemp - P[i][j]);
                // P[i][j] = (1.0 - omega)*P[i][j] + omega * pTemp; //same
            }
        }
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
            U[i][j] = U[i][j] - (dt/rho)*(P[i+1][j] - P[i][j])/Dx1[i];
            V[i][j] = V[i][j] - (dt/rho)*(P[i][j+1] - P[i][j])/Dy1[i];

            // uChange += abs(Uo[i][j] - U[i][j]);
            if (abs(Uo[i][j] - U[i][j]) > uChange)
            {
                uChange = abs(Uo[i][j] - U[i][j]);
            }
            // if (abs(Vo[i][j] - V[i][j]) > vChange)
            // {
            //     vChange = abs(Vo[i][j] - U[i][j]);
            // }
            Uo[i][j] = U[i][j];
            Vo[i][j] = V[i][j];
        }
    }
    // uChange = uChange / static_cast<double>(nx * ny);
    // vChange = vChange / static_cast<double>(nx * ny);
}

void navierFVD()
{
    initial();
    velBoundary();
    pressBoundary();
    filerAllSol(t);
    progressFileCreator();
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
        progressFiler();
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
