/*
 * navierFVD.cpp
 *
 *  Created on: 2018-01-26
 *      Author: Syed Ahmad Raza
 */

#include <cmath>        // math functions
// #include <algorithm>    // functions for ranges of elements (arrays)
#include <iostream>

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
double wRun         = 0.0;  // used in momentum and momentumTimeStep

double U[nx][ny][nz]    = { {0.0} };
double V[nx][ny][nz]    = { {0.0} };
double W[nx][ny][nz]    = { {0.0} };
double Ur[nx][ny][nz]   = { {0.0} };
double Vr[nx][ny][nz]   = { {0.0} };
double Wr[nx][ny][nz]   = { {0.0} };
double Uo[nx][ny][nz]   = { {0.0} };
double Vo[nx][ny][nz]   = { {0.0} };
double Wo[nx][ny][nz]   = { {0.0} };
double FU[nx][ny][nz]   = { {0.0} };
double FV[nx][ny][nz]   = { {0.0} };
double FW[nx][ny][nz]   = { {0.0} };
double FUo[nx][ny][nz]  = { {0.0} };
double FVo[nx][ny][nz]  = { {0.0} };
double FWo[nx][ny][nz]  = { {0.0} };
double P[nx][ny][nz]    = { {0.0} };

void initial()
{
    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < ny; ++j)
        {
            for (int k = 0; k < nz; ++k)
            {
                U[i][j][k] = ubegin;
                V[i][j][k] = vbegin;
                W[i][j][k] = wbegin;
                P[i][j][k] = pbegin;
            }
        }
    }
}

void velBoundary()
{
    // Applying the velocity boundary conditions
    // South boundary
    for (int i = 0; i < nx; ++i)
    {
        for (int k = 0; k < nz; ++k)
        {
            U[i][0][k] = -U[i][1][k];       // Dirichlet
            V[i][0][k] = 0.0;               // Dirichlet
            W[i][0][k] = -W[i][1][k];       // Dirichlet
        }
    }
    // North boundary
    for (int i = 0; i < nx; ++i)
    {
        for (int k = 0; k < nz; ++k)
        {
            U[i][ny-1][k] = -U[i][ny-2][k]; // Dirichlet
            V[i][ny-1][k] = 0.0;            // Dirichlet
            W[i][ny-1][k] = -W[i][ny-2][k]; // Dirichlet
        }
    }
    // West boundary
    for (int j = 0; j < ny; ++j)
    {
        for (int k = 0; k < nz; ++k)
        {
            U[0][j][k] = uin;               // Dirichlet
            V[0][j][k] = -V[1][j][k];       // Dirichlet
            W[0][j][k] = -W[1][j][k];       // Dirichlet
        }
    }
    // East boundary
    for (int j = 0; j < ny; ++j)
    {
        for (int k = 0; k < nz; ++k)
        {
            U[nx-1][j][k] = U[nx-2][j][k];  // Neumann
            V[nx-1][j][k] = V[nx-2][j][k];  // Neumann
            W[nx-1][j][k] = V[nx-2][j][k];  // Neumann
        }
    }
    // Bottom boundary
    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < ny; ++j)
        {
            U[i][0][k] = -U[i][1][k];       // Dirichlet
            V[i][0][k] = -V[i][1][k];       // Dirichlet
            W[i][0][k] = 0.0;               // Dirichlet
        }
    }
    // Top boundary
    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < ny; ++j)
        {
            U[i][ny-1][k] = -U[i][ny-2][k]; // Dirichlet
            V[i][ny-1][k] = -W[i][ny-2][k]; // Dirichlet
            W[i][ny-1][k] =  0.0;           // Dirichlet
        }
    }
}

void pressBoundary()
{
    // Applying the pressure boundary conditions

    // West boundary
    for (int j = 0; j < ny; ++j)
    {
        for (int k = 0; k < nz; ++k)
        {
            // P[0][j][k] = pin;               // Dirichlet
            P[0][j][k] = P[1][j][k];        // Neumann
        }
    }
    // East boundary
    for (int j = 0; j < ny; ++j)
    {
        for (int k = 0; k < nz; ++k)
        {
            P[nx-1][j][k] = pout;           // Dirichlet
            // P[nx-1][j][k] = P[nx-2][j][k];  // Neumann
        }
    }
    // South boundary
    for (int i = 0; i < nx; ++i)
    {
        for (int k = 0; k < nz; ++k)
        {
            P[i][0][k] = P[i][1][k];        // Neumann
        }
    }
    // North boundary
    for (int i = 0; i < nx; ++i)
    {
        for (int k = 0; k < nz; ++k)
        {
            P[i][ny-1][k] = P[i][ny-2][k];  // Neumann
        }
    }
    // Bottom boundary
    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < ny; ++j)
        {
            P[i][j][0] = P[i][j][1];        // Neumann
        }
    }
    // Top boundary
    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < ny; ++j)
        {
            P[i][j][nz-1] = P[i][j][nz-2];  // Neumann
        }
    }
}

void momentum(string velocityScheme)
{
    double ue, uw, un, us, uev, uwv, ve, vw, vn, vs, vnu, vsu;
    for (int i = 1; i < nx-1; ++i)
    {
        for (int j = 1; j < ny-1; ++j)
        {
            // Upwind scheme
            else if (velocityScheme == "up")
            {
                if (U[i][j] > 0.0)
                {
                    ue  = U[i][j];
                    uw  = U[i-1][j];
                    uev = U[i][j];
                    uwv = U[i-1][j];
                }
                else
                {
                    ue  = U[i+1][j];
                    uw  = U[i][j];
                    uev = U[i][j+1];
                    uwv = U[i-1][j+1];
                }
                if (V[i][j] > 0.0)
                {
                    vn  = V[i][j];
                    vs  = V[i][j-1];
                    vnu = V[i][j];
                    vsu = V[i][j-1];
                }
                else
                {
                    vn  = V[i][j+1];
                    vs  = V[i][j];
                    vnu = V[i+1][j];
                    vsu = V[i+1][j-1];
                }
            }
            else
            {
                printf("\nUndefined velocity scheme selected\n");
                // initialize ue, uw, uev, uwv, vn, vs, vnu, vsu to avoid
                // warning "=Wmaybe-uninitialized"
                ue  = 0.0;
                uw  = 0.0;
                uev = 0.0;
                uwv = 0.0;
                vn  = 0.0;
                vs  = 0.0;
                vnu = 0.0;
                vsu = 0.0;
            }

            un  = U[i][j] + ((U[i][j+1] - U[i][j]) / Dys[j+1]) * (Dy[j]/2.0);
            us  = U[i][j-1] + ((U[i][j] - U[i][j-1]) / Dys[j]) * (Dy[j-1]/2.0);
            ve  = V[i][j] + ((V[i+1][j] - V[i][j]) / Dxs[i+1]) * (Dx[i]/2.0);
            vw  = V[i-1][j] + ((V[i][j] - V[i-1][j]) / Dxs[i]) * (Dx[i-1]/2.0);

            FU[i][j] = - (ue*ue - uw*uw)/Dxs[i+1] - (un*vnu - us*vsu)/Dy[j]
                       + nu
                       * (
                          (
                           (U[i+1][j] - U[i][j])/Dx[i+1]
                           - (U[i][j] - U[i-1][j])/Dx[i]
                          )/Dxs[i+1]
                          +
                          (
                           (U[i][j+1] - U[i][j])/Dys[j+1]
                           - (U[i][j] - U[i][j-1])/Dys[j]
                          )/Dy[j]
                         );
            FV[i][j] = - (uev*ve - uwv*vw)/Dx[i] - (vn*vn - vs*vs)/Dys[j+1]
                       + nu
                       * (
                          (
                           (V[i+1][j] - V[i][j])/Dxs[i+1]
                           - (V[i][j] - V[i-1][j])/Dxs[i]
                          )/Dx[i]
                          +
                          (
                           (V[i][j+1] - V[i][j])/Dy[j+1]
                           - (V[i][j] - V[i][j-1])/Dy[j]
                          )/Dys[j+1]
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
                     - P[i+1][j]*Dy[j]/Dxs[i+1] - P[i-1][j]*Dy[j]/Dxs[i]
                     - P[i][j+1]*Dx[i]/Dys[j+1] - P[i][j-1]*Dx[i]/Dys[j]
                     + ( (U[i][j] - U[i-1][j])*Dy[j]
                        + (V[i][j] - V[i][j-1])*Dx[i] )
                     *(rho / dt)
                    )
                    /
                    (- Dy[j]/Dxs[i+1] - Dy[j]/Dxs[i]
                     - Dx[i]/Dys[j+1] - Dx[i]/Dys[j]);

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
            U[i][j] = U[i][j] - (dt/rho)*(P[i+1][j] - P[i][j])/Dxs[i+1];
            V[i][j] = V[i][j] - (dt/rho)*(P[i][j+1] - P[i][j])/Dys[j+1];

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
        // if ((t <= 10) || (t % 1000 == 0))
        // {
        //     filerAllSol(t);
        // }
        // End of timestep
        ++t;
    }
    --t;
    filerAllSol();
}
