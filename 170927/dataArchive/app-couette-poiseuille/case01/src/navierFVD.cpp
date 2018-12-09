/*
 * navierFVD.cpp
 *
 *  Created on: 2017-09-27
 *      Author: Syed Ahmad Raza
 */

#include <cmath>        // math functions
#include <iostream>     // functions for input and output to console
#include <string>       // for file name conversion
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
double Uo[nx][ny] = { {0.0} };
double Vo[nx][ny] = { {0.0} };
double FU[nx][ny] = { {0.0} };
double FV[nx][ny] = { {0.0} };
double FUo[nx][ny] = { {0.0} };
double FVo[nx][ny] = { {0.0} };
double P[nx][ny] = { {0.0} };

int t = 0;
double uMaxChange = 1.0;

void initial()
{
    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < ny; ++j)
        {
            U[i][j] = uin;
            V[i][j] = vin;
            P[i][j] = pconst;
        }
    }
}

void velBoundary()
{
    // Applying the velocity boundary conditions
    // Top boundary
    for (int i = 0; i < nx; ++i)
    {
        // U[i][ny-1] = -U[i][ny-2];       // Dirichlet
        U[i][ny-1] = -U[i][ny-2] + 2*uin;       // Dirichlet
        V[i][ny-1] = 0.0;               // Dirichlet
    }
    // Bottom boundary
    for (int i = 0; i < nx; ++i)
    {
        U[i][0] = -U[i][1];             // Dirichlet
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
        P[0][j] = P[1][j];              // Neumann
        // P[0][j] = pconst;               // Dirichlet
    }
    // Right boundary
    for (int j = 0; j < ny; ++j)
    {
        P[nx-1][j] = pconst;            // Dirichlet
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

void momentum(string velScheme)
{
    for (int i = 1; i < nx-1; ++i)
    {
        for (int j = 1; j < ny-1; ++j)
        {
            double ue, uw, un, us, ve, vw, vn, vs;

            if (velScheme == "central")
            {
                ue = (U[i][j] + U[i+1][j]) / 2.0;
                uw = (U[i-1][j] + U[i][j]) / 2.0;
                vn = (V[i][j] + V[i][j+1]) / 2.0;
                vs = (V[i][j-1] + V[i][j]) / 2.0;
            }
            else if (velScheme == "upwind")
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
            else if (velScheme == "quick")
            {
                if (U[i][j] > 0.0)
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
                else
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
                if (V[i][j] > 0.0)
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
                else
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
                cout << '\n' << "Undefined velocity scheme selected" << '\n';
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
    double pMaxChange = 1.0;
    while (m <= maxPressIterations && pMaxChange > pResidual)
    {
        pMaxChange = 0.0;
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

                if (abs(pTemp - P[i][j]) > pMaxChange)
                {
                    pMaxChange = abs(pTemp - P[i][j]);
                }

                // Applying the SOR method through following formula:
                P[i][j] = P[i][j] + omega * (pTemp - P[i][j]);
                // P[i][j] = (1.0 - omega)*P[i][j] + omega * pTemp; //same
            }
        }
        pressBoundary();
        ++m;
    }
    cout << "Number of pressure iterations\t= " << --m << '\n';
    cout << "Maximum change in pressure\t= " << pMaxChange << '\n';
}

void velUpdater()
{
    for (int i = 1; i < nx-1; ++i)
    {
        for (int j = 1; j < ny-1; ++j)
        {
            U[i][j] = U[i][j] - (dt / rho) * (P[i+1][j] - P[i][j])/dx;
            V[i][j] = V[i][j] - (dt / rho) * (P[i][j+1] - P[i][j])/dy;

            if (abs(Uo[i][j] - U[i][j]) > uMaxChange)
            {
                uMaxChange = abs(Uo[i][j] - U[i][j]);
            }
            // if (abs(Vo[i][j] - V[i][j]) > vMaxChange)
            // {
            //     vMaxChange = abs(Vo[i][j] - V[i][j]);
            // }
            Uo[i][j] = U[i][j];
            Vo[i][j] = V[i][j];
        }
    }
}

void navierFVD()
{
    gridder();

    cout << "\n\n" << "Time step\t\t\t= " << t << '\n';

    initial();
    velBoundary();
    pressBoundary();
    filerAllSol("../data/navierFVD", t);

    cout << '\n'; // End of time-step
    t = 1;

    while   (t <= 10 || (t <= maxTimesteps && uMaxChange > uResidual))
    {
        uMaxChange = 0.0;

        cout << "Time step\t\t\t= " << t << '\n';

        momentum("quick");
        momentumTimeStep();
        velBoundary();
        pressure();
        velUpdater();
        velBoundary();

        cout << "Maximum change in U velocity\t= " << uMaxChange << '\n';

        if ((t <= 2) || (t % 1000 == 0))
        {
            filerAllSol("../data/navierFVD", t);
        }

        cout << '\n'; // End of time-step
        ++t;
    }
    cout << "\nFinal time-step is " << --t << '\n';
    filerAllSol("../data/navierFVD");
}
