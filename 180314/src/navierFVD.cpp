/*
 * navierFVD.cpp
 *
 *  Created on: 2017-09-27
 *      Author: Syed Ahmad Raza
 */

#include <cmath>        // math functions
// #include <algorithm>    // functions for ranges of elements (arrays)
#include <iostream>     // functions for input and output to console
#include <ctime>        // to time the script
#include <chrono>       // to measure and display the elapsed time

#include "constants.h"
#include "gridder.h"
#include "filers.h"
#include "printers.h"
#include "navierFVD.h"

using namespace std;

int t               = 0;
int pIter           = -1;
double uChangeMax   = 1.0;
double vChangeMax   = 1.0;
double pChangeMax   = 1.0;
double mChangeMax   = 1.0;
int mChangeMax_i    = -1;
int mChangeMax_j    = -1;
int pChangeMax_i    = -1;
int pChangeMax_j    = -1;
int uChangeMax_i    = -1;
int uChangeMax_j    = -1;
int vChangeMax_i    = -1;
int vChangeMax_j    = -1;

double scriptRunningTime = 0.0;

double U[nx][ny]    = { {0.0} };
double V[nx][ny]    = { {0.0} };
double Uo[nx][ny]   = { {0.0} };
double Vo[nx][ny]   = { {0.0} };
double FU[nx][ny]   = { {0.0} };
double FV[nx][ny]   = { {0.0} };
double FU1[nx][ny]  = { {0.0} };
double FV1[nx][ny]  = { {0.0} };
double FU2[nx][ny]  = { {0.0} };
double FV2[nx][ny]  = { {0.0} };
double P[nx][ny]    = { {0.0} };
double MC[nx][ny]   = { {0.0} };
double PC[nx][ny]   = { {0.0} };

void initial()
{
    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < ny; ++j)
        {
            U[i][j] = uInitial;
            V[i][j] = vInitial;
            P[i][j] = pInitial;
        }
    }
}

void velBoundary()
{
    // Applying the velocity boundary conditions
    for (int j = 0; j < ny; ++j)
    {
        // Left boundary
        U[0][j] = 0.0;                  // Dirichlet
        V[0][j] = -V[1][j];             // Dirichlet
        // V[0][j] = -V[1][j] + 2*vIn;  // Dirichlet

        // Right boundary
        U[nx-1][j] = 0.0;               // Dirichlet
        V[nx-1][j] = -V[nx-2][j];       // Dirichlet
        // U[nx-1][j] = U[nx-2][j];        // Neumann
        // V[nx-1][j] = V[nx-2][j];        // Neumann
    }
    for (int i = 0; i < nx; ++i)
    {
        // Top boundary
        U[i][ny-1] = -U[i][ny-2] + 2*uIn;   // Dirichlet
        V[i][ny-1] = 0.0;                   // Dirichlet
        // U[i][ny-1] = uIn;                   // Dirichlet
        // U[i][ny-1] = -U[i][ny-2];           // Dirichlet

        // Bottom boundary
        U[i][0] = -U[i][1];             // Dirichlet
        V[i][0] = 0.0;                  // Dirichlet
        // U[i][0] = -U[i][1] + 2*uIn;  // Dirichlet
    }
}

void pressBoundary()
{
    // Applying the pressure boundary conditions
    for (int j = 0; j < ny; ++j)
    {
        // Left boundary
        P[0][j] = P[1][j];              // Neumann
        // P[0][j] = pIn;                  // Dirichlet

        // Right boundary
        P[nx-1][j] = P[nx-2][j];        // Neumann
        // P[nx-1][j] = pOut;              // Dirichlet
    }
    for (int i = 0; i < nx; ++i)
    {
        // Top boundary
        P[i][ny-1] = P[i][ny-2];        // Neumann
        // P[0][ny-1] = pInitial;          // fixed pressure point

        // Bottom boundary
        P[i][0] = P[i][1];              // Neumann
    }
}

void momentum(string velocityScheme)
{
    double ue, uw, un, us, uev, uwv, ve, vw, vn, vs, vnu, vsu;
    for (int i = 1; i < nx-1; ++i)
    {
        for (int j = 1; j < ny-1; ++j)
        {
            // // Central scheme
            if (velocityScheme == "ct")
            {
                ue  = 0.5*( U[i][j] + U[i+1][j] );
                uw  = 0.5*( U[i-1][j] + U[i][j] );
                vnu = 0.5*( Dx[i+1]*V[i][j]   + Dx[i]*V[i+1][j]   ) / Dxs[i+1];
                vsu = 0.5*( Dx[i+1]*V[i][j-1] + Dx[i]*V[i+1][j-1] ) / Dxs[i+1];
                ve  = 0.5*V[i][j] + ((V[i+1][j] - V[i][j]) / Dxs[i+1])*Dx[i];
                vw  = 0.5*V[i-1][j] + ((V[i][j] - V[i-1][j]) / Dxs[i])*Dx[i-1];

                vn  = 0.5*( V[i][j] + V[i][j+1] );
                vs  = 0.5*( V[i][j-1] + V[i][j] );
                uev = 0.5*( Dy[j+1]*U[i][j]   + Dy[j]*U[i][j+1]   ) / Dys[j+1];
                uwv = 0.5*( Dy[j+1]*U[i-1][j] + Dy[j]*U[i-1][j+1] ) / Dys[j+1];
                un  = 0.5*U[i][j] + ((U[i][j+1] - U[i][j]) / Dys[j+1])*Dy[j];
                us  = 0.5*U[i][j-1] + ((U[i][j] - U[i][j-1]) / Dys[j])*Dy[j-1];
            }
            // Upwind scheme
            else if (velocityScheme == "up")
            {
                if (U[i][j] >= 0.0)
                {
                    ue  = U[i][j];
                    uw  = U[i-1][j];
                    uev = U[i][j];
                    uwv = U[i-1][j];
                    un  = U[i][j];
                    us  = U[i][j-1];
                }
                else
                {
                    ue  = U[i+1][j];
                    uw  = U[i][j];
                    uev = U[i][j+1];
                    uwv = U[i-1][j+1];
                    un  = U[i][j+1];
                    us  = U[i][j];
                }
                if (V[i][j] >= 0.0)
                {
                    vn  = V[i][j];
                    vs  = V[i][j-1];
                    vnu = V[i][j];
                    vsu = V[i][j-1];
                    ve  = V[i][j];
                    vw  = V[i-1][j];
                }
                else
                {
                    vn  = V[i][j+1];
                    vs  = V[i][j];
                    vnu = V[i+1][j];
                    vsu = V[i+1][j-1];
                    ve  = V[i+1][j];
                    vw  = V[i][j];
                }
            }
            // QUICK scheme
            else if (velocityScheme == "qk")
            {
                if (U[i][j] >= 0.0)     // positive
                {
                    if (i == 1 || i == nx-2 || j == 1 || j == ny-2)
                    {
                        ue  = U[i][j];
                        uw  = U[i-1][j];
                        uev = U[i][j];
                        uwv = U[i-1][j];
                        un  = U[i][j];
                        us  = U[i][j-1];
                    }
                    else
                    {
                        ue  = 0.5*( U[i][j]     + U[i+1][j]   )
                              - 0.125 * Dx[i+1] * Dx[i+1]       / Dxs[i+1]
                              * ( ( U[i+1][j]   - U[i][j]     ) / Dx[i+1]
                                - ( U[i][j]     - U[i-1][j]   ) / Dx[i]     );

                        uw  = 0.5*( U[i-1][j]   + U[i][j]     )
                              - 0.125 * Dx[i]   * Dx[i]         / Dxs[i]
                              * ( ( U[i][j]     - U[i-1][j]   ) / Dx[i]
                                - ( U[i-1][j]   - U[i-2][j]   ) / Dx[i-1]   );

                        uev = 0.5*( U[i][j]     + U[i][j+1]   )
                              - 0.125 * Dys[j+1]* Dys[j+1]      / Dy[j]
                              * ( ( U[i][j+1]   - U[i][j]     ) / Dys[j+1]
                                - ( U[i][j]     - U[i][j-1]   ) / Dys[j]    );

                        uwv = 0.5*( U[i-1][j]   + U[i-1][j+1]   )
                              - 0.125 * Dys[j+1]* Dys[j+1]      / Dy[j]
                              * ( ( U[i-1][j+1] - U[i-1][j]   ) / Dys[j+1]
                                - ( U[i-1][j]   - U[i-1][j-1] ) / Dys[j]    );

                        un  = 0.5*( U[i][j]     + U[i][j+1]   )
                              - 0.125 * Dys[j+1]* Dys[j+1]      / Dy[j]
                              * ( ( U[i][j+1]   - U[i][j]     ) / Dys[j+1]
                                - ( U[i][j]     - U[i][j-1]   ) / Dys[j]    );

                        us  = 0.5*( U[i][j-1]   + U[i][j]     )
                              - 0.125 * Dys[j]  * Dys[j]        / Dy[j-1]
                              * ( ( U[i][j]     - U[i][j-1]   ) / Dys[j]
                                - ( U[i][j-1]   - U[i][j-2]   ) / Dys[j-1]  );
                    }
                }
                else                    // negative
                {
                    if (i == 1 || i == nx-2 || j == 1 || j == ny-2)
                    {
                        ue  = U[i+1][j];
                        uw  = U[i][j];
                        uev = U[i][j+1];
                        uwv = U[i-1][j+1];
                        un  = U[i][j+1];
                        us  = U[i][j];
                    }
                    else
                    {
                        ue  = 0.5*( U[i][j]     + U[i+1][j]   )
                              - 0.125 * Dx[i+1] * Dx[i+1]       / Dxs[i+2]
                              * ( ( U[i+2][j]   - U[i+1][j]   ) / Dx[i+2]
                                - ( U[i+1][j]   - U[i][j]     ) / Dx[i+1]   );

                        uw  = 0.5*( U[i-1][j]   + U[i][j]     )
                              - 0.125 * Dx[i]   * Dx[i]         / Dxs[i+1]
                              * ( ( U[i+1][j]   - U[i][j]     ) / Dx[i+1]
                                - ( U[i][j]     - U[i-1][j]   ) / Dx[i]     );

                        uev = 0.5*( U[i][j]     + U[i][j+1]   )
                              - 0.125 * Dys[j+1]* Dys[j+1]      / Dy[j+1]
                              * ( ( U[i][j+2]   - U[i][j+1]   ) / Dys[j+2]
                                - ( U[i][j+1]   - U[i][j]     ) / Dys[j+1]  );

                        uwv = 0.5*( U[i-1][j]   + U[i-1][j+1] )
                              - 0.125 * Dys[j+1]* Dys[j+1]      / Dy[j+1]
                              * ( ( U[i-1][j+2] - U[i-1][j+1] ) / Dys[j+2]
                                - ( U[i-1][j+1] - U[i-1][j]   ) / Dys[j+1]  );

                        un  = 0.5*( U[i][j]     + U[i][j+1]   )
                              - 0.125 * Dys[j+1]* Dys[j+1]      / Dy[j+1]
                              * ( ( U[i][j+2]   - U[i][j+1]   ) / Dys[j+2]
                                - ( U[i][j+1]   - U[i][j]     ) / Dys[j+1]  );

                        us  = 0.5*( U[i][j-1]   + U[i][j]     )
                              - 0.125 * Dys[j]  * Dys[j]        / Dy[j]
                              * ( ( U[i][j+1]   - U[i][j]     ) / Dys[j+1]
                                - ( U[i][j]     - U[i][j-1]   ) / Dys[j]    );
                    }
                }
                if (V[i][j] >= 0.0)     // positive
                {
                    if (i == 1 || i == nx-2 || j == 1 || j == ny-2)
                    {
                        vn  = V[i][j];
                        vs  = V[i][j-1];
                        vnu = V[i][j];
                        vsu = V[i][j-1];
                        ve  = V[i][j];
                        vw  = V[i-1][j];
                    }
                    else
                    {
                        vn  = 0.5*( V[i][j]     + V[i][j+1]   )
                              - 0.125 * Dy[j+1] * Dy[j+1]       / Dys[j+1]
                              * ( ( V[i][j+1]   - V[i][j]     ) / Dy[j+1]
                                - ( V[i][j]     - V[i][j-1]   ) / Dy[j]     );

                        vs  = 0.5*( V[i][j-1]   + V[i][j]     )
                              - 0.125 * Dy[j]   * Dy[j]         / Dys[j]
                              * ( ( V[i][j]     - V[i][j-1]   ) / Dy[j]
                                - ( V[i][j-1]   - V[i][j-2]   ) / Dy[j-1]   );

                        vnu = 0.5*( V[i][j]     + V[i+1][j]   )
                              - 0.125 * Dxs[i+1]* Dxs[i+1]      / Dx[i]
                              * ( ( V[i+1][j]   - V[i][j]     ) / Dxs[i+1]
                                - ( V[i][j]     - V[i-1][j]   ) / Dxs[i]    );

                        vsu = 0.5*( V[i][j-1]   + V[i+1][j-1] )
                              - 0.125 * Dxs[i+1]* Dxs[i+1]      / Dx[i]
                              * ( ( V[i+1][j-1] - V[i][j-1]   ) / Dxs[i+1]
                                - ( V[i][j-1]   - V[i-1][j-1] ) / Dxs[i]    );

                        ve  = 0.5*( V[i][j]     + V[i+1][j]   )
                              - 0.125 * Dxs[i+1]* Dxs[i+1]      / Dx[i]
                              * ( ( V[i+1][j]   - V[i][j]     ) / Dxs[i+1]
                                - ( V[i][j]     - V[i-1][j]   ) / Dxs[i]    );

                        vw  = 0.5*( V[i-1][j]   + V[i][j]     )
                              - 0.125 * Dxs[i]  * Dxs[i]        / Dx[i-1]
                              * ( ( V[i][j]     - V[i-1][j]   ) / Dxs[i]
                                - ( V[i-1][j]   - V[i-2][j]   ) / Dxs[i-1]  );
                    }
                }
                else                    // negative
                {
                    if (i == 1 || i == nx-2 || j == 1 || j == ny-2)
                    {
                        vn  = V[i][j+1];
                        vs  = V[i][j];
                        vnu = V[i+1][j];
                        vsu = V[i+1][j-1];
                        ve  = V[i+1][j];
                        vw  = V[i][j];
                    }
                    else
                    {
                        vn  = 0.5*( V[i][j]     + V[i][j+1]   )
                              - 0.125 * Dy[j+1] * Dy[j+1]       / Dys[j+2]
                              * ( ( V[i][j+2]   - V[i][j+1]   ) / Dy[j+2]
                                - ( V[i][j+1]   - V[i][j]     ) / Dy[j+1]   );

                        vs  = 0.5*( V[i][j-1]   + V[i][j]     )
                              - 0.125 * Dy[j]   * Dy[j]         / Dys[j+1]
                              * ( ( V[i][j+1]   - V[i][j]     ) / Dy[j+1]
                                - ( V[i][j]     - V[i][j-1]   ) / Dy[j]     );

                        vnu = 0.5*( V[i][j]     + V[i+1][j]   )
                              - 0.125 * Dxs[i+1]* Dxs[i+1]      / Dx[i+1]
                              * ( ( V[i+2][j]   - V[i+1][j]   ) / Dxs[i+2]
                                - ( V[i+1][j]   - V[i][j]     ) / Dxs[i+1]  );

                        vsu = 0.5*( V[i][j-1]   + V[i+1][j-1] )
                              - 0.125 * Dxs[i+1]* Dxs[i+1]      / Dx[i+1]
                              * ( ( V[i+2][j-1] - V[i+1][j-1] ) / Dxs[i+2]
                                - ( V[i+1][j-1] - V[i][j-1]   ) / Dxs[i+1]  );

                        ve  = 0.5*( V[i][j]     + V[i+1][j]   )
                              - 0.125 * Dxs[i+1]* Dxs[i+1]      / Dx[i+1]
                              * ( ( V[i+2][j]   - V[i+1][j]   ) / Dxs[i+2]
                                - ( V[i+1][j]   - V[i][j]     ) / Dxs[i+1]  );

                        vw  = 0.5*( V[i-1][j]   + V[i][j]     )
                              - 0.125 * Dxs[i]  * Dxs[i]        / Dx[i]
                              * ( ( V[i+1][j]   - V[i][j]     ) / Dxs[i+1]
                                - ( V[i][j]     - V[i-1][j]   ) / Dxs[i]    );
                    }
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
                ve  = 0.0;
                vw  = 0.0;

                vn  = 0.0;
                vs  = 0.0;
                vnu = 0.0;
                vsu = 0.0;
                un  = 0.0;
                us  = 0.0;
            }

            FU[i][j]
            =   - ( ue  * ue  - uw  * uw  )     / Dxs[i+1]
                - ( un  * vnu - us  * vsu )     / Dy[j]
                + nu
                * (
                   (  ( U[i+1][j] - U[i][j]   ) / Dx[i+1]
                    - ( U[i][j]   - U[i-1][j] ) / Dx[i]     ) / Dxs[i+1]
                   +
                   (  ( U[i][j+1] - U[i][j]   ) / Dys[j+1]
                    - ( U[i][j]   - U[i][j-1] ) / Dys[j]    ) / Dy[j]
                  );

            FV[i][j]
            =   - ( uev * ve  - uwv * vw  )     / Dx[i]
                - ( vn  * vn  - vs  * vs  )     / Dys[j+1]
                + nu
                * (
                   (  ( V[i+1][j] - V[i][j]   ) / Dxs[i+1]
                    - (V[i][j]    - V[i-1][j] ) / Dxs[i]    ) / Dx[i]
                   +
                   (  ( V[i][j+1] - V[i][j]   ) / Dy[j+1]
                    - ( V[i][j]   - V[i][j-1] ) / Dy[j]     ) / Dys[j+1]
                  );
        }
    }
}

void timeStepNSchemer()
{
    for (int i = 1; i < nx-1; ++i)
    {
        for (int j = 1; j < ny-1; ++j)
        {
            FU[i][j] = dt * FU[i][j];
            FV[i][j] = dt * FV[i][j];

            if (t == 1)         // Euler scheme for the first timestep
            {
                U[i][j] += FU[i][j];
                V[i][j] += FV[i][j];
            }
            else if (t == 2)    // Adams-Bashforth 2nd order for 2nd timestep
            {
                U[i][j] += 1.5*FU[i][j] - 0.5*FU1[i][j];
                V[i][j] += 1.5*FV[i][j] - 0.5*FV1[i][j];
            }
            else                // Adams-Bashforth 3rd order subsequently
            {
                U[i][j] += (23.0*FU[i][j] -16.0*FU1[i][j] +5.0*FU2[i][j])/12.0;
                V[i][j] += (23.0*FV[i][j] -16.0*FV1[i][j] +5.0*FV2[i][j])/12.0;
            }
            FU2[i][j] = FU1[i][j];
            FU1[i][j] = FU[i][j];
            FV2[i][j] = FV1[i][j];
            FV1[i][j] = FV[i][j];
        }
    }
}

void pressure()
{
    int h = 0;
    double pNew = 0.0;
    pChangeMax = 1.0;
    while (h <= maxPressIters && pChangeMax > pResidual)
    {
        pChangeMax = 0.0;       // maximum pressure residual during timestep
        mChangeMax = 0.0;       // maximum mass residual during timestep
        double pChange = 0.0;   // pressure residual during current iteration
        double mChange = 0.0;   // mass residual during current iteration
        for (int i = 1; i < nx-1; ++i)
        {
            for (int j = 1; j < ny-1; ++j)
            {
                mChange =   ( U[i][j] - U[i-1][j] ) * Dy[j]
                          + ( V[i][j] - V[i][j-1] ) * Dx[i];

                if (abs(mChange) > mChangeMax)
                {
                  mChangeMax = abs(mChange);
                  mChangeMax_i = i;
                  mChangeMax_j = j;
                }

                pNew =
                    (
                     - P[i+1][j] * Dy[j] / Dxs[i+1] - P[i-1][j] * Dy[j] / Dxs[i]
                     - P[i][j+1] * Dx[i] / Dys[j+1] - P[i][j-1] * Dx[i] / Dys[j]
                     + rho * mChange / dt
                    )
                    /
                    (- Dy[j] / Dxs[i+1] - Dy[j] / Dxs[i]
                     - Dx[i] / Dys[j+1] - Dx[i] / Dys[j] );

                pChange = abs(pNew - P[i][j]);
                if (pChange > pChangeMax)
                {
                    pChangeMax = pChange;
                    pChangeMax_i = i;
                    pChangeMax_j = j;
                }

                // Applying the SOR method through following formula:
                P[i][j] = P[i][j] + omega * (pNew - P[i][j]);
                // P[i][j] = (1.0 - omega)*P[i][j] + omega * pNew; //same
            }
        }
        pressBoundary();
        ++h;
    }
    pIter = --h;
}

void velUpdater()
{
    for (int i = 1; i < nx-1; ++i)
    {
        for (int j = 1; j < ny-1; ++j)
        {
            U[i][j] = U[i][j] - (dt/rho) * ( P[i+1][j] - P[i][j] ) / Dxs[i+1];
            V[i][j] = V[i][j] - (dt/rho) * ( P[i][j+1] - P[i][j] ) / Dys[j+1];

            // uChangeMax += abs(Uo[i][j] - U[i][j]);
            // vChangeMax += abs(Vo[i][j] - V[i][j]);
            if (abs(Uo[i][j] - U[i][j]) > uChangeMax)
            {
                uChangeMax = abs(Uo[i][j] - U[i][j]);
                uChangeMax_i = i;
                uChangeMax_j = j;
            }
            if (abs(Vo[i][j] - V[i][j]) > vChangeMax)
            {
                vChangeMax = abs(Vo[i][j] - V[i][j]);
                vChangeMax_i = i;
                vChangeMax_j = j;
            }
            Uo[i][j] = U[i][j];
            Vo[i][j] = V[i][j];
        }
    }
    // uChangeMax = uChangeMax / static_cast<double>(nx * ny);
    // vChangeMax = vChangeMax / static_cast<double>(nx * ny);
}

void navierFVD()
{
    chrono::steady_clock::time_point tStart = chrono::steady_clock::now();
    initial();
    velBoundary();
    pressBoundary();
    filerAllSol(t);
    progressFileCreator();
    // End of timestep
    t = 1;
    while   (t <= 10 || (t <= maxTimesteps && uChangeMax > uResidual))
    {
        uChangeMax = 0.0;
        vChangeMax = 0.0;

        momentum(velScheme);
        timeStepNSchemer();
        velBoundary();
        pressure();
        velUpdater();
        velBoundary();
        progressPrinter();
        progressFiler();
        if (t % 100000 == 0)
        {
            filerAllSol(t);
        }
        // End of timestep
        ++t;
    }
    --t;
    filerAllSol();
    chrono::steady_clock::time_point tEnd = chrono::steady_clock::now();
    scriptRunningTime =
        chrono::duration_cast<chrono::minutes>(tEnd - tStart).count();
}
