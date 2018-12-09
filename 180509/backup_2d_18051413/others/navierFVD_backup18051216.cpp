/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*\
|                ____       |                                                 |
|       //\     |    \      |  Filename:    navierFVD.cpp                     |
|      //  \    |___ /      |  Created on:  2017-09-27                        |
|     //====\   |  \        |  Last update: 2018-04-25                        |
|    //      \  |   \_      |  Author:      Syed Ahmad Raza                   |
\*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


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
using namespace Numeric_lib;

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

Matrix<double,2> U(nx,ny);
Matrix<double,2> V(nx,ny);
Matrix<double,2> Uo(nx,ny);
Matrix<double,2> Vo(nx,ny);
Matrix<double,2> FU(nx,ny);
Matrix<double,2> FV(nx,ny);
Matrix<double,2> FU1(nx,ny);
Matrix<double,2> FV1(nx,ny);
Matrix<double,2> FU2(nx,ny);
Matrix<double,2> FV2(nx,ny);
Matrix<double,2> P(nx,ny);
Matrix<double,2> MC(nx,ny);
Matrix<double,2> PC(nx,ny);

void initial()
{
    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < ny; ++j)
        {
            U(i,j) = uInitial;
            V(i,j) = vInitial;
            P(i,j) = pInitial;
        }
    }
}

void velBoundary()
{
    // Applying the velocity boundary conditions
    for (int j = 0; j <= ny-1; ++j)
    {
        // Left boundary
        U(1,j) = 0.0;
        U(0,j) = U(1,j);
        V(1,j) = -V(2,j);
        V(0,j) = V(1,j);

        // V(0,j) = -V(1,j);               // last correct
        // V(0,j) = -V(1,j) + 2*vIn;       // moving boundary

        // Right boundary
        U(nx-3,j) = 0.0;
        U(nx-2,j) = U(nx-3,j);
        // U(nx-1,j) = U(nx-2,j);
        V(nx-2,j) = -V(nx-3,j);
        V(nx-1,j) = V(nx-2,j);
        // V(nx-1,j) = -V(nx-2,j);         // last correct
        // U(nx-1,j) = U(nx-2,j);          // Neumann
        // V(nx-1,j) = V(nx-2,j);          // Neumann
    }
    for (int i = 0; i <= nx-1; ++i)
    {
        // Bottom boundary
        U(i,1) = -U(i,2);
        U(i,0) = U(i,1);
        V(i,1) = 0.0;
        V(i,0) = V(i,1);
        // U(i,0) = -U(i,1);               // last correct
        // U(i,0) = -U(i,1) + 2*uIn;       // Dirichlet

        // Top boundary
        U(i,ny-2) = -U(i,ny-3) + 2*uIn;
        U(i,ny-1) = U(i,ny-2);
        V(i,ny-3) = 0.0;
        V(i,ny-2) = V(i,ny-3);

        // V(i,ny-1) = V(i,ny-2);
        // U(i,ny-1) = -U(i,ny-2) + 2*uIn; // last correct
        // U(i,ny-1) = -U(i,ny-2);         // stationary boundary
    }
}

void pressBoundary()
{
    // Applying the pressure boundary conditions
    for (int j = 1; j <= ny-2; ++j)
    {
        // Left boundary
        P(1,j) = P(2,j);                // Neumann
        // P(0,j) = pIn;                   // Dirichlet

        // Right boundary
        P(nx-2,j) = P(nx-3,j);          // Neumann
        // P(nx-1,j) = pOut;               // Dirichlet
    }
    for (int i = 1; i <= nx-2; ++i)
    {
        // Top boundary
        P(i,ny-2) = P(i,ny-3);          // Neumann
        // P(0,ny-1) = pInitial;           // fixed pressure point

        // Bottom boundary
        P(i,1) = P(i,2);                // Neumann
    }
}

// QUICK scheme
void quick()
{
    double ue, uw, un, us, ve, vw, vn, vs,
           uev, uwv, vnu, vsu;
    // double uex, uwx, vnx, vsx, vny, vsy, uey, uwy;
    for (int i = 2; i <= nx-3; ++i)
    {
        for (int j = 2; j <= ny-3; ++j)
        {
            // ue
            // uex = 0.5 * (  U(i,j)   * Dx(i+1) / Dxs(i+1)
            //              + U(i+1,j) * Dx(i)   / Dxs(i+1) );
            // if (uex > 0.0)          // positive
            if (U(i,j) > 0.0)          // positive
            {
                ue  = 0.5*( U(i,j)     + U(i+1,j)   )
                      // - 0.125 * Dx(i+1) * Dx(i+1)       / Dxs(i+1)
                      // * ( ( U(i+1,j)   - U(i,j)     ) / Dx(i+1)
                      //   - ( U(i,j)     - U(i-1,j)   ) / Dx(i)     );
                      ;
            }
            else
            {
                ue  = 0.5*( U(i,j)     + U(i+1,j)   )
                      // - 0.125 * Dx(i+1) * Dx(i+1)       / Dxs(i+2)
                      // * ( ( U(i+2,j)   - U(i+1,j)   ) / Dx(i+2)
                      //   - ( U(i+1,j)   - U(i,j)     ) / Dx(i+1)   );
                      ;
            }
            // uw
            // uwx = 0.5 * (  U(i,j)   * Dx(i-1) / Dxs(i)
            //              + U(i-1,j) * Dx(i)   / Dxs(i) );
            // if (uwx > 0.0)          // positive
            if (U(i,j) > 0.0)          // positive
            {
                uw  = 0.5*( U(i-1,j)   + U(i,j)     )
                      // - 0.125 * Dx(i)   * Dx(i)         / Dxs(i)
                      // * ( ( U(i,j)     - U(i-1,j)   ) / Dx(i)
                      //   - ( U(i-1,j)   - U(i-2,j)   ) / Dx(i-1)   );
                      ;
            }
            else
            {
                uw  = 0.5*( U(i-1,j)   + U(i,j)     )
                      // - 0.125 * Dx(i)   * Dx(i)         / Dxs(i+1)
                      // * ( ( U(i+1,j)   - U(i,j)     ) / Dx(i+1)
                      //   - ( U(i,j)     - U(i-1,j)   ) / Dx(i)     );
                      ;
            }
            // un
            // vnx = 0.5 * ( V[i][j] + V[i+1][j] );
            // if (vnx > 0.0)          // positive
            if (V(i,j) > 0.0)          // positive
            {
                un  = 0.5*( U(i,j)     + U(i,j+1)   )
                      // - 0.125 * Dys(j+1)* Dys(j+1)      / Dy(j)
                      // * ( ( U(i,j+1)   - U(i,j)     ) / Dys(j+1)
                      //   - ( U(i,j)     - U(i,j-1)   ) / Dys(j)    );
                      ;
            }
            else
            {
                un  = 0.5*( U(i,j)     + U(i,j+1)   )
                      // - 0.125 * Dys(j+1)* Dys(j+1)      / Dy(j+1)
                      // * ( ( U(i,j+2)   - U(i,j+1)   ) / Dys(j+2)
                      //   - ( U(i,j+1)   - U(i,j)     ) / Dys(j+1)  );
                      ;
            }
            // us
            // vsx = 0.5 * ( V[i][j-1] + V[i+1][j-1] );
            // if (vsx > 0.0)          // positive
            if (V(i,j) > 0.0)          // positive
            {
                us  = 0.5*( U(i,j-1)   + U(i,j)     )
                      // - 0.125 * Dys(j)  * Dys(j)        / Dy(j-1)
                      // * ( ( U(i,j)     - U(i,j-1)   ) / Dys(j)
                      //   - ( U(i,j-1)   - U(i,j-2)   ) / Dys(j-1)  );
                      ;
            }
            else
            {
                us  = 0.5*( U(i,j-1)   + U(i,j)     )
                      // - 0.125 * Dys(j)  * Dys(j)        / Dy(j)
                      // * ( ( U(i,j+1)   - U(i,j)     ) / Dys(j+1)
                      //   - ( U(i,j)     - U(i,j-1)   ) / Dys(j)    );
                      ;
            }

            // Old variables
            // vnu
            // if (0.5*(U[i][j] + U[i][j+1]) >= 0.0)       // positive
            if (U(i,j) > 0.0)          // positive
            {
                vnu = 0.5*( V(i,j)     + V(i+1,j)   )
                      // - 0.125 * Dxs(i+1)* Dxs(i+1)      / Dx(i)
                      // * ( ( V(i+1,j)   - V(i,j)     ) / Dxs(i+1)
                      //   - ( V(i,j)     - V(i-1,j)   ) / Dxs(i)    );
                      ;
            }
            else
            {
                vnu = 0.5*( V(i,j)     + V(i+1,j)   )
                      // - 0.125 * Dxs(i+1)* Dxs(i+1)      / Dx(i+1)
                      // * ( ( V(i+2,j)   - V(i+1,j)   ) / Dxs(i+2)
                      //   - ( V(i+1,j)   - V(i,j)     ) / Dxs(i+1)  );
                      ;
            }
            // vsu
            // if (0.5*(U[i][j-1] + U[i][j]) >= 0.0)       // positive
            if (U(i,j) > 0.0)          // positive
            {
                vsu = 0.5*( V(i,j-1)   + V(i+1,j-1) )
                      // - 0.125 * Dxs(i+1)* Dxs(i+1)      / Dx(i)
                      // * ( ( V(i+1,j-1) - V(i,j-1)   ) / Dxs(i+1)
                      //   - ( V(i,j-1)   - V(i-1,j-1) ) / Dxs(i)    );
                      ;
            }
            else
            {
                vsu = 0.5*( V(i,j-1)   + V(i+1,j-1) )
                      // - 0.125 * Dxs(i+1)* Dxs(i+1)      / Dx(i+1)
                      // * ( ( V(i+2,j-1) - V(i+1,j-1) ) / Dxs(i+2)
                      //   - ( V(i+1,j-1) - V(i,j-1)   ) / Dxs(i+1)  );
                      ;
            }

            // vn
            // if (0.5*(V[i][j] + V[i][j+1]) >= 0.0)       // old condition
            // vny = 0.5 * (  V(i,j)   * Dy(j+1) / Dys(j+1)
            //              + V(i,j+1) * Dy(j)   / Dys(j+1) );
            // if (vny > 0.0)          // positive
            if (V(i,j) > 0.0)          // positive
            {
                vn  = 0.5*( V(i,j)     + V(i,j+1)   )
                      // - 0.125 * Dy(j+1) * Dy(j+1)       / Dys(j+1)
                      // * ( ( V(i,j+1)   - V(i,j)     ) / Dy(j+1)
                      //   - ( V(i,j)     - V(i,j-1)   ) / Dy(j)     );
                      ;
            }
            else
            {
                vn  = 0.5*( V(i,j)     + V(i,j+1)   )
                      // - 0.125 * Dy(j+1) * Dy(j+1)       / Dys(j+2)
                      // * ( ( V(i,j+2)   - V(i,j+1)   ) / Dy(j+2)
                      //   - ( V(i,j+1)   - V(i,j)     ) / Dy(j+1)   );
                      ;
            }
            // vs
            // vsy = 0.5 * (  V(i,j)   * Dy(j-1) / Dys(j)
            //              + V(i,j-1) * Dy(j)   / Dys(j) );
            // if (vsy > 0.0)          // positive
            if (V(i,j) > 0.0)          // positive
            {
                vs  = 0.5*( V(i,j-1)   + V(i,j)     )
                      // - 0.125 * Dy(j)   * Dy(j)         / Dys(j)
                      // * ( ( V(i,j)     - V(i,j-1)   ) / Dy(j)
                      //   - ( V(i,j-1)   - V(i,j-2)   ) / Dy(j-1)   );
                      ;
            }
            else
            {
                vs  = 0.5*( V(i,j-1)   + V(i,j)     )
                      // - 0.125 * Dy(j)   * Dy(j)         / Dys(j+1)
                      // * ( ( V(i,j+1)   - V(i,j)     ) / Dy(j+1)
                      //   - ( V(i,j)     - V(i,j-1)   ) / Dy(j)     );
                      ;
            }
            // ve
            // uey = 0.5 * ( U[i][j] + U[i][j+1] );
            // if (uey > 0.0)          // positive
            if (U(i,j) > 0.0)          // positive
            {
                ve  = 0.5*( V(i,j)     + V(i+1,j)   )
                      // - 0.125 * Dxs(i+1)* Dxs(i+1)      / Dx(i)
                      // * ( ( V(i+1,j)   - V(i,j)     ) / Dxs(i+1)
                      //   - ( V(i,j)     - V(i-1,j)   ) / Dxs(i)    );
                      ;
            }
            else
            {
                ve  = 0.5*( V(i,j)     + V(i+1,j)   )
                      // - 0.125 * Dxs(i+1)* Dxs(i+1)      / Dx(i+1)
                      // * ( ( V(i+2,j)   - V(i+1,j)   ) / Dxs(i+2)
                      //   - ( V(i+1,j)   - V(i,j)     ) / Dxs(i+1)  );
                      ;
            }
            // vw
            // uwy = 0.5 * ( U[i-1][j] + U[i-1][j+1] );
            // if (uwy > 0.0)          // positive
            if (U(i,j) > 0.0)          // positive
            {
                vw  = 0.5*( V(i-1,j)   + V(i,j)     )
                      // - 0.125 * Dxs(i)  * Dxs(i)        / Dx(i-1)
                      // * ( ( V(i,j)     - V(i-1,j)   ) / Dxs(i)
                      //   - ( V(i-1,j)   - V(i-2,j)   ) / Dxs(i-1)  );
                      ;
            }
            else
            {
                vw  = 0.5*( V(i-1,j)   + V(i,j)     )
                      // - 0.125 * Dxs(i)  * Dxs(i)        / Dx(i)
                      // * ( ( V(i+1,j)   - V(i,j)     ) / Dxs(i+1)
                      //   - ( V(i,j)     - V(i-1,j)   ) / Dxs(i)    );
                      ;
            }

            // Old variables
            // uev
            // if (0.5*(V[i][j] + V[i+1][j]) >= 0.0)       // positive
            if (V(i,j) > 0.0)          // positive
            {
                uev = 0.5*( U(i,j)     + U(i,j+1)   )
                      // - 0.125 * Dys(j+1)* Dys(j+1)      / Dy(j)
                      // * ( ( U(i,j+1)   - U(i,j)     ) / Dys(j+1)
                      //   - ( U(i,j)     - U(i,j-1)   ) / Dys(j)    );
                      ;
            }
            else
            {
                uev = 0.5*( U(i,j)     + U(i,j+1)   )
                      // - 0.125 * Dys(j+1)* Dys(j+1)      / Dy(j+1)
                      // * ( ( U(i,j+2)   - U(i,j+1)   ) / Dys(j+2)
                      //   - ( U(i,j+1)   - U(i,j)     ) / Dys(j+1)  );
                      ;
            }
            // uwv
            // if (0.5*(V[i-1][j] + V[i][j]) >= 0.0)       // positive
            if (V(i,j) > 0.0)          // positive
            {
                uwv = 0.5*( U(i-1,j)   + U(i-1,j+1)   )
                      // - 0.125 * Dys(j+1)* Dys(j+1)      / Dy(j)
                      // * ( ( U(i-1,j+1) - U(i-1,j)   ) / Dys(j+1)
                      //   - ( U(i-1,j)   - U(i-1,j-1) ) / Dys(j)    );
                      ;
            }
            else
            {
                uwv = 0.5*( U(i-1,j)   + U(i-1,j+1) )
                      // - 0.125 * Dys(j+1)* Dys(j+1)      / Dy(j+1)
                      // * ( ( U(i-1,j+2) - U(i-1,j+1) ) / Dys(j+2)
                      //   - ( U(i-1,j+1) - U(i-1,j)   ) / Dys(j+1)  );
                      ;
            }

            FU(i,j)
            =   - ( ue  * ue  - uw  * uw  )     / Dxs(i+1)
                - ( un  * vnu - us  * vsu )     / Dy(j)
                + nu
                * (
                   (  ( U(i+1,j) - U(i,j)   ) / Dx(i+1)
                    - ( U(i,j)   - U(i-1,j) ) / Dx(i)     ) / Dxs(i+1)
                   +
                   (  ( U(i,j+1) - U(i,j)   ) / Dys(j+1)
                    - ( U(i,j)   - U(i,j-1) ) / Dys(j)    ) / Dy(j)
                  );

            FV(i,j)
            =   - ( vn  * vn  - vs  * vs  )     / Dys(j+1)
                - ( ve  * uev - vw  * uwv )     / Dx(i)
                + nu
                * (
                   (  ( V(i,j+1) - V(i,j)   ) / Dy(j+1)
                    - ( V(i,j)   - V(i,j-1) ) / Dy(j)     ) / Dys(j+1)
                   +
                   (  ( V(i+1,j) - V(i,j)   ) / Dxs(i+1)
                    - ( V(i,j)   - V(i-1,j) ) / Dxs(i)    ) / Dx(i)
                  );
        }
    }
}

void timeStepNSchemer()
{
    for (int i = 2; i <= nx-3; ++i)
    {
        for (int j = 2; j <= ny-3; ++j)
        {
            FU(i,j) = dt * FU(i,j);
            FV(i,j) = dt * FV(i,j);

            if (t == 1)         // Euler scheme for the 1st timestep
            {
                U(i,j) += FU(i,j);
            }
            else if (t == 2)    // Adams-Bashforth 2nd order for 2nd timestep
            {
                U(i,j) += 1.5*FU(i,j) - 0.5*FU1(i,j);
            }
            else                // Adams-Bashforth 3rd order subsequently
            {
                U(i,j) += (23.0*FU(i,j) - 16.0*FU1(i,j) + 5.0*FU2(i,j)) / 12.0;
            }

            if (t == 1)         // Euler scheme for the 1st timestep
            {
                V(i,j) += FV(i,j);
            }
            else if (t == 2)    // Adams-Bashforth 2nd order for 2nd timestep
            {
                V(i,j) += 1.5*FV(i,j) - 0.5*FV1(i,j);
            }
            else                // Adams-Bashforth 3rd order subsequently
            {
                V(i,j) += (23.0*FV(i,j) - 16.0*FV1(i,j) + 5.0*FV2(i,j)) / 12.0;
            }

            FU2(i,j) = FU1(i,j);
            FU1(i,j) = FU(i,j);
            FV2(i,j) = FV1(i,j);
            FV1(i,j) = FV(i,j);
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
        for (int i = 2; i <= nx-3; ++i)
        {
            for (int j = 2; j <= ny-3; ++j)
            {
                mChange =   ( U(i,j) - U(i-1,j) ) * Dy(j)
                          + ( V(i,j) - V(i,j-1) ) * Dx(i);

                if (abs(mChange) > mChangeMax)
                {
                  mChangeMax = abs(mChange);
                  mChangeMax_i = i;
                  mChangeMax_j = j;
                }

                pNew =
                    (
                     - P(i+1,j) * Dy(j) / Dxs(i+1) - P(i-1,j) * Dy(j) / Dxs(i)
                     - P(i,j+1) * Dx(i) / Dys(j+1) - P(i,j-1) * Dx(i) / Dys(j)
                     + rho * mChange / dt
                    )
                    /
                    (- Dy(j) / Dxs(i+1) - Dy(j) / Dxs(i)
                     - Dx(i) / Dys(j+1) - Dx(i) / Dys(j) );

                pChange = abs(pNew - P(i,j));
                if (pChange > pChangeMax)
                {
                    pChangeMax = pChange;
                    pChangeMax_i = i;
                    pChangeMax_j = j;
                }

                // Applying the SOR method through following formula:
                P(i,j) = P(i,j) + omega * (pNew - P(i,j));
                // P(i,j) = (1.0 - omega)*P(i,j) + omega * pNew; //same
            }
        }
        pressBoundary();
        ++h;
    }
    pIter = --h;
}

void velUpdater()
{
    for (int i = 2; i <= nx-3; ++i)
    {
        for (int j = 2; j <= ny-3; ++j)
        {
            U(i,j) = U(i,j) - (dt/rho) * ( P(i+1,j) - P(i,j) ) / Dxs(i+1);

            // uChangeMax += abs(Uo(i,j) - U(i,j));
            if (abs(Uo(i,j) - U(i,j)) > uChangeMax)
            {
                uChangeMax = abs(Uo(i,j) - U(i,j));
                uChangeMax_i = i;
                uChangeMax_j = j;
            }
            Uo(i,j) = U(i,j);

            V(i,j) = V(i,j) - (dt/rho) * ( P(i,j+1) - P(i,j) ) / Dys(j+1);

            // vChangeMax += abs(Vo(i,j) - V(i,j));
            if (abs(Vo(i,j) - V(i,j)) > vChangeMax)
            {
                vChangeMax = abs(Vo(i,j) - V(i,j));
                vChangeMax_i = i;
                vChangeMax_j = j;
            }
            Vo(i,j) = V(i,j);
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
    while   (t <= 10 || (t <= maxTimesteps
                         && (uChangeMax > uResidual || vChangeMax > vResidual)))
    {
        uChangeMax = 0.0;
        vChangeMax = 0.0;

        quick();
        timeStepNSchemer();
        velBoundary();
        pressure();
        velUpdater();
        velBoundary();
        progressPrinter();
        progressFiler();
        if (t % 50000 == 0)
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
