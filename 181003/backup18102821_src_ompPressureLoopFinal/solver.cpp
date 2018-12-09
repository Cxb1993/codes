/*
 * solver.cpp
 *
 *     Project: Finite Volume Navier-Stokes Solver
 *      Author: Syed Ahmad Raza
 */

#include <cmath>        // math functions
#include <ctime>        // to time the script
#include <chrono>       // to measure and display the elapsed time
#include <omp.h>        // openMP header
#include <iostream>     // functions for input and output to console

#include "constants.h"
#include "filers.h"
#include "printers.h"
#include "velschemers.h"
#include "solver.h"

using namespace Numeric_lib;

double
    ue, uw, un, us, uf, ub, vnu, vsu, wfu, wbu,
    ve, vw, vn, vs, vf, vb, uev, uwv, wfv, wbv,
    we, ww, wn, ws, wf, wb, uew, uww, vnw, vsw;

int t               = 0;
int pIter           = -1;
double uChangeMax   = 1.0;
double vChangeMax   = 1.0;
double wChangeMax   = 1.0;
double pChangeMax   = 1.0;
double mChangeMax   = 1.0;
int mChangeMax_i    = -1;
int mChangeMax_j    = -1;
int mChangeMax_k    = -1;
int pChangeMax_i    = -1;
int pChangeMax_j    = -1;
int pChangeMax_k    = -1;
int uChangeMax_i    = -1;
int uChangeMax_j    = -1;
int uChangeMax_k    = -1;
int vChangeMax_i    = -1;
int vChangeMax_j    = -1;
int vChangeMax_k    = -1;
int wChangeMax_i    = -1;
int wChangeMax_j    = -1;
int wChangeMax_k    = -1;

double scriptRunningTime = 0.0;

Matrix<double,3> U(nx,ny,nz);
Matrix<double,3> V(nx,ny,nz);
Matrix<double,3> W(nx,ny,nz);
Matrix<double,3> Uo(nx,ny,nz);
Matrix<double,3> Vo(nx,ny,nz);
Matrix<double,3> Wo(nx,ny,nz);
Matrix<double,3> FU(nx,ny,nz);
Matrix<double,3> FV(nx,ny,nz);
Matrix<double,3> FW(nx,ny,nz);
Matrix<double,3> FU1(nx,ny,nz);
Matrix<double,3> FV1(nx,ny,nz);
Matrix<double,3> FW1(nx,ny,nz);
Matrix<double,3> FU2(nx,ny,nz);
Matrix<double,3> FV2(nx,ny,nz);
Matrix<double,3> FW2(nx,ny,nz);
Matrix<double,3> P(nx,ny,nz);
Matrix<double,3> MC(nx,ny,nz);
Matrix<double,3> PC(nx,ny,nz);

int numOfThreads = 0;

void initial()
{
    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < ny; ++j)
        {
            for (int k = 0; k < nz; ++k)
            {
                U(i,j,k) = uInitial;
                V(i,j,k) = vInitial;
                W(i,j,k) = wInitial;
                P(i,j,k) = pInitial;
            }
        }
    }
}

void velBoundary()
{
    // Applying the velocity boundary conditions

    for (int i = 0; i <= nx-1; ++i)
    {
        for (int j = 0; j <= ny-1; ++j)
        {
            // Back boundary
            U(i,j,1) = -U(i,j,2);
            U(i,j,0) = U(i,j,1);
            V(i,j,1) = -V(i,j,2);
            V(i,j,0) = V(i,j,1);
            W(i,j,1) = 0.0;
            W(i,j,0) = W(i,j,1);

            // Front boundary
            U(i,j,nz-2) = -U(i,j,nz-3);
            U(i,j,nz-1) = U(i,j,nz-2);
            V(i,j,nz-2) = -V(i,j,nz-3);
            V(i,j,nz-1) = V(i,j,nz-2);
            W(i,j,nz-3) = 0.0;
            W(i,j,nz-2) = W(i,j,nz-3);
        }
    }
    for (int j = 0; j <= ny-1; ++j)
    {
        for (int k = 0; k <= nz-1; ++k)
        {
            // West boundary
            U(1,j,k) = 0.0;
            U(0,j,k) = U(1,j,k);
            V(1,j,k) = -V(2,j,k);
            V(0,j,k) = V(1,j,k);
            W(1,j,k) = -W(2,j,k);
            W(0,j,k) = W(1,j,k);

            // East boundary
            U(nx-3,j,k) = 0.0;
            U(nx-2,j,k) = U(nx-3,j,k);
            V(nx-2,j,k) = -V(nx-3,j,k);
            V(nx-1,j,k) = V(nx-2,j,k);
            W(nx-2,j,k) = -W(nx-3,j,k);
            W(nx-1,j,k) = W(nx-2,j,k);
        }
    }
    for (int i = 0; i <= nx-1; ++i)
    {
        for (int k = 0; k <= nz-1; ++k)
        {
            // South boundary
            U(i,1,k) = -U(i,2,k);
            U(i,0,k) = U(i,1,k);
            V(i,1,k) = 0.0;
            V(i,0,k) = V(i,1,k);
            W(i,1,k) = -W(i,2,k);
            W(i,0,k) = W(i,1,k);

            // North boundary
            U(i,ny-2,k) = -U(i,ny-3,k) + 2*uIn;
            U(i,ny-1,k) = U(i,ny-2,k);
            V(i,ny-3,k) = 0.0;
            V(i,ny-2,k) = V(i,ny-3,k);
            W(i,ny-2,k) = -W(i,ny-3,k);
            W(i,ny-1,k) = W(i,ny-2,k);
        }
    }
}

void pressBoundary()
{
    // Applying the pressure boundary conditions

    for (int i = 1; i <= nx-2; ++i)
    {
        for (int j = 1; j <= ny-2; ++j)
        {
            // Back boundary
            P(i,j,1) = P(i,j,2);                // Neumann

            // Front boundary
            P(i,j,nz-2) = P(i,j,nz-3);          // Neumann
        }
    }
    for (int j = 1; j <= ny-2; ++j)
    {
        for (int k = 1; k <= nz-2; ++k)
        {
            // West boundary
            P(1,j,k) = P(2,j,k);                // Neumann

            // East boundary
            P(nx-2,j,k) = P(nx-3,j,k);          // Neumann
        }
    }
    for (int i = 1; i <= nx-2; ++i)
    {
        for (int k = 1; k <= nz-2; ++k)
        {
            // South boundary
            P(i,1,k) = P(i,2,k);                // Neumann

            // North boundary
            P(i,ny-2,k) = P(i,ny-3,k);          // Neumann
        }
    }
}

void momentum(void (*velschemer)(int i, int j, int k))
{
    ue = uw = un = us = uf = ub = vnu = vsu = wfu = wbu =
    ve = vw = vn = vs = vf = vb = uev = uwv = wfv = wbv =
    we = ww = wn = ws = wf = wb = uew = uww = vnw = vsw = 0.0;

    for (int i = 2; i <= nx-3; ++i)
    {
        for (int j = 2; j <= ny-3; ++j)
        {
            for (int k = 2; k <= nz-3; ++k)
            {
                // Call the selected velocity scheme function
                (*velschemer)(i,j,k);

                FU(i,j,k)
                =   - ( ue  * ue  - uw  * uw  )   / Dxs(i)
                    - ( un  * vnu - us  * vsu )   / Dy(j)
                    - ( uf  * wfu - ub  * wbu )   / Dz(k)
                    + nu
                    * (
                       (  ( U(i+1,j,k) - U(i,j,k)   ) / Dx(i+1)
                        - ( U(i,j,k)   - U(i-1,j,k) ) / Dx(i)    ) / Dxs(i)
                       +
                       (  ( U(i,j+1,k) - U(i,j,k)   ) / Dys(j)
                        - ( U(i,j,k)   - U(i,j-1,k) ) / Dys(j-1) ) / Dy(j)
                       +
                       (  ( U(i,j,k+1) - U(i,j,k)   ) / Dzs(k)
                        - ( U(i,j,k)   - U(i,j,k-1) ) / Dzs(k-1) ) / Dz(k)
                      );

                FV(i,j,k)
                =   - ( uev * ve  - uwv * vw  )   / Dx(i)
                    - ( vn  * vn  - vs  * vs  )   / Dys(j)
                    - ( vf  * wfv - vb  * wbv )   / Dz(k)
                    + nu
                    * (
                       (  ( V(i+1,j,k) - V(i,j,k)   ) / Dxs(i)
                        - ( V(i,j,k)   - V(i-1,j,k) ) / Dxs(i-1) ) / Dx(i)
                       +
                       (  ( V(i,j+1,k) - V(i,j,k)   ) / Dy(j+1)
                        - ( V(i,j,k)   - V(i,j-1,k) ) / Dy(j)    ) / Dys(j)
                       +
                       (  ( V(i,j,k+1) - V(i,j,k)   ) / Dzs(k)
                        - ( V(i,j,k)   - V(i,j,k-1) ) / Dzs(k-1) ) / Dz(k)
                      );

                FW(i,j,k)
                =   - ( uew * we  - uww * ww  )   / Dx(i)
                    - ( vnw * wn  - vsw * ws  )   / Dy(j)
                    - ( wf  * wf  - wb  * wb  )   / Dzs(k)
                    + nu
                    * (
                       (  ( W(i+1,j,k) - W(i,j,k)   ) / Dxs(i)
                        - ( W(i,j,k)   - W(i-1,j,k) ) / Dxs(i-1) ) / Dx(i)
                       +
                       (  ( W(i,j+1,k) - W(i,j,k)   ) / Dys(j)
                        - ( W(i,j,k)   - W(i,j-1,k) ) / Dys(j-1) ) / Dy(j)
                       +
                       (  ( W(i,j,k+1) - W(i,j,k)   ) / Dz(k+1)
                        - ( W(i,j,k)   - W(i,j,k-1) ) / Dz(k)    ) / Dzs(k)
                      );
            }
        }
    }
}

void timeStepNSchemer()
{
    for (int i = 2; i <= nx-3; ++i)
    {
        for (int j = 2; j <= ny-3; ++j)
        {
            for (int k = 2; k <= nz-3; ++k)
            {
                FU(i,j,k) = dt * FU(i,j,k);
                FV(i,j,k) = dt * FV(i,j,k);
                FW(i,j,k) = dt * FW(i,j,k);

                if (t == 1)         // Euler scheme
                {
                    U(i,j,k) += FU(i,j,k);
                    V(i,j,k) += FV(i,j,k);
                    W(i,j,k) += FW(i,j,k);
                }
                else if (t == 2)    // Adams-Bashforth 2nd order
                {
                    U(i,j,k) += 1.5 * FU(i,j,k) - 0.5 * FU1(i,j,k);
                    V(i,j,k) += 1.5 * FV(i,j,k) - 0.5 * FV1(i,j,k);
                    W(i,j,k) += 1.5 * FW(i,j,k) - 0.5 * FW1(i,j,k);
                }
                else                // Adams-Bashforth 3rd order
                {
                    U(i,j,k) += ( 23.0 * FU(i,j,k) - 16.0 * FU1(i,j,k)
                                  + 5.0 * FU2(i,j,k) ) / 12.0;

                    V(i,j,k) += ( 23.0 * FV(i,j,k) - 16.0 * FV1(i,j,k)
                                  + 5.0 * FV2(i,j,k) ) / 12.0;

                    W(i,j,k) += ( 23.0 * FW(i,j,k) - 16.0 * FW1(i,j,k)
                                  + 5.0 * FW2(i,j,k) ) / 12.0;
                }

                FU2(i,j,k) = FU1(i,j,k);
                FU1(i,j,k) = FU(i,j,k);
                FV2(i,j,k) = FV1(i,j,k);
                FV1(i,j,k) = FV(i,j,k);
                FW2(i,j,k) = FW1(i,j,k);
                FW1(i,j,k) = FW(i,j,k);
            }
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

        #pragma omp parallel for default(none)\
            shared(P, U, V, W, Dx, Dxs, Dy, Dys, Dz, Dzs, numOfThreads,\
                mChangeMax, mChangeMax_i, mChangeMax_j, mChangeMax_k,\
                pChangeMax, pChangeMax_i, pChangeMax_j, pChangeMax_k)\
            private(mChange, pChange, pNew)
        for (int i = 2; i <= nx-3; ++i)
        {
            for (int j = 2; j <= ny-3; ++j)
            {
                for (int k = 2; k <= nz-3; ++k)
                {
                    mChange =   ( U(i,j,k) - U(i-1,j,k) ) * Dy(j) * Dz(k)
                              + ( V(i,j,k) - V(i,j-1,k) ) * Dx(i) * Dz(k)
                              + ( W(i,j,k) - W(i,j,k-1) ) * Dx(i) * Dy(j);

                    if (std::abs(mChange) > mChangeMax)
                    {
                      mChangeMax = std::abs(mChange);
                      mChangeMax_i = i;
                      mChangeMax_j = j;
                      mChangeMax_k = k;
                    }

                    pNew =
                        (
                         - P(i+1,j,k) * Dy(j) * Dz(k) / Dxs(i)
                         - P(i-1,j,k) * Dy(j) * Dz(k) / Dxs(i-1)
                         - P(i,j+1,k) * Dx(i) * Dz(k) / Dys(j)
                         - P(i,j-1,k) * Dx(i) * Dz(k) / Dys(j-1)
                         - P(i,j,k+1) * Dx(i) * Dy(j) / Dzs(k)
                         - P(i,j,k-1) * Dx(i) * Dy(j) / Dzs(k-1)
                         + rho * mChange / dt
                        )
                        /
                        (- Dy(j) * Dz(k) / Dxs(i) - Dy(j) * Dz(k) / Dxs(i-1)
                         - Dx(i) * Dz(k) / Dys(j) - Dx(i) * Dz(k) / Dys(j-1)
                         - Dx(i) * Dy(j) / Dzs(k) - Dx(i) * Dy(j) / Dzs(k-1) );

                    pChange = std::abs(pNew - P(i,j,k));
                    if (pChange > pChangeMax)
                    {
                        pChangeMax = pChange;
                        pChangeMax_i = i;
                        pChangeMax_j = j;
                        pChangeMax_k = k;
                    }

                    // Applying the SOR method through following formula:
                    P(i,j,k) += omega * (pNew - P(i,j,k));
                    // P(i,j) = (1.0 - omega)*P(i,j) + omega * pNew; //same

                    #ifdef _OPENMP
                    numOfThreads = omp_get_num_threads();
                    #endif
                }
            }
        }/*-- End of omp parallel for --*/
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
            for (int k = 2; k <= nz-3; ++k)
            {
                U(i,j,k) = U(i,j,k)
                           - (dt/rho) * ( P(i+1,j,k) - P(i,j,k) ) / Dxs(i);

                // uChangeMax += abs(Uo(i,j,k) - U(i,j,k));
                if (std::abs(Uo(i,j,k) - U(i,j,k)) > uChangeMax)
                {
                    uChangeMax = std::abs(Uo(i,j,k) - U(i,j,k));
                    uChangeMax_i = i;
                    uChangeMax_j = j;
                    uChangeMax_k = k;
                }
                Uo(i,j,k) = U(i,j,k);

                V(i,j,k) = V(i,j,k)
                           - (dt/rho) * ( P(i,j+1,k) - P(i,j,k) ) / Dys(j);

                // vChangeMax += std::abs(Vo(i,j,k) - V(i,j,k));
                if (std::abs(Vo(i,j,k) - V(i,j,k)) > vChangeMax)
                {
                    vChangeMax = std::abs(Vo(i,j,k) - V(i,j,k));
                    vChangeMax_i = i;
                    vChangeMax_j = j;
                    vChangeMax_k = k;
                }
                Vo(i,j,k) = V(i,j,k);

                W(i,j,k) = W(i,j,k)
                           - (dt/rho) * ( P(i,j,k+1) - P(i,j,k) ) / Dzs(k);

                // wChangeMax += std::abs(Wo(i,j,k) - W(i,j,k));
                if (std::abs(Wo(i,j,k) - W(i,j,k)) > wChangeMax)
                {
                    wChangeMax = std::abs(Wo(i,j,k) - W(i,j,k));
                    wChangeMax_i = i;
                    wChangeMax_j = j;
                    wChangeMax_k = k;
                }
                Wo(i,j,k) = W(i,j,k);
            }
        }
    }
    // uChangeMax = uChangeMax / static_cast<double>(nx * ny * nz);
    // vChangeMax = vChangeMax / static_cast<double>(nx * ny * nz);
    // wChangeMax = wChangeMax / static_cast<double>(nx * ny * nz);
}

void mainSolver()
{
    std::chrono::steady_clock::time_point tStart
        = std::chrono::steady_clock::now();
    initial();
    velBoundary();
    pressBoundary();
    // filerAllSol(t);
    // filerCreateProgress();

    // End of timestep

    #ifdef _OPENMP
    omp_set_num_threads(20);
    #endif

    t = 1;
    while   (t <= 10 || (t <= maxTimesteps
                         && (   uChangeMax > uResidual
                             || vChangeMax > vResidual
                             || wChangeMax > wResidual  )))
    {
        uChangeMax = 0.0;
        vChangeMax = 0.0;
        wChangeMax = 0.0;

        if (velScheme == "upwind")
            momentum(&upwind);
        else if (velScheme == "quick")
            momentum(&quick);
        timeStepNSchemer();
        velBoundary();
        pressure();
        velUpdater();
        velBoundary();
        printerProgress();
        // filerProgress();
        // if (t % 30000 == 0)
        // {
        //     filerAllSol(t);     // file intermittent solution
        // }
        #ifdef _OPENMP
        std::cout << numOfThreads << std::endl;
        #endif

        // End of timestep

        ++t;
    }
    --t;
    std::chrono::steady_clock::time_point tEnd
        = std::chrono::steady_clock::now();
    scriptRunningTime =
        std::chrono::duration_cast<std::chrono::minutes>(tEnd - tStart).count();
}
