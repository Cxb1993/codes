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

int t               = 0;
int tPrConvrg       = -1;
int pIter           = -1;
double uChangeMax   = 1.0;
double vChangeMax   = 1.0;
//double wChangeMax   = 1.0;
double pChangeMax   = 1.0;
double mChangeMax   = 1.0;
double TF_X         = 0.0; //Total force in x direction
double TF_Y         = 0.0; //Total force in y direction
double CD_Temp      = 0.0; //drag coefficient
double CL_Temp      = 0.0; //lift coefficient

double scriptRunningTime = 0.0;

Matrix<double,2> Eta(nx,ny);
Matrix<double,2> U(nx,ny);
Matrix<double,2> V(nx,ny);
//Matrix<double,2> W(nx,ny);
Matrix<double,2> Uo(nx,ny);
Matrix<double,2> Vo(nx,ny);
//Matrix<double,2> Wo(nx,ny);
Matrix<double,2> FU(nx,ny);
Matrix<double,2> FV(nx,ny);
//Matrix<double,2> FW(nx,ny);
Matrix<double,2> FU1(nx,ny);
Matrix<double,2> FV1(nx,ny);
//Matrix<double,3> FW1(nx,ny,nz);
Matrix<double,2> FU2(nx,ny);
Matrix<double,2> FV2(nx,ny);
//atrix<double,3> FW2(nx,ny,nz);
Matrix<double,2> P(nx,ny);
Matrix<double,2> MC(nx,ny);
Matrix<double,2> PC(nx,ny);
Matrix<double,2> U_DS(nx,ny); // U double stars
Matrix<double,2> V_DS(nx,ny); // V double stars
Matrix<double,2> VF_X(nx,ny); //virtual force x direction
Matrix<double,2> VF_Y(nx,ny); // virtual force y direction
Matrix<double,1> F_X(ny);
Matrix<double,1> F_Y(nx);
Matrix<double,2> STR(nx,ny); //streamline
Matrix<double,2> VORT(nx,ny); //vorticity

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
        U(1,j) = 1.0;
        U(0,j) = U(1,j);
        V(1,j) = 0.0 ;
        V(0,j) = V(1,j);

        // Right boundary
        U(nx-3,j) = U(nx-4,j);
        U(nx-2,j) = U(nx-3,j);
        V(nx-2,j) = V(nx-3,j);
        V(nx-1,j) = V(nx-2,j);
    }
    for (int i = 0; i <= nx-1; ++i)
    {
        // Bottom boundary
        U(i,1) = 1.0;
        U(i,0) = U(i,1);
        V(i,1) = 0.0;
        V(i,0) = V(i,1);

        // Top boundary
        U(i,ny-2) = 1.0;
        U(i,ny-1) = U(i,ny-2);
        V(i,ny-3) = 0.0;
        V(i,ny-2) = V(i,ny-3);
    }

}

void pressBoundary()
{
    // Applying the pressure boundary conditions

    for (int j = 1; j <= ny-2; ++j)
    {
        // Left boundary
        P(1,j) = 0.0;                // Neumann
       // P(0,j) = pIn;                 // Dirichlet
         // P(0,j) = 0.0;

        // Right boundary
        P(nx-2,j) = P(nx-3,j);          // Neumann
        // P(nx-1,j) = pOut;               // Dirichlet
        //P(nx-1,j) = 0.0;
    }
    for (int i = 1; i <= nx-2; ++i)
    {
        // Top boundary
        P(i,ny-2) = 0.0;          // Neumann
        // P(0,ny-1) = pInitial;           // fixed pressure point

        // Bottom boundary
        P(i,1) = 0.0;                // Neumann
    }
}

void detEta()
{
//determine Eta

    #pragma omp parallel for\
         default(none)\
         shared(X, Y, Dxs, Dys, Eta)


    for (int j = 2; j < ny-3; ++j)
    {
        for (int i = 2; i < nx-3; ++i)
        {
    double c;                   // distance between center cylinder and position of eta
    double dc;
    double dsl ;
    double cc  ;
    int     n  = 20;            //number of sub grid
    double sx [n], sy [n];
    int   sEta = 0;    // sub-eta
    //double  Xc = 6.5 ; // position of center cylinder x-direction
    //double  Yc = 3.5 ; // position of center cylinder y-direction
           c    = sqrt( ( (((X(i+1)+X(i)) / 2)- Xc) * (((X(i+1)+X(i)) / 2)- Xc))
                + (( ((Y(j+1)+Y(j)) / 2)- Yc) * (((Y(j+1)+Y(j)) / 2) - Yc)) );

           dsl  = sqrt( ( Dxs(i) * (Dxs(i))) + (( Dys(j)) * (Dys(j))) );
           dc   = fabs(r - c);
           if  ( c <= r && dc > dsl )
           {
               Eta(i,j) = 1;
           }
           else if ( c > r && dc > dsl )
           {
               Eta(i,j) = 0;
           }
           else
           {
               sx[0] = ((X(i+1)+X(i)) / 2)  - (Dxs(i) * 0.5) + (0.5 * (Dxs(i)/n)) ;
               sy[0] = ((Y(j+1)+Y(j)) / 2)  - (Dys(j) * 0.5) + (0.5 * (Dys(j)/n));
               for (int k = 1; k < n; ++k)
               {
                   sx [k] = sx [k-1] + (Dxs(i) / n);
                   sy [k] = sy [k-1] + (Dys(j) / n);

                for (int l = 0; l < n; ++l)
                {
                    for (int m = 0; m < n; ++m)
                    {
                        cc    = sqrt( (( sx[l]-Xc) * (sx[l]-Xc)) + (( sy[m]-Yc) * (sy[m]-Yc)) );
                        if  ( cc <= r  )
                        {
                            sEta = sEta + 1 ;
                        }
                        else
                        {
                            sEta = sEta + 0;
                        }
                    }
                }
                Eta(i,j) = sEta / (n*n*1.0);
                sEta = 0;
            }
            }
        }
    }
}
void momentum()
{
    #pragma omp parallel for\
        default(none)\
        shared(FU, FV, U, V, Dx, Dxs, Dy, Dys)
    for (int i = 2; i <= nx-3; ++i)
    {
        for (int j = 2; j <= ny-3; ++j)
        {
                double
                    ue, uw, un, us, vnu, vsu,
                    ve, vw, vn, vs, uev, uwv;


                ue = uw = un = us = vnu = vsu =
                ve = vw = vn = vs = uev = uwv = 0.0;
                // Call the selected velocity scheme function
                quick(i, j,\
                    ue, uw, un, us, vnu, vsu,\
                    ve, vw, vn, vs, uev, uwv);

                FU(i,j)
                =   - ( ue  * ue  - uw  * uw  )   / Dxs(i)
                    - ( un  * vnu - us  * vsu )   / Dy(j)
                    + nu
                    * (
                       (  ( U(i+1,j) - U(i,j)   ) / Dx(i+1)
                        - ( U(i,j)   - U(i-1,j) ) / Dx(i)    ) / Dxs(i)
                       +
                       (  ( U(i,j+1) - U(i,j)   ) / Dys(j)
                        - ( U(i,j)   - U(i,j-1) ) / Dys(j-1) ) / Dy(j)
                       );

                FV(i,j)
                =   - ( uev * ve  - uwv * vw  )   / Dx(i)
                    - ( vn  * vn  - vs  * vs  )   / Dys(j)
                    + nu
                    * (
                       (  ( V(i+1,j) - V(i,j)   ) / Dxs(i)
                        - ( V(i,j)   - V(i-1,j) ) / Dxs(i-1) ) / Dx(i)
                       +
                       (  ( V(i,j+1) - V(i,j)   ) / Dy(j+1)
                        - ( V(i,j)   - V(i,j-1) ) / Dy(j)    ) / Dys(j)
                        );

        }
    }/*-- End of omp parallel for --*/
}

void timeStepNSchemer()
{
    const int tTemp = t;
    #pragma omp parallel for\
        default(none)\
        shared(U, V, FU, FU1, FU2, FV, FV1, FV2)
    for (int i = 2; i <= nx-3; ++i)
    {
        for (int j = 2; j <= ny-3; ++j)
        {

                FU(i,j) = dt * FU(i,j);
                FV(i,j) = dt * FV(i,j);


                if (tTemp == 1)         // Euler scheme
                {
                    U(i,j) += FU(i,j);
                    V(i,j) += FV(i,j);

                }
                else if (tTemp == 2)    // Adams-Bashforth 2nd order
                {
                    U(i,j) += 1.5 * FU(i,j) - 0.5 * FU1(i,j);
                    V(i,j) += 1.5 * FV(i,j) - 0.5 * FV1(i,j);

                }
                else                // Adams-Bashforth 3rd order
                {
                    U(i,j) += ( 23.0 * FU(i,j) - 16.0 * FU1(i,j)
                                  + 5.0 * FU2(i,j) ) / 12.0;

                    V(i,j) += ( 23.0 * FV(i,j) - 16.0 * FV1(i,j)
                                  + 5.0 * FV2(i,j) ) / 12.0;

                }

                FU2(i,j) = FU1(i,j);
                FU1(i,j) = FU(i,j);
                FV2(i,j) = FV1(i,j);
                FV1(i,j) = FV(i,j);

        }
    }/*-- End of omp parallel for --*/
}

void pressure()
{
    int h = 0;
    double pNew = 0.0;      // calculation of pressure for the next timestep
    double pChange = 0.0;   // pressure residual at a particular cell
    double mChange = 0.0;   // mass residual at a particular cell
    pChangeMax = 1.0;

    while (h <= maxPrIters && pChangeMax > pResidual)
    {
        pChangeMax = 0.0;   // maximum pr. residual in current iteration
        mChangeMax = 0.0;   // maximum mass residual in current iteration

        #pragma omp parallel for\
            default(none)\
            shared(P, U, V, Dx, Dxs, Dy, Dys)\
            private(mChange, pChange, pNew)\
            reduction(max:mChangeMax, pChangeMax)
        for (int i = 2; i <= nx-3; ++i)
        {
            for (int j = 2; j <= ny-3; ++j)
            {
                    mChange =   ( U(i,j) - U(i-1,j) ) * Dy(j)
                              + ( V(i,j) - V(i,j-1) ) * Dx(i);

                    pNew =
                        (
                         - P(i+1,j) * Dy(j) / Dxs(i)
                         - P(i-1,j) * Dy(j) / Dxs(i-1)
                         - P(i,j+1) * Dx(i) / Dys(j)
                         - P(i,j-1) * Dx(i) / Dys(j-1)
                         + rho * mChange / dt
                        )
                        /
                        (- Dy(j) / Dxs(i) - Dy(j) / Dxs(i-1)
                         - Dx(i) / Dys(j) - Dx(i) / Dys(j-1));

                    pChange = std::abs(pNew - P(i,j));

                    // Applying the SOR method through following formula:
                    P(i,j) += omega * (pNew - P(i,j));
                    // P(i,j) = (1.0 - omega)*P(i,j) + omega * pNew; //same

                    if (std::abs(mChange) > mChangeMax)
                    {
                        mChangeMax = std::abs(mChange);
                    }
                    if (pChange > pChangeMax)
                    {
                        pChangeMax = pChange;
                    }

            }
        }/*-- End of omp parallel for --*/
        pressBoundary();
        ++h;
    }
    pIter = --h;
    if (pIter < maxPrIters && tPrConvrg == -1)
    {
        tPrConvrg = t;
    }
}

void velUpdater()
{
    #pragma omp parallel for\
        default(none)\
        shared(U, Uo, V, Vo, P, Dxs, Dys)\

    for (int i = 2; i <= nx-3; ++i)
    {
        for (int j = 2; j <= ny-3; ++j)
        {
                U(i,j) = U(i,j)
                           - (dt/rho) * ( P(i+1,j) - P(i,j) ) / Dxs(i);



                V(i,j) = V(i,j)
                           - (dt/rho) * ( P(i,j+1) - P(i,j) ) / Dys(j);


        }
    }/*-- End of omp parallel for --*/
    // uChangeMax = uChangeMax / static_cast<double>(nx * ny * nz);
    // vChangeMax = vChangeMax / static_cast<double>(nx * ny * nz);
    // wChangeMax = wChangeMax / static_cast<double>(nx * ny * nz);
}

void velUpdater_Eta()
{

    #pragma omp parallel for\
        default(none)\
        shared(U, Uo, V, Vo, U_DS, V_DS, Eta, Dxs, Dys)\
        reduction(max:uChangeMax, vChangeMax)
    for (int i = 2; i <= nx-3; ++i)
    {
        for (int j = 2; j <= ny-3; ++j)
        {
            double Us = 0;
            U_DS(i,j)  = U(i,j);
            V_DS(i,j)  = V(i,j) ;

            if ( Eta (i,j) == 1.0 )
            {
                U(i,j) = 0.0;
                V(i,j) = 0.0;
            }

            else if ( Eta(i,j) == 0.0 )
            {
                U(i,j) = U_DS(i,j) ;
                V(i,j) = V_DS(i,j) ;
            }

            else
            {
                U(i,j)     = (Eta (i,j) * Us)  +( (1 - Eta(i,j) ) *  U_DS(i,j) );
                V(i,j)     = (Eta (i,j) * Us)  +( (1 - Eta(i,j) ) *  V_DS(i,j) );
            }


            // uChangeMax += abs(Uo(i,j,k) - U(i,j,k));
                if (std::abs(Uo(i,j) - U(i,j)) > uChangeMax)
                {
                    uChangeMax = std::abs(Uo(i,j) - U(i,j));
                }
                Uo(i,j) = U(i,j);

            // vChangeMax += std::abs(Vo(i,j,k) - V(i,j,k));
                if (std::abs(Vo(i,j) - V(i,j)) > vChangeMax)
                {
                    vChangeMax = std::abs(Vo(i,j) - V(i,j));
                }
                Vo(i,j) = V(i,j);
        }
    }
    // uChangeMax = uChangeMax / static_cast<double>(nx * ny);
    // vChangeMax = vChangeMax / static_cast<double>(nx * ny);
}

void virtualforce()
{

    #pragma omp parallel for\
        default(none)\
        shared(U, V, U_DS, V_DS, Eta, VF_X, VF_Y)

    for (int i = 2; i <= nx-3; ++i)
    {

        for (int j = 2; j <= ny-3; ++j)
        {

              if ( Eta (i,j) == 0.0 )
            {
                VF_X(i,j) = 0.0;
                VF_Y(i,j) = 0.0;
            }

            else
            {
                VF_X(i,j)  = Eta (i,j) * (U(i,j) - U_DS(i,j)) / dt ;
                VF_Y(i,j)  = Eta (i,j) * (V(i,j) - V_DS(i,j)) / dt ;
            }

        }
    }

    // calculate draq coefficient using simpson's 1/3 rule

    for (int j = 2; j <= ny-3; ++j)
    {

        F_X(j)    = 0.0;
        for (int i = 2; i <= nx-3; ++i)
        {
            //F_X(j)      = F_X(j) + (VF_X(i,j) * 0.066);
            F_X(j)  = F_X(j) + (((VF_X(i-1,j) * Dxs(i-1) ) + (4.0 * VF_X(i,j) *
                      ((Dxs(i+1)+ Dxs(i-1))/2.0)) + (VF_X(i+1,j) * Dxs(i+1))) / 6.0)  ;
        }
    }

    TF_X = 0.0;
    for (int j = 2; j <= ny-3; ++j)
    {
        //TF_X    = TF_X + (F_X(j) * 0.066);
        TF_X  = TF_X + (((F_X(j-1) * Dys(j-1) ) + (4.0 * F_X(j) *
                ((Dys(j+1)+ Dys(j-1))/2.0)) + (F_X(j+1) * Dys(j+1))) / 6.0)  ;
    }

    CD_Temp = (-2.0) * TF_X;
    //CD      = CD + CD_Temp;

    // calculate lift coefficient using simpson's 1/3 rule

    for (int i = 1; i <= nx-3; ++i)
    {
        F_Y(i)  = 0.0;
        for (int j = 2; j <= ny-3; ++j)
        {
            F_Y(i)  = F_Y(i) + (((VF_Y(i,j-1) * Dy(j-1) ) + (4.0 * VF_Y(i,j) *
                      ((Dy(j+1)+ Dy(j-1))/2.0)) + (VF_Y(i,j+1) * Dy(j+1))) / 6.0)  ;
        }
    }

    TF_Y = 0.0;
    for (int n = 2; n <= nx-3; ++n)
    {
        TF_Y  = TF_Y + (((F_Y(n-1) * Dx(n-1) ) + (4.0 * F_Y(n) *
                ((Dx(n+1)+ Dx(n-1))/2.0)) + (F_Y(n+1) * Dx(n+1))) / 6.0)  ;
    }

    CL_Temp = (-2.0) * TF_Y;
    //CL      = CL + CL_Temp;
    // uChangeMax = uChangeMax / static_cast<double>(nx * ny);
    // vChangeMax = vChangeMax / static_cast<double>(nx * ny);
}

void STREAM()
// calculate stream funcion and vorticity
{
    int i = 1;
    int j = 1;
    STR (i,j)= 0.0;
    VORT(i,j)= 0.0;
    #pragma omp parallel for\
        default(none)\
        shared(U, V, X, Y, STR, VORT, Eta, VF_X, VF_Y)\

    for (int i = 2; i <= nx-3; ++i)
    {
        STR (i,1)= 0.0;
        VORT(i,1)= 0.0;

        for (int j = 2; j <= ny-3; ++j)
        {
          STR (i,j) = STR (i,j-1) + (Y (j+1) - Y(j)) * U(i,j) ;
          VORT(i,j) = -(U(i-1,j) - U(i-1,j-1)) / (0.5 * (Y(j+1)-Y(j-1))) + (V(i,j-1)
                      - V(i-1,j-1)) / (0.5 * (X(i+1)-X(i-1)));
        }
    }

}


void mainSolver()
{
    std::chrono::steady_clock::time_point tStart
        = std::chrono::steady_clock::now();
    initial();
    velBoundary();
    pressBoundary();
    filerCreateCD();
    filerCreateCL();
    // filerAllSol(t);
    // filerCreateProgress();

    // End of timestep

    #ifdef _OPENMP
    omp_set_num_threads(numOfThreads);
    #endif

    t = 1;
    while   (t <= 10 || (t <= maxTimesteps
                         && (   uChangeMax > uResidual
                             || vChangeMax > vResidual )))
    {
        uChangeMax = 0.0;
        vChangeMax = 0.0;


        momentum();
        timeStepNSchemer();
        velBoundary();
        detEta();
        pressure();
        velUpdater();
        velUpdater_Eta();
        virtualforce();
        STREAM();
        velBoundary();
        printerProgress();


         if (t % 20000 == 0)
         {
             filerAllSol(t);
             // file intermittent solution
         }
         if (t % 2000 == 0)
         {
             filerCD();
             filerCL();
             // file intermittent solution
         }

        // End of timestep

        ++t;
    }
    --t;
    std::chrono::steady_clock::time_point tEnd
        = std::chrono::steady_clock::now();
    scriptRunningTime =
        std::chrono::duration_cast<std::chrono::minutes>(tEnd - tStart).count();
}
