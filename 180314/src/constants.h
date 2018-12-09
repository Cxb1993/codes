/*
 * constants.h
 *
 *  Created on: 2017-09-27
 *      Author: Syed Ahmad Raza
 */

#ifndef CONSANTS_H
#define CONSTANTS_H

#include <string>       // for sring operations

const double pi(3.14159265);

// Domain, range
const int x_y_ratio = 1;
const double ly = 1.0;
const double lx = ly * static_cast<double>(x_y_ratio);

// Dimensions
const double D = ly;
const double L = lx;

// Number of nodes on the x and y axes
// Compile-time constants needed for arrays
const int ny = 100;
const int nx = ny * x_y_ratio;

// Intial and/or boundary conditions
// const double uin = Uin;
const double uIn = 1.0;
const double vIn = 0.0;
// const double pIn = 1.00016e6;
// const double pOut = 1.0e6;
    // specify pIn and pOut correctly for Dirichlet pressure conditions
const double uInitial = 0.0;
const double vInitial = 0.0;
const double pInitial = 0.0;
// const double pconst = pconst_dim / (rho * Uin * Uin);

// Properties
const double nu = 1.0e-2;               // kinematic viscosity (hypothetical)
const double rho = 1.0;                 // density
const double mu = nu * rho;             // dynamic viscosity

// Algorithmic constants
const double omega = 1.2;               // relaxation factor
const int maxTimesteps = 1e6;           // max time iterations
const int maxPressIters = 1e3;          // max pressure iterations
const double uResidual = 1.0e-12;       // residual criteria for U velocity
const double vResidual = uResidual;     // residual criteria for V velocity
const double pResidual = 1.0e-6;        // residual criteria for pressure
const double dt = 1.0e-3;               // time-step

// String variables
extern std::string velScheme;           // where "up" is upwind, "qk" is QUICK
                                        // and "ct" is central scheme
const std::string fileExt = ".dat";
const std::string filePath = "../data/";
const std::string fileUniqueName = "_test";
                                        // optional unique identifier
const std::string comments = "Test";
                                        // comments for characteristics file

// Intervals on x and y axes
// const double dx = L / static_cast<double>(nx-2);
// const double dy = D / static_cast<double>(ny-2);

// Nondimensional variables
const double Re = uIn * D / nu;         // Reynolds number

// Variables for cross-checking of case validity
// const double entryLength = 0.05 * Re * D;   // 2D channel flow
// const double CFL = (uin * dt / dx) + (vin * dt / dy);
const double CFL = 0.0;


// Nondimensionalization attempt
// const double D = 0.01;
// const double Uin = Re * nu / D;
// const double T = t * D / Uin;
// const double pconst_dim = 100000.0;

// Other global variables
// (defined in navierFVD.cpp)
extern int t;                           // timer variable
extern int pIter;
extern double mChangeMax;
extern double pChangeMax;
extern double uChangeMax;
extern double vChangeMax;
extern int mChangeMax_i;
extern int mChangeMax_j;
extern int pChangeMax_i;
extern int pChangeMax_j;
extern int uChangeMax_i;
extern int uChangeMax_j;
extern int vChangeMax_i;
extern int vChangeMax_j;
// (defined in navierAnalytical.cpp)
// extern double l2Norm;       // 2D channel flow
// extern double pavg_inlet;   // 2D channel flow
// extern double pavg_outlet;  // 2D channel flow
// extern double deltaP;       // 2D channel flow
// extern double uAvg;         // 2D channel flow

extern double scriptRunningTime;

// Coordinates arrays (defined in gridder.cpp)
extern double X[nx+1];
extern double Y[ny+1];
extern double Dx[nx+1];
extern double Dxs[nx+1];
extern double Dy[ny+1];
extern double Dys[ny+1];

// Numerical solution arrays (defined in navierFVD.cpp)
extern double U[nx][ny];
extern double V[nx][ny];
extern double Uo[nx][ny];
extern double Vo[nx][ny];
extern double FU[nx][ny];
extern double FV[nx][ny];
extern double FU1[nx][ny];
extern double FV1[nx][ny];
extern double FU2[nx][ny];
extern double FV2[nx][ny];
extern double P[nx][ny];
extern double MC[nx][ny];
extern double PC[nx][ny];

// Analytical solution arrays (defined in navierAnalytical.cpp)
extern double Ua[ny];       // 2D channel flow
extern double Va[ny];       // 2D channel flow

#endif
