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
const int x_y_ratio = 20;
const double ly = 0.01;
const double lx = ly * static_cast<double>(x_y_ratio);

// Dimensions
const double D = ly;
const double L = lx;

// Number of nodes on the x and y axes
// Compile-time constants needed for arrays
const int ny = 51;
const int nx = ny * x_y_ratio;

// Intial and/or boundary conditions
// const double uin = Uin;
const double uin = 0.067;
const double vin = 0.0;
const double pin = 1226800;
const double pout = 1226640;
    // specify pin an pout correctly for Dirichlet pressure conditions
const double ubegin = uin;
const double vbegin = vin;
const double pbegin = 1.0e5;
// const double pconst = pconst_dim / (rho * Uin * Uin);

// Properties
const double nu = 1.0e-4;               // kinematic viscosity (hypothetical)
const double rho = 1000.0;              // density
const double mu = nu * rho;             // dynamic viscosity

// Algorithmic constants
const double omega = 0.8;               // relaxation factor
const int maxTimesteps = 1e6;           // max time iterations
const int maxPressIters = 1000;         // max pressure iterations
const double uResidual = 1.0e-6;        // residual criteria for U velocity
const double vResidual = uResidual;     // residual criteria for V velocity
const double pResidual = 1.0e-7;        // residual criteria for pressure
const double dt = 1.0e-6;               // time-step

// String variables
extern std::string velScheme;
    // where "up" is upwind, "qk" is QUICK scheme and "ct" is central
const std::string fileExt = ".dat";
const std::string filePath = "../data/";

// Intervals on x and y axes
const double dx = L / static_cast<double>(nx - 3);
const double dy = D / static_cast<double>(ny - 3);
const double dx1 = dx;
const double dx2 = dx;
const double dy1 = dy;
const double dy2 = dy;

// Nondimensional variables
const double Re = uin * D / nu;         // Reynolds number

// Variables for cross-checking of case validity
const double entryLength = 0.05 * Re * D;
const double CFL = (uin * dt / dx) + (vin * dt / dy);


// Nondimensionalization attempt
// const double D = 0.01;
// const double Uin = Re * nu / D;
// const double T = t * D / Uin;
// const double pconst_dim = 100000.0;

// Other global variables
// (defined in navierFVD.cpp)
extern int t;                           // timer variable
extern int pIter;
extern double uChange;
extern double pChange;
// (defined in navierAnalytical.cpp)
extern double l2Norm;
extern double pavg_inlet;
extern double pavg_outlet;
extern double deltaP;
extern double uAvg;

// Coordinates arrays (defined in gridder.cpp)
extern double X[2*nx-5];
extern double Y[2*ny-5];
// extern double XS[nx];
// extern double YS[ny];

// Numerical solution arrays (defined in navierFVD.cpp)
extern double U[nx-1][ny-1];
extern double V[nx-1][ny-1];
extern double Ur[nx-1][ny-1];
extern double Vr[nx-1][ny-1];
extern double Uo[nx-1][ny-1];
extern double Vo[nx-1][ny-1];
extern double FU[nx-1][ny-1];
extern double FV[nx-1][ny-1];
extern double FUo[nx-1][ny-1];
extern double FVo[nx-1][ny-1];
extern double P[nx-1][ny-1];

// Analytical solution arrays (defined in navierAnalytical.cpp)
extern double Ua[ny-1];
extern double Va[ny-1];

#endif
