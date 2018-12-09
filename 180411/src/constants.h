/*
 * constants.h
 *
 *  Created on: 2017-09-27
 *      Author: Syed Ahmad Raza
 */

#ifndef CONSANTS_H
#define CONSTANTS_H

#include <string>       // for sring operations
#include "Matrix11.h"   // for Matrix definitions

const double pi(3.14159265);

// Domain, range
const int x_y_ratio = 1;
const double ly = 1.0;
const double lx = ly * static_cast<double>(x_y_ratio);

// Dimensions
const double D = ly;
const double L = lx;

// Number of nodes on the x and y axes
// Compile-time constants needed for native C-style arrays
const int nxy = 41;                     // grid dimension
const int ghostCells = 4;               // total ghost cells in an axis
const int ny = nxy + ghostCells;
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
const double pResidual = 1.0e-8;        // residual criteria for pressure
const double dt = 1.0e-4;               // time-step

// String variables
extern std::string velScheme;           // where "up" is upwind, "qk" is QUICK
                                        // and "ct" is central scheme
const std::string fileExt = ".dat";
const std::string filePath = "../data/";
const std::string fileUniqueName = "_050213";
                                        // optional unique identifier
const std::string comments = "Average velocity conditions";
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
extern Numeric_lib::Matrix<double,1> X;
extern Numeric_lib::Matrix<double,1> Y;
extern Numeric_lib::Matrix<double,1> Xa;
extern Numeric_lib::Matrix<double,1> Ya;
extern Numeric_lib::Matrix<double,1> Dx;
extern Numeric_lib::Matrix<double,1> Dxs;
extern Numeric_lib::Matrix<double,1> Dy;
extern Numeric_lib::Matrix<double,1> Dys;

// Numerical solution arrays (defined in navierFVD.cpp)
extern Numeric_lib::Matrix<double,2> U;
extern Numeric_lib::Matrix<double,2> V;
extern Numeric_lib::Matrix<double,2> Uo;
extern Numeric_lib::Matrix<double,2> Vo;
extern Numeric_lib::Matrix<double,2> FU;
extern Numeric_lib::Matrix<double,2> FV;
extern Numeric_lib::Matrix<double,2> FU1;
extern Numeric_lib::Matrix<double,2> FV1;
extern Numeric_lib::Matrix<double,2> FU2;
extern Numeric_lib::Matrix<double,2> FV2;
extern Numeric_lib::Matrix<double,2> P;
extern Numeric_lib::Matrix<double,2> MC;
extern Numeric_lib::Matrix<double,2> PC;

// // 2D channel flow
// // Analytical solution arrays (defined in navierAnalytical.cpp)
// // extern Numeric_lib::Matrix<double,1> Ua;
// // extern Numeric_lib::Matrix<double,1> Va;
// extern double Ua[ny];
// extern double Va[ny];

#endif
