/*
 * constants.h
 *
 *  Created on: 2017-09-27
 *      Author: Syed Ahmad Raza
 */

#ifndef CONSANTS_H
#define CONSTANTS_H

const double pi(3.14159265);

// Timer variable
extern int t;

// Domain, range
const int x_y_ratio = 20;
const double ly = 0.01;
const double lx = ly * static_cast<double>(x_y_ratio);

// Dimensions
const double D = ly;
const double L = lx;

// Number of nodes on the x and y axes
// Compile-time constants needed for arrays
constexpr int ny = 60;                  // although "const" is also sufficient; also, it must be an even number for compatibility with navierAnaltical()
constexpr int nx = ny * x_y_ratio;      // although "const" is also sufficient

// Intervals on x and y axes
const double dx = L / static_cast<double>(nx-1);
const double dy = D / static_cast<double>(ny-1);
const double dx1 = dx;
const double dx2 = dx;
const double dy1 = dy;
const double dy2 = dy;

// Intial and/or boundary conditions
// const double uin = Uin;
const double uin = 0.1;
const double vin = 0.0;
// const double pconst = pconst_dim / (rho * Uin * Uin);
const double pconst = 1.0e5;

// Properties
const double nu = 1.0e-4;               // kinematic viscosity (hypothetical)
const double rho = 1000.0;              // density
const double mu = nu * rho;             // dynamic viscosity

// Algorithmic constants
const double omega = 1.8;               // relaxation factor
const int maxTimesteps = 200000;        // max time iterations
const int maxPressIterations = 2000;    // max pressure iterations
const double uResidual = 1.0e-6;        // residual criteria for U velocity
const double vResidual = 1.0e-6;        // residual criteria for V velocity
const double pResidual = 1.0e-4;        // residual criteria for pressure
const double dt = 1.0e-5;               // time-step

// Nondimensional variables
const double Re = uin * D / nu;        // Reynolds number

// Nondimensionalization attempt
// const double D = 0.01;
// const double Uin = Re * nu / D;
// const double T = t * D / Uin;
// const double pconst_dim = 100000.0;

// Comparison variables
extern double uMaxChange;
extern double vMaxChange;
extern double l2Norm;

// Coordinates arrays (defined in gridder.cpp)
extern double X[nx];                    // x coordinates
extern double Y[ny];                    // y coordinates

// Numerical solution arrays (defined in navierFVD.cpp)
extern double U[nx][ny];
extern double V[nx][ny];
extern double Ur[nx][ny];
extern double Vr[nx][ny];
extern double Uo[nx][ny];
extern double Vo[nx][ny];
extern double FU[nx][ny];
extern double FV[nx][ny];
extern double FUo[nx][ny];
extern double FVo[nx][ny];
extern double P[nx][ny];

// Analytical solution arrays (defined in navierAnalytical.cpp)
extern double Pavg[ny];
extern double Ulast[ny];
extern double Ua[ny];
extern double Va[ny];

#endif
