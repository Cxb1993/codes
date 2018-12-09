/*
 * constants.h
 *
 *  Created on: 2017-06-13
 *      Author: Syed Ahmad Raza
 */

#ifndef CONSANTS_H
#define CONSTANTS_H

const double pi(3.14159265);

// Domain, range
const int x_y_ratio = 20;
const double ly = 1.0;
const double lx = ly * x_y_ratio * 1.0;   // lx = 20

// Number of nodes on the x and y axes
// Compile-time constants needed for arrays
constexpr int ny = 30;              // although "const" is also sufficient
constexpr int nx = ny * x_y_ratio;  // although "const" is also sufficient

// Intervals on x and y axes
const double dx = lx / static_cast<double>(nx-1);
const double dy = ly / static_cast<double>(ny-1);
const double dx1 = dx;
const double dx2 = dx;
const double dy1 = dy;
const double dy2 = dy;

// Boundary conditions
const double uin = 1.0;
const double vin = 0.0;
const double pconst = 10000.0;

// Other constants
const double Re = 1;                    // Reynolds number
const double omega = 1.8;               // relaxation factor
const int maxTimesteps = 10000;         // max time iterations
const int maxPressIterations = 1000;    // max pressure iterations
const double pResidual = 1.0e-4;        // residual criteria for pressure
const double uResidual = 1.0e-6;        // residual criteria for U velocity
const double vResidual = 1.0e-6;        // residual criteria for V velocity
const double dt = 1.0e-5;               // time-step

extern double uMaxChange;
extern double vMaxChange;

// Timer variable
extern int t;                           // timer

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

extern double Ua[nx][ny];
extern double Va[nx][ny];
extern double Pavg[ny];

// Analytical solution arrays (defined in navierAnalytical.cpp)
extern double Ua[nx][ny];
extern double Va[nx][ny];

// Properties
const double nu = 8.9e-7;               // kinematic viscosity
const double rho = 1000.0;              // density
const double mu = nu * rho;             // dynamic viscosity

// Dimensional variables
const double D = 0.01;
const double Uin = 0.1;
const double T = D/Uin;
const double L = lx * D;
const double pconst_dim = pconst * Uin * Uin * rho;

#endif
