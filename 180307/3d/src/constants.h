/*
 * constants.h
 *
 *  Created on: 2018-01-26
 *      Author: Syed Ahmad Raza
 */

#ifndef CONSANTS_H
#define CONSTANTS_H

#include <string>       // for sring operations

const double pi(3.14159265);

// Domain, range
const int x_y_ratio = 20;
const double ly = 0.01;
const double lz = ly;
const double lx = ly * static_cast<double>(x_y_ratio);

// Dimensions
const double H = ly;
const double B = lz;
const double L = lx;

// Number of nodes on the x and y axes
// Compile-time constants needed for arrays
const int ny = 10;
const int nz = ny;
const int nx = ny * x_y_ratio;

// Intial and/or boundary conditions
// const double uin = Uin;
const double uin = 0.07;
const double vin = 0.0;
const double win = 0.0;
// const double pin = 1.00016e6;
const double pout = 1.0e6;
    // specify pin and pout correctly for Dirichlet pressure conditions
const double ubegin = uin;
const double vbegin = vin;
const double wbegin = vin;
const double pbegin = 1.0e5;

// Properties
const double nu = 1.0e-4;               // kinematic viscosity (hypothetical)
const double rho = 1000.0;              // density
const double mu = nu * rho;             // dynamic viscosity

// Algorithmic constants
const double omega = 0.8;               // relaxation factor
const int maxTimesteps = 1e6;           // max time iterations
const int maxPressIters = 10000;         // max pressure iterations
const double uResidual = 1.0e-6;        // residual criteria for u velocity
const double vResidual = uResidual;     // residual criteria for v velocity
const double wResidual = uResidual;     // residual criteria for w velocity
const double pResidual = 1.0e-6;        // residual criteria for pressure
const double dt = 1.0e-5;               // time-step

// String variables
extern std::string velScheme;
    // where "up" is upwind, "qk" is QUICK scheme and "ct" is central
const std::string fileExt = ".dat";
const std::string filePath = "../data/";

// Nondimensional variables
const double Re = uin * H / nu;         // Reynolds number

// Variables for cross-checking of case validity
const double entryLength = 0.05 * Re * H;
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
extern double uChange;
extern double pChange;
// (defined in navierAnalytical.cpp)
extern double l2Norm;
extern double pavg_inlet;
extern double pavg_outlet;
extern double deltaP;
extern double uAvg;

const int mx = nx;
const int my = ny;
const int mz = nz;

// Coordinates arrays (defined in gridder.cpp)
extern double X[nx+1];
extern double Y[ny+1];
extern double Z[nz+1];

// Intervals on x, y and z axes
extern double Dx[nx+1];
extern double Dxs[nx+1];
extern double Dy[ny+1];
extern double Dys[ny+1];
extern double Dz[nz+1];
extern double Dzs[nz+1];

// Numerical solution arrays (defined in navierFVD.cpp)
extern double U[nx][ny][nz];
extern double V[nx][ny][nz];
extern double W[nx][ny][nz];
extern double Ur[nx][ny][nz];
extern double Vr[nx][ny][nz];
extern double Wr[nx][ny][nz];
extern double Uo[nx][ny][nz];
extern double Vo[nx][ny][nz];
extern double Wo[nx][ny][nz];
extern double FU[nx][ny][nz];
extern double FV[nx][ny][nz];
extern double FW[nx][ny][nz];
extern double FUo[nx][ny][nz];
extern double FVo[nx][ny][nz];
extern double FWo[nx][ny][nz];
extern double P[nx][ny][nz];

// Analytical solution arrays (defined in navierAnalytical.cpp)
extern double Ua[ny];
extern double Va[ny];

#endif
