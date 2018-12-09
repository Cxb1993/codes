/*
 * constants.h
 *
 *  Created on: 2017-04-17
 *      Author: Syed Ahmad Raza
 */

#ifndef CONSANTS_H
#define CONSTANTS_H

const double pi(3.14159265);

// Domain and range
const int x_y_ratio = 20;
const double ly = 0.01;
const double lx = ly * x_y_ratio * 1.0;   // lx = 0.20

// Number of nodes on the x and y axes
// Compile-time constants needed for arrays
constexpr int ny = 30;              // although "const" is also sufficient
constexpr int nx = ny * x_y_ratio;  // although "const" is also sufficient

// Intervals on x and y axes
const double dx = lx/static_cast<double>(nx-1);
const double dy = ly/static_cast<double>(ny-1);
const double dx1 = dx;
const double dx2 = dx;
const double dy1 = dy;
const double dy2 = dy;

// Boundary conditions
const double uin = 0.1;
const double vin = 0.0;
const double pconst = 1.0e5;


// Other constants
const double nu = 8.9e-4;               // kinematic viscosity
const double rho = 1000;                // density
const double omega = 1.8;               // relaxation factor
const int maxTimesteps = 1000000;       // max time iterations
const int maxPressIterations = 2000;    // max pressure iterations
const double pResidual = 1.0e-4;        // residual criteria for pressure
const double dt = 1.0e-5;               // time-step

const double uResidual = 1.0e-4;
const double vResidual = 1.0e-4;

extern double uMaxChange;               // residual criteria for U velocity
extern double vMaxChange;               // residual criteria for V velocity

// Timer variable
extern int t;                           // timer

// Coordinates arrays
extern double X[nx];    // x coordinates
extern double Y[ny];    // y coordinates

// Solution arrays
extern double U[nx][ny];
extern double V[nx][ny];
extern double Ur[nx][ny];
extern double Vr[nx][ny];
extern double FU[nx][ny];
extern double FV[nx][ny];
extern double FUo[nx][ny];
extern double FVo[nx][ny];
extern double P[nx][ny];


#endif
