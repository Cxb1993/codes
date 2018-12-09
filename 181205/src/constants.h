/*
 * constants.h
 *
 *     Project: Finite Volume Navier-Stokes Solver
 *      Author: Syed Ahmad Raza
 */

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <string>       // for string operations
#include "Matrix11.h"   // for Matrix definitions

const double pi(3.14159265);

// Dimensions of the computational domain and solid objects
const double L = 1.0;
const double B = 1.0;
const double H = 1.0;
const double lx = L;
const double ly = B;
const double lz = H;

const double R = 0.25;                  // radius of spherical solid body
const double Xc = 0.5;                  // x-coordinate of solid body center
const double Yc = 0.5;                  // y-coordinate of solid body center
const double Zc = 0.5;                  // z-coordinate of solid body center

// Number of nodes on the x, y and z axes
// Compile-time constants needed for native C-style arrays
const int nd = 40;                      // grid dimension
const int ghostCells = 4;               // total ghost cells in an axis
const int nx = nd + ghostCells;
const int ny = static_cast<int>(nx * ly/lx);
const int nz = static_cast<int>(nx * lz/lx);
const int nSubGrids = 20;

// Intial and/or boundary conditions
// const double uin = Uin;
const double uIn = 1.0;
const double vIn = 0.0;
const double wIn = 0.0;
// const double pIn = 1.00016e6;
// const double pOut = 1.0e6;
    // specify pIn and pOut correctly for Dirichlet pressure conditions
const double pInitial = 0.0;
const double uInitial = 0.0;
const double vInitial = 0.0;
const double wInitial = 0.0;

// Solid body velocity
const double uS = 0.0;
const double vS = 0.0;
const double wS = 0.0;

// Properties
// const double nu = 1.0e-2;               // kinematic viscosity (hypothetical)
// const double rho = 1.0;                 // density
// const double mu = nu * rho;             // dynamic viscosity

// Algorithmic constants
const double omega = 1.2;               // relaxation factor
const int maxTimesteps = 1e6;           // max time iterations
const int maxPrIters = 1e3;             // max pressure iterations
const double uResidual = 1.0e-12;        // residual criteria for U velocity
const double vResidual = uResidual;     // residual criteria for V velocity
const double wResidual = uResidual;     // residual criteria for V velocity
const double pResidual = 1.0e-8;        // residual criteria for pressure
const double dt = 1.0e-5;               // time-step
const int numOfThreads = 4;             // number of threads for omp

// Ranges of x, y and z direction indexes for solid object
const int iBegVOS = 2;
const int iEndVOS = nx-3;
const int jBegVOS = 2;
const int jEndVOS = ny-3;
const int kBegVOS = 2;
const int kEndVOS = nz-3;

// String variables
const std::string velScheme = "quick";  // options are: "upwind", "quick"
const std::string filePath = "../data/";
const std::string fileUniqueName = "test_Re1";
                                        // optional unique identifier
const std::string comments =
"Testing very low Renolds numbers";
                                        // comments for characteristics file

// Nondimensional variables
const double Re = 10;                   // Reynolds number

// Variables for cross-checking of case validity
// const double CFL = (uin * dt / dx) + (vin * dt / dy);
const double CFL = 0.0;

// Other global variables
// (defined in solver.cpp)
extern int t;                           // timer variable
extern int pIter;                       // number of pr. iterations
extern int tPrConvrg;                   // 1st timestep when pIter < maxPrIters
extern double mChangeMax;
extern double pChangeMax;
extern double uChangeMax;
extern double vChangeMax;
extern double wChangeMax;

extern double scriptRunningTime;

// Virtual force and force coefficients
extern double fxTotal;
extern double fyTotal;
extern double fzTotal;

// Coordinates arrays (defined in gridder.cpp)
extern Numeric_lib::Matrix<double,1> X;
extern Numeric_lib::Matrix<double,1> Y;
extern Numeric_lib::Matrix<double,1> Z;
extern Numeric_lib::Matrix<double,1> Xa;
extern Numeric_lib::Matrix<double,1> Ya;
extern Numeric_lib::Matrix<double,1> Za;
extern Numeric_lib::Matrix<double,1> Xs;
extern Numeric_lib::Matrix<double,1> Ys;
extern Numeric_lib::Matrix<double,1> Zs;
extern Numeric_lib::Matrix<double,1> Dx;
extern Numeric_lib::Matrix<double,1> Dxs;
extern Numeric_lib::Matrix<double,1> Dy;
extern Numeric_lib::Matrix<double,1> Dys;
extern Numeric_lib::Matrix<double,1> Dz;
extern Numeric_lib::Matrix<double,1> Dzs;
extern Numeric_lib::Matrix<double,3> ETA;

// Numerical solution arrays (defined in solver.cpp)
extern Numeric_lib::Matrix<double,3> U;
extern Numeric_lib::Matrix<double,3> V;
extern Numeric_lib::Matrix<double,3> W;
extern Numeric_lib::Matrix<double,3> Uo;
extern Numeric_lib::Matrix<double,3> Vo;
extern Numeric_lib::Matrix<double,3> Wo;
extern Numeric_lib::Matrix<double,3> FU;
extern Numeric_lib::Matrix<double,3> FV;
extern Numeric_lib::Matrix<double,3> FW;
extern Numeric_lib::Matrix<double,3> FU1;
extern Numeric_lib::Matrix<double,3> FV1;
extern Numeric_lib::Matrix<double,3> FW1;
extern Numeric_lib::Matrix<double,3> FU2;
extern Numeric_lib::Matrix<double,3> FV2;
extern Numeric_lib::Matrix<double,3> FW2;
extern Numeric_lib::Matrix<double,3> P;
extern Numeric_lib::Matrix<double,3> FX;
extern Numeric_lib::Matrix<double,3> FY;
extern Numeric_lib::Matrix<double,3> FZ;
extern Numeric_lib::Matrix<double,1> FXY;
extern Numeric_lib::Matrix<double,2> FXZ;
extern Numeric_lib::Matrix<double,1> FYY;
extern Numeric_lib::Matrix<double,2> FYZ;
extern Numeric_lib::Matrix<double,1> FZY;
extern Numeric_lib::Matrix<double,2> FZZ;

#endif
