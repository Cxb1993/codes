/*
 * burgers2D.cpp
 *
 *  Created on: Mar 10, 2017
 *      Author: Ahmad
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

using namespace std;

const double pi(3.14159);

// Parameters for the steady solution
const double a1 = 1.3e13;
const double a2 = 1.3e13;
const double a3 = 0.0;
const double a4 = 0.0;
const double a5 = 1.0;
const double lm = 25.0;
const double x0 = 1.0;
const double re = 50;

// Number of grid points on the x and y axes
const int nx = 50;
const int ny = 50;

// Start and end values for x and y
const double xs = -1.0;
const double xe = 1.0;
const double ys = 0.0;
const double ye = pi/6/lm;

// Size of intervals on the x and y axes
const double dx = (xe - xs) / (nx*1.0 - 1.0);
const double dy = (ye - ys) / (ny*1.0 - 1.0);

// Delaring arrays
double U[nx][ny] = {{0.0}};
double V[nx][ny] = {{0.0}};
double Ur[nx][ny] = {{0.0}};
double Vr[nx][ny] = {{0.0}};

double X[nx] = {0.0};
double Y[ny] = {0.0};

// Initial values of velocities
const double ui = 0.1;
const double vi = 0.0;

// CFL variable
double CFL = 0.0;

// Time values
int t = 0; // time iterator
double dt = 0.000001;

void writeCoordinates()
{
    ofstream fileX; fileX.open("data/coordinatesX.dat");
    ofstream fileY; fileY.open("data/coordinatesY.dat");

    for (int i = 0; i < nx; ++i)
    {
        X[i] = xs + i*dx;
        fileX << X[i] << '\n';
    }
    for (int j = 0; j < ny; ++j)
    {
        Y[j] = ys + j*dy;
        fileY << Y[j] << '\n';
    }
}

void writeFile(string fileNameU, string fileNameV)
{
    ofstream fileU; fileU.open(fileNameU);
    ofstream fileV; fileV.open(fileNameV);

    for (int j = 0; j < ny; ++j)
    {
        for (int i = 0; i < nx; ++i)
        {
            fileU << U[i][j] << '\t';
            fileV << V[i][j] << '\t';
        }
        fileU << '\n';
        fileV << '\n';
    }
}

void analytical()
{
    double x = 0.0;
    double y = 0.0;

    for (int j = 0; j < ny; ++j)
    {
        y = Y[j];
        for (int i = 0; i < nx; ++i)
        {
            x = X[i];

            U[i][j] = (-2*(a2 + a4*y + lm*a5
                      *(exp(lm*(x - x0)) - exp(-lm*(x - x0)))*cos(lm*y)))
                      /
                      (re*(a1 + a2*x + a3*y + a4*x*y + a5
                      *(exp(lm*(x - x0)) + exp(-lm*(x - x0)))*cos(lm*y)));

            V[i][j] = (-2*(a3 + a4*x - lm*a5
                      *(exp(lm*(x - x0)) + exp(-lm*(x - x0)))*sin(lm*y)))
                      /
                      (re*(a1 + a2*x + a3*y + a4*x*y + a5
                      *(exp(lm*(x - x0)) + exp(-lm*(x - x0)))*cos(lm*y)));
        }
    }
}

void updater()
{
    for (int j = 1; j < ny - 1; ++j)
    {
        for (int i = 1; i < nx - 1; ++i)
        {
            U[i][j] += Ur[i][j]*dt;
            V[i][j] += Vr[i][j]*dt;
        }
    }
}

void numerical()
{
    // Boundary conditions are applied implicitly by leaving the existing
    // values from the analytical solution in the array

    // CFL condition
    if (dx <= dy)
        CFL = ui*dt/dx + vi*dt/dx;
    else
        CFL = ui*dt/dy + vi*dt/dy;
    ofstream fileCFL;
    fileCFL.open("data/CFL.dat");
    fileCFL << "CFL = " << CFL << endl;

    // Applying intial conditions
    for (int j = 1; j < ny - 1; ++j)
    {
        for (int i = 1; i < nx - 1; ++i)
        {
            U[i][j] = ui;
            V[i][j] = vi;
        }
    }

    // Computing the numerical solution
    for (t = 1; t <= 100000; ++t)
    {
        for (int j = 1; j < ny - 1; ++j)
        {
            for (int i = 1; i < nx - 1; ++i)
            {
                Ur[i][j] = 1.0/re
                           *( (U[i+1][j] - 2.0*U[i][j] + U[i-1][j])/dx/dx
                            + (U[i][j+1] - 2.0*U[i][j] + U[i][j-1])/dy/dy )
                            - U[i][j]*(U[i+1][j] - U[i-1][j])/2.0/dx
                            - V[i][j]*(U[i][j+1] - U[i][j-1])/2.0/dy;
                Vr[i][j] = 1.0/re
                           *( (V[i+1][j] - 2.0*V[i][j] + V[i-1][j])/dx/dx
                            + (V[i][j+1] - 2.0*V[i][j] + V[i][j-1])/dy/dy )
                            - U[i][j]*(V[i+1][j] - V[i-1][j])/2.0/dx
                            - V[i][j]*(V[i][j+1] - V[i][j-1])/2.0/dy;
            }
        }
        updater();
        cout << t << endl;
        if (t % 10000 == 0)
        {
            ostringstream tStr;
            tStr << t*dt;
            string Uname = ("data/numericalU_" + tStr.str() + ".dat").c_str();
            string Vname = ("data/numericalV_" + tStr.str() + ".dat").c_str();
            writeFile(Uname, Vname);
        }
    }
}

int main()
{
    writeCoordinates();
    analytical();
    writeFile("data/analyticalU.dat", "data/analyticalV.dat");
    numerical();
    return 0;
}

// For single value of x across all y
//double x = -0.766044;
//double y = 0.0;
//
//ofstream fileVt; fileVt.open("data/analyticalVt.dat");
//ofstream fileYt; fileYt.open("data/coordinateYt.dat");
//for (int j = 0; j < ny - 1; ++j)
//{
//    Yt[j] = ys + j*dy;
//    y = Yt[j];
//    Ut[j] = (-2*(a2 + a4*y + lm*a5
//            *(exp(lm*(x - x0)) - exp(-lm*(x - x0)))*cos(lm*y)))
//            /
//            (re*(a1 + a2*x + a3*y + a4*x*y + a5
//            *(exp(lm*(x - x0)) + exp(-lm*(x - x0)))*cos(lm*y)));
//
//    Vt[j] = (-2*(a3 + a4*x - lm*a5
//            *(exp(lm*(x - x0)) + exp(-lm*(x - x0)))*sin(lm*y)))
//            /
//            (re*(a1 + a2*x + a3*y + a4*x*y + a5
//            *(exp(lm*(x - x0)) + exp(-lm*(x - x0)))*cos(lm*y)));
// fileVt << Vt[j] << '\n';
// fileYt << Yt[j] << '\n';
//}
