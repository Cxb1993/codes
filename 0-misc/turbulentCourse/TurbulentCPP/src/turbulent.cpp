/*
 * turbulent.cpp
 *
 *  Created on: Mar 7, 2017
 *      Author: Ahmad
 */

#include <iostream>
#include <fstream>
#include <sstream>  // for int to string conversion
#include <cmath>
#include <iomanip>

using namespace std;

const int nx = 48;
const int ny = 48;
int t = 0;
int maxTimesteps = 2000;

const double lx = 0.010;
const double ly = 0.005;
const double uin = 1.0;
const double vin = 0.0;
const double rho = 1.2;
const double pin = 0.1e6;
const double mu = 2e-5;
const double omega = 1.8;
const double mMax = 1000;
const double epsilon = 1e-2;
const double dt = 0.00001;
const double w = 0.05;

double CFL = 0.0;

const double dx = lx / ((nx - 1) * 1.0);
const double dy = ly / ((ny - 1) * 1.0);

double U[nx][ny] = {0.0};   //array for velocity u
double V[nx][ny] = {0.0};   //array for velocity v

double Ur[nx][ny] = {0.0};  //rhs array for velocity u
double Vr[nx][ny] = {0.0};  //rhs array for velocity v

double P[nx][ny] = {0.0};
double F[nx][ny] = {0.0};

double data[nx][ny] = {0.0};    // array for storing temporary data

void dataDump()
{
    ostringstream tStr; tStr << t * dt;
    ofstream fileU; ofstream fileV;
    fileU.open(("./data/U_" + tStr.str() + ".dat").c_str());
    fileV.open(("./data/V_" + tStr.str() + ".dat").c_str());
    for (int i = 0; i < ny; ++i)
    {
        for (int j = 0; j < nx; ++j)
        {
            fileU << U[i][j] << '\t';
            fileV << V[i][j] << '\t';
        }
        fileU << '\n';
        fileV << '\n';
    }
}

void boundary()
{
    // Left boundary
    int j = 0;
    for (int i = 0; i < ny; ++i)
    {
        U[i][j] = uin;
        V[i][j] = vin;
    }

    // Right boundary
    j = nx - 1;
    for (int i = 0; i < ny; ++i)
    {
        U[i][j] = U[i][j-1];
        V[i][j] = V[i][j-1];
    }

    // Top boundary
    int i = ny - 1;
    for (int j = 0; j < nx; ++j)
    {
        U[i][j] = U[i-1][j];
        V[i][j] = V[i-1][j];
    }

    // Bottom boundary
    i = 0;
    for (int j = 0; j < nx; ++j)
    {
        U[i][j] = U[i+1][j];
        V[i][j] = V[i+1][j];
    }
}

// Object in the flow
void object()
{
    for (int i = 21; i <= 26; ++i)
    {
        for (int j = 21; j <= 31; ++j)
        {
            U[i][j] = 0.0;
            V[i][j] = 0.0;
        }
    }
}

void initial()
{
    for (int i = 0; i < ny; ++i)
    {
        for (int j = 0; j < nx; ++j)
        {
            U[i][j] = uin;
            V[i][j] = vin;
            P[i][j] = pin;
        }
    }
    boundary();
    object();
    dataDump();

    // CFL condition
    if (dx <= dy)
        CFL = dt / (dx / uin);
    else
        CFL = dt / (dy / uin);
    cout << CFL << endl;
}

void momentum()
{
    for (int i = 1; i < (ny - 1); ++i)
    {
        for (int j = 1; j < (nx - 1); ++j)
        {
            Ur[i][j] = mu / rho
                    * ( (U[i][j+1] - 2.0 * U[i][j] + U[i][j-1]) / dx / dx
                      + (U[i+1][j] - 2.0 * U[i][j] + U[i-1][j]) / dy / dy )
                    - U[i][j] * ( (U[i][j+1] - U[i][j-1]) / dx / 2.0 )
                    - V[i][j] * ( (U[i+1][j] - U[i-1][j]) / dy / 2.0 );

            Vr[i][j] = mu / rho
                    * ( (V[i][j+1] - 2.0 * V[i][j] + V[i][j-1]) / dx / dx
                      + (V[i+1][j] - 2.0 * V[i][j] + V[i-1][j]) / dy / dy )
                    - U[i][j] * ( (V[i][j+1] - V[i][j-1]) / dx / 2.0 )
                    - V[i][j] * ( (V[i+1][j] - V[i-1][j]) / dy / 2.0 );

            U[i][j] += Ur[i][j] * dt;
            V[i][j] += Vr[i][j] * dt;
        }
    }
}

void pressure()
{
    // Assign initial values to F and P arrays
    for (int i = 1; i < ny - 1; ++i)
    {
        for (int j = 1; j < nx - 1; ++j)
        {
            F[i][j] = rho * ( (U[i][j+1] - U[i][j-1]) / 2.0 / dx
                            + (V[i+1][j] - V[i-1][j]) / 2.0 / dy ) / dt;
            P[i][j] = pin;
        }
    }

    // Computing Poisson's equation using SOR method.
    double pave1 = 0.0;
    double pave2 = 0.0;
    for (int m = 1; m <= mMax; m++)
    {
        for (int i = 1; i < ny - 1; i++)
        {
            for (int j = 1; j < nx - 1; j++)
            {
                P[i][j] = (1.0 - omega) * P[i][j] + omega
                            * ((P[i + 1][j] + P[i - 1][j]) / dy / dy
                            + (P[i][j + 1] + P[i][j - 1]) / dx / dx - F[i][j])
                            * (dx * dx * dy * dy) / 2.0
                            / (dx * dx + dy * dy);
            }
        }
        pave2 = pave1;
        pave1 = 0.0;
        for (int i = 1; i < ny - 1; i++)
        {
            for (int j = 1; j < nx - 1; j++)
            {
                pave1 = pave1 + P[i][j];
            }
        }
        pave1 = pave1 / (nx * ny * 1.0);
        if (abs(pave1 - pave2) < epsilon)
        {
            break;
        }
    }

    // Computing the corrected velocity field using the calculated pressures
    for (int i = 1; i < ny - 1; ++i)
    {
        for (int j = 1; j < nx - 1; ++j)
        {
            U[i][j] -= dt / rho * (P[i][j+1] - P[i][j-1]) / 2.0 / dx;
            V[i][j] -= dt / rho * (P[i+1][j] - P[i-1][j]) / 2.0 / dy;
        }
    }
}

void filter()
{
    // U filter
    for (int i = 1; i < ny - 1; ++i)
    {
        for (int j = 1; j < nx - 1; ++j)
        {
            data[i][j] = 0.25 * (U[i][j-1] + U[i][j+1]
                               + U[i-1][j] + U[i+1][j]);
        }
    }
    for (int i = 1; i < ny - 1; ++i)
    {
        for (int j = 1; j < nx - 1; ++j)
        {
            U[i][j] = w * data[i][j] + (1.0 - w) * U[i][j];
        }
    }
    // V filter
    for (int i = 1; i < ny - 1; ++i)
    {
        for (int j = 1; j < nx - 1; ++j)
        {
            data[i][j] = 0.25
                    * (V[i][j-1] + V[i][j+1]
                     + V[i-1][j] + V[i+1][j]);
        }
    }
    for (int i = 1; i < ny - 1; ++i)
    {
        for (int j = 1; j < nx - 1; ++j)
        {
            V[i][j] = w * data[i][j] + (1.0 - w) * V[i][j];
        }
    }
    // P filter
    for (int i = 1; i < ny - 1; ++i)
    {
        for (int j = 1; j < nx - 1; ++j)
        {
            data[i][j] = 0.25
                    * (P[i][j-1] + P[i][j+1]
                     + P[i-1][j] + P[i+1][j]);
        }
    }
    for (int i = 1; i < ny - 1; ++i)
    {
        for (int j = 1; j < nx - 1; ++j)
        {
            P[i][j] = w * data[i][j] + (1.0 - w) * P[i][j];
        }
    }
}

int main()
{
    initial();
    for (t = 1; t < maxTimesteps; ++t)
    {
        filter();
        momentum();
        boundary();
        object();
        pressure();
        boundary();
        object();
        if (t % 1 == 0)
        {
            dataDump();
        }
        cout << t << endl;
    }
    return 0;
}
