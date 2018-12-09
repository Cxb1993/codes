/*
 * burgers2D.cpp
 *
 *  Created on: Mar 3, 2017
 *      Author: Syed Ahmad Raza
 */

#include <iostream>
#include <fstream>
#include <sstream>  // for int to string conversion
#include <cmath>


using namespace std;

const double pi(3.14159);

//// Length and width of the specimen
//const double L = 1.0;
//const double W = 1.0;

// Number of intervals on the x and y axes
const int nx = 50;
const int ny = 50;

// Delaring arrays
double U[ny + 1][nx + 1] = {0.0};
double V[ny + 1][nx + 1] = {0.0};
double Ur[ny + 1][nx + 1] = {0.0};
double Vr[ny + 1][nx + 1] = {0.0};

// Initial values of velocities
const double u_i = 0.1;
const double v_i = 0.1;

// Parameters for the steady solution
const double a1 = 1.3e13;
const double a2 = 1.3e13;
const double a3 = 0.0;
const double a4 = 0.0;
const double a5 = 1.0;
const double lm = 25.0;
const double x0 = 1.0;
const double re = 500;

// Start and end values for x and y
const double xs = -1.0;
const double xe = 1.0;
const double ys = 0.0;
const double ye = 2.0;

// Interval width on the x axis
const double deltaX = (xe - xs) / (nx * 1.0);
const double deltaY = (ye - ys) / (ny * 1.0);

// Time step
const double deltaT = 0.001;

void solverAnalytical()
{
    ofstream fileU;
    ofstream fileV;
    fileU.open("./data/analyticalU.dat");
    fileV.open("./data/analyticalV.dat");

    int x = 0.0; // intializing x variable
    int y = 0.0; // intializing y variable
    for (int i = 0; i <= ny; ++i)
    {
        y = ys + i * deltaY;
        for (int j = 0 ; j <= nx; ++j)
        {
            x = xs + j * deltaX;
            U[i][j] =
                    (
                    -2 * (a2 + a4 * y + lm * a5 * (exp(lm * (x - x0))
                            - exp(-lm * (x - x0))) * cos(lm * y))
                    )
                    /
                    (
                    re * (a1 + a2 * x + a3 * y + a4 * x * y + a5
                            * (exp(lm * (x - x0)) + exp(-lm * (x - x0)))
                            * cos(lm * y) )
                    );
            V[i][j] =
                    (
                    -2 * (a3 + a4 * x - lm * a5 * (exp(lm * (x - x0))
                            + exp(-lm * (x - x0))) * sin(lm * y))
                    )
                    /
                    (
                    re * (a1 + a2 * x + a3 * y + a4 * x * y
                            + a5 * (exp(lm * (x - x0)) + exp(-lm * (x - x0)))
                            * cos(lm * y) )
                    );
            fileU << U[i][j] << '\t';
            fileV << V[i][j] << '\t';
        }
        fileU << '\n';
        fileV << '\n';
    }
}

void solverNumerical()
{
    // Applying the initial conditions
    for (int i = 1; i < ny; ++i)
    {
        for (int j = 1; j < nx; ++j)
        {
            U[i][j] = u_i; V[i][j] = v_i;
        }
    }
    // Main loop for numerical solution
    for (int t = 1; t <= 5000000; ++t)
    {
        for (int i = 1; i < ny; ++i)
        {
            for (int j = 1; j < nx; ++j)
            {
                Ur[i][j] = (1.0 / re)
                 * ((U[i][j+1] - 2.0 * U[i][j] + U[i][j-1]) / (deltaX * deltaX)
                  + (U[i+1][j] - 2.0 * U[i][j] + U[i-1][j]) / (deltaY * deltaY))
                 - U[i][j] * ((U[i][j+1] - U[i][j-1]) / (2.0 * deltaX))
                 - V[i][j] * ((U[i+1][j] - U[i-1][j]) / (2.0 * deltaY));

                Vr[i][j] = (1.0 / re)
                 * ((V[i][j+1] - 2.0 * V[i][j] + V[i][j-1]) / (deltaX * deltaX)
                  + (V[i+1][j] - 2.0 * V[i][j] + V[i-1][j]) / (deltaY * deltaY))
                 - U[i][j] * ((V[i][j+1] - V[i][j-1]) / (2.0 * deltaX))
                 - V[i][j] * ((V[i+1][j] - V[i-1][j]) / (2.0 * deltaY));
            }
        }
        for (int i = 0; i <= ny; ++i)
        {
            for (int j = 0; j <= nx; ++j)
            {
                U[i][j] += deltaT * Ur[i][j];
                V[i][j] += deltaT * Vr[i][j];
            }
        }

        // Filing the result at various values of time
        if ( ( t % 1000000 ) == 0 )
        {
            ostringstream tStr; tStr << t * deltaT;
            ofstream fileU; ofstream fileV;
            fileU.open(("./data/numericalU_" + tStr.str() + ".dat").c_str());
            fileV.open(("./data/numericalV_" + tStr.str() + ".dat").c_str());
            for (int i = 0; i <= ny; ++i)
            {
                for (int j = 0; j <= nx; ++j)
                {
                    fileU << U[i][j] << '\t'; fileV << V[i][j] << '\t';
                }
                fileU << '\n'; fileV << '\n';
            }
        }
        cout << t << '\n';
    }
}

int main()
{
    solverAnalytical();
    solverNumerical();
    return 0;
}
