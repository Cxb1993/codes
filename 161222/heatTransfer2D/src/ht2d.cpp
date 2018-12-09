/*
 * ht2d.cpp
 *
 *  Created on: Dec 22, 2016
 *      Author: Syed Ahmad Raza
 */

#include <iostream>
#include <fstream>
#include <sstream>  // for int to string conversion
#include <cmath>
#include <iomanip>

using namespace std;

//// boundary conditions
//const double T_xi = 25.0;
//const double T_xf = 50.0;
//const double T_yi = 25.0;
//const double T_yf = 50.0;
//
//// initial conditions
//const double T_i = 25.0;

// length and width of the specimen
const double length = 1.0;
const double width = 1.0;

// number of intervals on the x and y axes and their respective values
const int nx = 200;
const int ny = 200;

// value of thermal diffusivity assumed for pure silver (99.9%)
const double alpha = 1.6563e-4;


const double deltaX = length / (nx * 1.0);
const double deltaY = width / (ny * 1.0);
const double deltaT = 0.001;

double T[nx + 1][ny + 1] = { 0.0 };
double Tlast[nx + 1] = { 0.0 };

void gridFiler()
{
    ofstream fileX; ofstream fileY;
    fileX.open("./data/coordinatesX.dat");
    fileY.open("./data/coordinatesY.dat");
    const double x = 0.0;
    const double y = 0.0;
    for (int i = 0; i <= nx; ++i)
    {
        fileX << x + (i * deltaX) << "\t";
    }
    for (int j = 0; j <= nx; ++j)
    {
        fileY << y + (j * deltaY) << "\t";
    }
}

void solFiler(string fileName)
{
    ofstream file;
    file.open(fileName);
    for (int i = 0; i <= nx; ++i)
    {
        for (int j = 0; j <= ny; ++j)
        {
            file << T[i][j] << "\t";
        }
        file << endl;
    }
}

void solver(double T_i, double T_xi, double T_xf, double T_yi, double T_yf,
        string fileNameStr)
{
    // Applying the boundary conditions
    // Two vertical boundaries
    for (int i = 0; i <= ny; ++i)
    {
        T[i][0] = T_xi;
        T[i][ny] = T_xf;
    }
    // Two horizontal boundaries
    for (int j = 0; j <= nx; ++j)
    {
        T[0][j] = T_yi;
        T[ny][j] = T_yf;
    }
    // Applying the initial condition
    for (int i = 1; i < nx; ++i)
    {
        for (int j = 1; j < ny; ++j)
            T[i][j] = T_i;
    }
    // Determining the numerical solution
    double rmsDiff = 1.0; // for storing the root mean square difference
    int tLast = 0;  // for storing last value of time iterator
    int t = 1;  // time iterator for the while loop
    while ((!(rmsDiff < 1e-10)) || (t <= 1000000))
    {
        for (int i = 1; i < nx; ++i)
        {
            for (int j = 1; j < ny; ++j)
            {
                T[i][j] = T[i][j] + alpha * deltaT *
                        (
                         ((T[i + 1][j] - 2 * T[i][j] + T[i - 1][j])
                                 / (deltaX * deltaX))
                         + ((T[i][j + 1] - 2 * T[i][j] + T[i][j - 1])
                                 / (deltaY * deltaY))
                        );
            }
        }
        // Checking the change in all the y-values at a particular x position
        // with respect to the last iteration
        double sumDiff = 0.0;
        for (int i = 0; i <= nx; ++i)
        {
            sumDiff += pow(abs(Tlast[i] - T[nx / 2][i]), 2);
            Tlast[i] = T[nx / 2][i];
        }
        rmsDiff = pow(sumDiff
                / ((nx - 0) + 1), 0.5);

        // Filing the result at various values of time
        if ((t == static_cast<int>(1e5)) || ((t % static_cast<int>(5e5)) == 0))
        {
            ostringstream tStr; tStr << t * deltaT;
            solFiler((fileNameStr + tStr.str() + ".dat").c_str());
        }
        tLast = t++;
    }
    cout << setprecision(5) << fixed << tLast << "\t" << (tLast) * deltaT;
    ostringstream tStr; tStr << setprecision(3) << fixed << (tLast) * deltaT;
    solFiler((fileNameStr + tStr.str() + ".dat").c_str());
}

int main()
{
    gridFiler();
    solver(25.0, 25.0, 25.0, 25.0, 50.0, "./data/ht2dCase01T");
//    solver(25.0, 25.0, 50.0, 25.0, 50.0, "./data/ht2dCase02T");
}
