/*
 * ht2d.cpp
 *
 *  Created on: Dec 24, 2016
 *      Author: Syed Ahmad Raza
 */

#include <iostream>
#include <fstream>
#include <sstream>  // for int to string conversion
#include <cmath>
#include <iomanip>

using namespace std;

const double pi(3.141592653589793238462643383279);

// Length and width of the specimen
const double L = 1.0;
const double W = 1.0;

// number of intervals on the x and y axes
const int nx = 50;
const int ny = 50;

// value of thermal diffusivity assumed for pure silver (99.9%)
const double alpha = 1.6563e-4;

const double deltaT = 0.001;        // time step

double X[nx+1] = { 0.0 };           // array for x coordinates
double Xd[nx] = { 0.0 };            // array for x intervals
double Y[ny+1] = { 0.0 };           // array for y coordinates
double Yd[ny] = { 0.0 };            // array for y intervals
double T[nx+1][ny+1] = { 0.0 };     // array for temperature values
double TlastN[ny+1] = { 0.0 };      // array for temperature values of
                                    // numerical solution at a particular
                                    // x cross-section
double TcompA[ny + 1] = { 0.0 };    // array for temperature values of
                                    // analytical solution at a particular
                                    // x cross-section

void gridderSimple()
{
    ofstream fileX;
    fileX.open("./data/coordinatesXSimple.dat");
    for (int x = 0; x <= nx; ++x)
    {
        X[x] = L * (sin( (pi / L) * (x * 1.0 / nx) ));
        fileX << X[x] << '\n';
    }
    for (int x = 0; x < nx; ++x)
    {
        Xd[x] = X[x+1] - X[x];
    }
    ofstream fileY;
    fileY.open("./data/coordinatesYSimple.dat");
    for (int y = 0; y <= ny; ++y)
    {
        Y[y] = W * (sin( (pi / W) * (y * 1.0 / ny) ));
        fileY << Y[y] << '\n';
    }
    for (int y = 0; y < ny; ++y)
    {
        Yd[y] = Y[y+1] - Y[y];
    }
}

void gridderCase01()
{
    ofstream fileX;
    fileX.open("./data/coordinatesXCase01.dat");
    for (int x = 0; x <= nx; ++x)
    {
        X[x] = (0.5 * L) * (1 + sin( pi * ((x * 1.0 / nx) - 0.5)));
        fileX << X[x] << '\n';
    }
    for (int x = 0; x < nx; ++x)
    {
        Xd[x] = X[x+1] - X[x];
    }
    ofstream fileY;
    fileY.open("./data/coordinatesYCase01.dat");
    for (int y = 0; y <= ny; ++y)
    {
        Y[y] = W * (sin( pi * (y * 1.0 / ny) * 0.5));
        fileY << Y[y] << '\n';
    }
    for (int y = 0; y < ny; ++y)
    {
        Yd[y] = Y[y+1] - Y[y];
    }
}

void gridderCase02()
{
    ofstream fileX;
    fileX.open("./data/coordinatesXCase02.dat");
    for (int x = 0; x <= nx; ++x)
    {
        X[x] = L * (sin( pi * (x * 1.0 / nx) * 0.5));
        fileX << X[x] << '\n';
    }
    for (int x = 0; x < nx; ++x)
    {
        Xd[x] = X[x+1] - X[x];
    }
    ofstream fileY;
    fileY.open("./data/coordinatesYCase02.dat");
    for (int y = 0; y <= ny; ++y)
    {
        Y[y] = W * (sin( pi * (y * 1.0 / ny) * 0.5));
        fileY << Y[y] << '\n';
    }
    for (int y = 0; y < ny; ++y)
    {
        Yd[y] = Y[y+1] - Y[y];
    }
}

void solFiler(string fileName)
{
    ofstream file;
    file.open(fileName);
    for (int x = 0; x <= nx; ++x)
    {
        for (int y = 0; y <= ny; ++y)
        {
            file << T[x][y] << '\t';
        }
        file << '\n';
    }
}

void solver(double T_i, double T_xi, double T_xf, double T_yi, double T_yf,
        void (*gridGenerator)(), string fileNameStr)
{
    // Generating the grid
    gridGenerator();

    // Applying the boundary conditions
    // Two vertical boundaries
    for (int y = 0; y <= ny; ++y)
    {
        T[y][0] = T_xi;
        T[y][nx] = T_xf;
    }
    // Two horizontal boundaries
    for (int x = 0; x <= nx; ++x)
    {
        T[0][x] = T_yi;
        T[ny][x] = T_yf;
    }
    // Applying the initial condition
    for (int x = 1; x < nx; ++x)
    {
        for (int y = 1; y < ny; ++y)
            T[y][x] = T_i;
    }
    // Determining the numerical solution until steady state is achieved
    double rmsDiff = 1.0; // for storing the root mean square difference
    int tLast = 0;  // for storing last value of time iterator
    int t = 1;  // time iterator for the while loop
    while ((!(rmsDiff < 0.0000000001)) || (t <= 1000000))
//    for (int t = 1; t <= 4000000; ++t)  // without steady state checking
    {
        for (int x = 1; x < nx; ++x)
        {
            for (int y = 1; y < ny; ++y)
            {
                T[y][x] = T[y][x] + deltaT * alpha *
                        (
                         (
                          ( 2 * (T[y][x+1] + T[y][x-1]) - 4 * T[y][x] )
                           / ( Xd[x] * Xd[x] + Xd[x-1] * Xd[x-1] )
                          - ( 2 * (T[y][x+1] - T[y][x-1]) * (Xd[x] - Xd[x-1]) )
                           / ( (Xd[x] + Xd[x-1])
                            * ( Xd[x] * Xd[x] + Xd[x-1] * Xd[x-1] ) )
                         )
                         +
                         (
                          ( 2 * (T[y+1][x] + T[y-1][x]) - 4 * T[y][x] )
                           / ( Yd[y] * Yd[y] + Yd[y-1] * Yd[y-1] )
                          - ( 2 * (T[y+1][x] - T[y-1][x]) * (Yd[y] - Yd[y-1]) )
                           / ( (Yd[y] + Yd[y-1])
                            * ( Yd[y] * Yd[y] + Yd[y-1] * Yd[y-1] ) )
                         )
                        );
            }
        }
        // Checking the change in all (or many of) the y-values at a
        // particular x position with respect to the last iteration
        double sumDiff = 0.0;
        for (int y = 0; y <= ny; ++y)
        {
            sumDiff += (TlastN[y] - T[y][nx / 2]) * (TlastN[y] - T[y][nx / 2]);
            TlastN[y] = T[y][nx / 2];
        }
        rmsDiff = sqrt(sumDiff / ((ny - 0) + 1));

        // Filing the result at various values of time
        if ( ( t % 100000 ) == 0 )
        {
            ostringstream tStr; tStr << t * deltaT;
            if (t < 1000000)
                solFiler((fileNameStr + "0" + tStr.str() + ".dat").c_str());
            else
                solFiler((fileNameStr + tStr.str() + ".dat").c_str());
        }
        tLast = t++;
    }
    ostringstream tStr; tStr << setprecision(3) << fixed << tLast * deltaT;
    solFiler((fileNameStr + tStr.str() + ".dat").c_str());
}

void solverAnalyticalCase01(double T_1, double T_2,
        void (*gridGenerator)(), string fileNameStr)
{
    gridGenerator();
    double T_3 = T_2 - T_1;

    // Applying the boundary conditions
    // Two vertical boundaries
    for (int y = 0; y <= ny; ++y)
    {
        T[y][0] = T_1;
        T[y][nx] = T_1;
    }
    // Two horizontal boundaries
    for (int x = 0; x <= nx; ++x)
    {
        T[0][x] = T_1;
        T[ny][x] = T_2;
    }

    for (int x = 0; x < nx; ++x)
    {
        for (int y = 0; y < ny; ++y)
        {
            double summation = 0;
            for (int n = 1; n <= 200; ++n)
            {
                double lambda = (n * pi) / L;
                summation += ( (pow(-1, n + 1) + 1) / n )
                        * ( sin(lambda * X[x]) * sinh(lambda * Y[y])
                                / sinh(lambda * W) );
            }
            T[y][x] = ( 2 * T_3 / pi) * summation + T_1;
        }
    }
    // Storing values of temperature for all y-values at the center of
    // x cross-section, for comparision with numerical solution
    for (int y = 0; y <= ny; ++y)
    {
        TcompA[y] = T[y][nx / 2];
    }
    solFiler(fileNameStr.c_str());
}
// Function to compare analytical solution with numerical solution; it must
// always be called after calling the numerical solver and analytical solver
void testSteady()
{
    // Comparing all the y-values in the center of the x cross-section to
    // output the root mean square difference
    double rmsDiff = 0.0; // for storing the root mean square difference
    double sumDiff = 0.0; // for storing the sum of the square differences
    for (int y = 0; y <= ny; ++y)
    {
        sumDiff += (TlastN[y] - TcompA[y]) * (TlastN[y] - TcompA[y]);
    }
    rmsDiff = sqrt(sumDiff / ((ny - 0) + 1));
    ofstream file;
    file.open("./data/steadyComparison.dat");
    file << setprecision(5) << fixed;
    file << "The root mean square difference between analytical and\n"\
            "numerical solution for the first case, at the mid cross-\n"\
            "section of the x-axis is " << rmsDiff << " degree celsius";
}

int main()
{
    solver(25.0, 25.0, 25.0, 25.0, 50.0, gridderCase01, "./data/ht2dCase01T");
    solverAnalyticalCase01(25.0, 50.0, gridderCase01,
            "./data/ht2dCase01Analytical.dat");
    testSteady();
    solver(25.0, 25.0, 50.0, 25.0, 50.0, gridderCase02, "./data/ht2dCase02T");
//    gridderSimple();
//    gridderCase01();
//    gridderCase02();
    return 0;
}

// Alternative main loops for better vectorization:
//for (int x = 1; x < nx; ++x)
//        {
//            for (int y = 1; y < ny; ++y)
//            {
//                Lt[y][x] =
//                        (
//                         (
//                          ( 2 * (T[y][x+1] + T[y][x-1]) - 4 * T[y][x] )
//                           / ( Xd[x] * Xd[x] + Xd[x-1] * Xd[x-1] )
//                          - ( 2 * (T[y][x+1] - T[y][x-1]) * (Xd[x] - Xd[x-1]) )
//                           / ( (Xd[x] + Xd[x-1])
//                            * ( Xd[x] * Xd[x] + Xd[x-1] * Xd[x-1] ) )
//                         )
//                         +
//                         (
//                          ( 2 * (T[y+1][x] + T[y-1][x]) - 4 * T[y][x] )
//                           / ( Yd[y] * Yd[y] + Yd[y-1] * Yd[y-1] )
//                          - ( 2 * (T[y+1][x] - T[y-1][x]) * (Yd[y] - Yd[y-1]) )
//                           / ( (Yd[y] + Yd[y-1])
//                            * ( Yd[y] * Yd[y] + Yd[y-1] * Yd[y-1] ) )
//                         )
//                        );
//            }
//        }
//        for (int x = 1; x < nx; ++x)
//        {
//            for (int y = 1; y < ny; ++y)
//            {
//                T[y][x] = T[y][x] + deltaT * alpha * Lt[y][x];
//            }
//        }
