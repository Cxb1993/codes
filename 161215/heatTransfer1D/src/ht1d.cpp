/*
 * ht1d.cpp
 *
 *  Created on: Dec 15, 2016
 *      Author: Syed Ahmad Raza
 */

#include <iostream>
#include <fstream>  // for file manipulation
#include <iomanip>  // for std::precision
#include <sstream>  // for integer-to-string conversion
#include <string>
#include <cmath>

using namespace std;

// value of thermal diffusivity assumed for pure silver (99.9%)
const double ALPHA = 1.6563e-4;

const double timeStep = 0.001;

// boundary conditions
const double T_xi = 25.0;
const double T_xf = 50.0;

// initial conditions
const double T_i = 25.0;

// length of the specimen
const double length = 1.0;

double solver(int t, double deltaT, int n, int returnIndex);
void solverFiler(int t, double deltaT, int n);

int main()
{
    int n = 158;
    solverFiler(1e2, timeStep, n);
    solverFiler(1e3, timeStep, n);
    solverFiler(1e4, timeStep, n);
    solverFiler(1e5, timeStep, n);
    solverFiler(1e6, timeStep, n);
    solverFiler(2.1e6, timeStep, n);
    solverFiler(3e6, timeStep, n);
    solverFiler(3.5e6, timeStep, n);
    return 0;
}

void solverFiler(int t, double deltaT, int n)
{
    const double deltaX = length / (n * 1.0);
    double T[n + 1] = { 0.0 };
    double x = 0.0;

    T[0] = T_xi;
    T[n] = T_xf;

    // Applying the initial condition to the array (first row of the array);
    // where i is index for columns, indicating location on x axis
    for (int i = 1; i < n; ++i)
    {
        T[i] = T_i;
    }

    for (int j = 1; j <= t; ++j)
    {
        for (int i = 1; i < n; ++i)  // first and last columns are bc's
        {
            T[i] = T[i] + (ALPHA * deltaT) / (deltaX * deltaX)
                            * (T[i + 1] - 2 * T[i] + T[i - 1]);
        }
    }
//    // Console output
//    for (int i = 0; i <= n; ++i)
//    {
//        cout << x + (i * deltaX) << "\t" << T[i] << "\n";
//    }
    // File writing
    ofstream file;
    ostringstream convertN; convertN << n;
    ostringstream convertT; convertT << (t * deltaT);
    file.open(("./data/ht1dn" + convertN.str() + "t" + convertT.str()
            + ".dat").c_str());
    for (int i = 0; i <= n; ++i)
    {
        file << x + (i * deltaX) << "\t" << T[i] << "\n";
    }
    file.close();
}

double solver(int t, double deltaT, int n, int returnIndex)
{
    const double deltaX = length / (n * 1.0);
    double T[n + 1] = { 0.0 };

    T[0] = T_xi;
    T[n] = T_xf;

    // Applying the initial condition to the array (first row of the array);
    // where i is index for columns, indicating location on x axis
    for (int i = 1; i < n; ++i)
    {
        T[i] = T_i;
    }

    for (int j = 1; j <= t; ++j)
    {
        for (int i = 1; i < n; ++i)  // first and last columns are bc's
        {
            T[i] = T[i] + (ALPHA * deltaT) / (deltaX * deltaX)
                            * (T[i + 1] - 2 * T[i] + T[i - 1]);
        }
    }
    return T[returnIndex];
}

// Grid independence is achieved at 160 intervals and steady state is achieved
// at almost 3.5e6 total time, meaning t = 3500

////Code for checking time taken to achieve steady state
//int n = 160;
//int nMid = n / 2;
//int initialTime = 1e6;
//double lastTemp = solver(initialTime, timeStep, n, nMid);
//int nextTime = initialTime + 500000;
//double nextTemp = solver(nextTime, timeStep, n, nMid);
//while (abs(nextTemp - lastTemp) >= 0.1)
//{
//    lastTemp = nextTemp;
//    nextTime += 500000;
//    nextTemp = solver(nextTime, timeStep, n, nMid);
//    cout << nextTime << "\t" << nextTemp << endl;
//}

////Code for checking grid independence
//int time = 3.5e6;
//int xInitialIntervals = 150;
//double lastTemp = solver(time, timeStep, xInitialIntervals,
//        xInitialIntervals - 1);
//int xNextIntervals = xInitialIntervals + 1;
//double nextTemp = solver(time, timeStep, xNextIntervals, xNextIntervals - 1);
//while (abs(nextTemp - lastTemp) >= 10e-4)
//{
//    lastTemp = nextTemp;
//    ++xNextIntervals;
//    nextTemp = solver(time, timeStep, xNextIntervals, xNextIntervals - 1);
//    cout << xNextIntervals << "\t" << nextTemp << endl;
//}
