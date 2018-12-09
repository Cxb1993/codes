/*
 * main.cpp
 *
 *  Created on: 2017-04-17
 *      Author: Syed Ahmad Raza
 */

#include <cmath>        // math functions
#include <iostream>     // functions for input and output to console
#include <fstream>      // functions for file input and output
#include <sstream>      // functions for string conversion
#include <algorithm>    // functions for ranges of elements (arrays)

#include "constants.h"
#include "gridder.h"
#include "filers.h"
#include "navierFVD.h"

using namespace std;

int t = 0;
double uMaxChange = 1.0;
double vMaxChange = 1.0;

int main()
{
    gridder();
    gridFiler(X, nx, "../data/coordinateX.dat");
    gridFiler(Y, ny, "../data/coordinateY.dat");

    // double time = t * dt;
    // cout << time << '\n';
    cout << t << '\n';
    initial();
    velBoundary();
    pressBoundary();
    solFiler(U, nx, ny, "../data/navierFVD_U_0.dat");
    solFiler(V, nx, ny, "../data/navierFVD_V_0.dat");
    solFiler(P, nx, ny, "../data/navierFVD_P_0.dat");

    t = 1;
    while   (   t <= 100 ||
             (
                 t <= maxTimesteps &&
                 (uMaxChange > uResidual || vMaxChange > vResidual)
             )
            )
    {
        uMaxChange = 0.0;
        vMaxChange = 0.0;
        // time = t * dt;
        // cout << time << '\n';
        cout << t << '\n';

        momentum();
        momentumTimeStep();
        velBoundary();
        pressure();
        pressBoundary();
        velUpdater();
        velBoundary();

        if ((t <= 2) || (t % 100 == 0))
        {
            ostringstream tStr;
            // tStr << time;
            tStr << t;
            solFiler(U, nx, ny, "../data/navierFVD_U_" + tStr.str() + ".dat");
            solFiler(V, nx, ny, "../data/navierFVD_V_" + tStr.str() + ".dat");
            solFiler(P, nx, ny, "../data/navierFVD_P_" + tStr.str() + ".dat");
        }
        cout << '\n' << uMaxChange << '\t' << vMaxChange;
        cout << "\n\n";
        ++t;
    }
    ostringstream tStr;
    // tStr << time;
    tStr << --t;
    solFiler(U, nx, ny, "../data/navierFVD_U_" + tStr.str() + ".dat");
    solFiler(V, nx, ny, "../data/navierFVD_V_" + tStr.str() + ".dat");
    solFiler(P, nx, ny, "../data/navierFVD_P_" + tStr.str() + ".dat");

    return 0;
}
