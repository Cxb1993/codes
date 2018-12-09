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

    for (t = 1; t <= maxTimesteps; ++t)
    {
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

        if (t % 1 == 0)
        {
            ostringstream tStr;
            // tStr << time;
            tStr << t;
            solFiler(U, nx, ny, "../data/navierFVD_U_" + tStr.str() + ".dat");
            solFiler(V, nx, ny, "../data/navierFVD_V_" + tStr.str() + ".dat");
            solFiler(P, nx, ny, "../data/navierFVD_P_" + tStr.str() + ".dat");
        }

        cout << "\n\n";
    }

    return 0;
}
