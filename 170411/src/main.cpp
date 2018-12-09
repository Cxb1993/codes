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

int main()
{
    gridder();
    gridFiler(X, nx, "../data/coordinateX.dat");
    gridFiler(Y, ny, "../data/coordinateY.dat");
    initial();
    boundary();
    solFiler(U, nx, ny, "../data/navierFVD_U_0.dat");
    solFiler(V, nx, ny, "../data/navierFVD_V_0.dat");
    for (int t = 1; t <= maxTimesteps; ++t)
    {
        cout << t*dt << '\n';
        boundary();
        momentum();
        momentumTimeStep();
        pressure();
        velocityCorrector();
        if (t % 100 == 0)
        {
            ostringstream tStr;
            tStr << t * dt;
            solFiler(U, nx, ny, "../data/navierFVD_U_" + tStr.str() + ".dat");
            solFiler(V, nx, ny, "../data/navierFVD_V_" + tStr.str() + ".dat");
        }
        cout << "\n\n";
    }
    return 0;
}
