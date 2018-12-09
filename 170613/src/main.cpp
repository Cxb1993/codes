/*
 * main.cpp
 *
 *  Created on: 2017-06-13
 *      Author: Syed Ahmad Raza
 */

#include <cmath>        // math functions
// #include <algorithm>    // functions for ranges of elements (arrays)

#include "constants.h"
#include "gridder.h"
#include "filers.h"
#include "navierFVD.h"
#include "navierAnalytical.h"

using namespace std;

// g++ -std=c++14 gridder.cpp filers.cpp navierFVD.cpp main.cpp -o main & ./main

int main()
{
    navierFVD();
    return 0;
}
