/*
 * main.cpp
 *
 *  Created on: 2017-09-27
 *      Author: Syed Ahmad Raza
 */

#include <string>
#include <iostream>

#include "constants.h"
#include "gridder.h"
#include "filers.h"
#include "printers.h"
#include "navierFVD.h"
// #include "navierAnalytical.h"

std::string velScheme;

/*
g++ -std=c++14 -O3 -Wall gridder.cpp filers.cpp printers.cpp navierFVD.cpp
main.cpp -o main && ./main up && ./main qk && ./main ct
*/

void run()
{
    gridder();
    infoFilerStart();
    navierFVD();
    // navierAnalytical();
    // navierComparison();
    infoFilerEnd();
    infoPrinter();
}

int main(int argc, char *argv[])
{
    char *input = argv[1];
    std::string str(input);
    velScheme = input;

    run();
    return 0;
}
