/*
 * main.cpp
 *
 *     Project: Finite Volume Navier-Stokes Solver
 *      Author: Syed Ahmad Raza
 */

#include <string>
#include <iostream>

#include "constants.h"
#include "gridder.h"
#include "filers.h"
#include "printers.h"
#include "navierFVD.h"

/*
g++ -std=c++14 -O3 -Wall gridder.cpp filers.cpp printers.cpp navierFVD.cpp
main.cpp -o main && ./main
*/

void run()
{
    gridder();

    // filerCoordinates();
    filerInfoStart();

    navierFVD();

    filerAllSol();
    filerInfoEnd();
    printerInfoEnd();
}

int main()
{
    run();
    return 0;
}

// int main(int argc, char *argv[])
// {
//     char *input = argv[1];
//     std::string str(input);
//     velScheme = input;
//
//     run();
//     return 0;
// }
