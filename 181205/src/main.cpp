/*
 * main.cpp
 *
 *     Project: Finite Volume Navier-Stokes Solver
 *      Author: Syed Ahmad Raza
 */

#include <omp.h>        // openMP header
// #include <string>
// #include <iostream>

#include "constants.h"
#include "gridder.h"
#include "filers.h"
#include "printers.h"
#include "solver.h"

/*
g++ -std=c++14 -Wall -O3 -fopenmp gridder.cpp filers.cpp printers.cpp
velschemers.cpp solver.cpp main.cpp -o main && ./main
*/

void run()
{
    #ifdef _OPENMP
    omp_set_num_threads(numOfThreads);
    #endif

    gridder();
    // volumeOfSolid();
    // filerCoordinates();
    filerInfoStart();
    mainSolver();
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
