/*
 * sequential_dotProduct.cpp
 *
 *     Project: Testing OpenMP
 *      Author: Syed Ahmad Raza
 * Description: This sequential program multiplies the individual elements of
 *              two arrays and saves the result in the variable sum; sum is a
 *              so-called reduction variable.
 */

// #include <cmath>        // math functions
#include <iostream>     // functions for input and output to console
#include <ctime>        // to time the script
#include <chrono>       // to measure and display the elapsed time
#include <omp.h>        // openMP header
#include "Matrix11.h"   // for Matrix definitions

using namespace std;

void runner();

int n = 1000;
double sum;
Numeric_lib::Matrix<double,2> a(n,n), b(n,n);

int main(int argc, char *argv[])
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            a[i] = i * j * 0.00001;
            b[i] = i * j * 0.000001;
        }
    }

    std::chrono::steady_clock::time_point
        tStart = std::chrono::steady_clock::now();

    #ifdef _OPENMP
    omp_set_num_threads(5);
    #endif

    #pragma omp parallel
    {
        // runner();
        int a = 0;
        cout << endl << omp_get_thread_num();
        while (a < 5)
        {
            runner();
            ++a;
            // cout << endl << omp_get_num_threads();
        }
    }

    std::chrono::steady_clock::time_point tEnd
        = std::chrono::steady_clock::now();
    double scriptRunningTime =
        std::chrono::duration_cast<std::chrono::microseconds>
        (tEnd - tStart).count();
    printf("\ntime = %f\n", scriptRunningTime);
}

// #ifdef _OPENMP
// void omp_set_num_threads(int num_threads);
// #endif

void runner()
{
    sum = 0;

    #pragma omp for reduction(+:sum)
    for (int i = 0; i < n; i++ )
    {
        for (int j = 0; j < n; j++)
        {
            sum = sum + a[i][j]*b[i][j];
        }
    }
    #pragma omp single
    {
        printf("\nsum = %f\n", sum);
    }
}
