/*
 * jacobi.cpp
 *
 *  Created on: Mar 27, 2017
 *      Author: Syed Ahmad Raza
 */

#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

// Indexers
int i = 0;
int j = 0;

int k = 0;

int n = 3;

double A[n][n] = {
                  {3.0, -0.1, -0.2},
                  {0.1, 7.0, -0.3},
                  {0.3, -0.2, 10.0}
                 };
double B[n] = {7.85, -19.3, 71.4};
double X[n] = {0.0};
double Xi[n] = {0.0};

void jacobi()
{
    for (k = 0; k < 10; ++k)
    {
        for (i = 0; i < n; ++i)
        {
            Xi[i] = 0.0;
            for (j = 0; j < 3; ++j)
            {
                if (i != j)
                    Xi[i] += A[i][j]*X[j];
            }
            Xi[i] = (B[i] - Xi[i]) / A[i][i];
        }
        for (i = 0; i < n; ++i)
        {
            X[i] = Xi[i];
            cout << X[i] << '\t';
        }
        cout << '\n';
    }
}
