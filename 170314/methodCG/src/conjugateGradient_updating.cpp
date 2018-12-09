/*
 * conjugateGradient.cpp
 *
 *  Created on: Mar 14, 2017
 *      Author: Syed Ahmad Raza
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

using namespace std;

const int ii = 3;
const int jj = ii;
const int kk = 1;

double C[ii][kk] = {{0.0}};

double A[ii][jj] = {{ 3.0, -0.1, -0.2},
                    { 0.1,  7.0, -0.3},
                    { 0.3, -0.2, 10.0}};
double b[ii][kk] = {{ 7.85},
                    {-19.3},
                    { 71.4}};

double xo[jj][kk] = {{0.0},
                     {0.0}};
double xn[jj][kk] = {{0.0}};

double ro[jj][kk] = {{0.0}};
double rn[jj][kk] = {{0.0}};
double po[jj][kk] = {{0.0}};
double pn[jj][kk] = {{0.0}};

double tm[jj][kk] = {{0.0}};
double tn[jj][kk] = {{0.0}};
double to[jj][kk] = {{0.0}};
double tp[jj][kk] = {{0.0}};

int i = 0; int j = 0; int k = 0;

void mProduct(double A[ii][jj], double B[jj][kk])
{
    for (i = 0; i < ii; ++i)
        for (j = 0; j < jj; ++j)
            for (k = 0; k < kk; ++k)
            {
                C[i][k] += A[i][j]*B[j][k];
            }
    for (i = 0; i < ii; ++i)
    {
        for (k = 0; k < kk; ++k)
        {
            cout << C[i][k] << '\t';
        }
    cout << endl;
    }
}

int notmain()
{
//    mProduct(A, b);

    // Calculation: A*xo
    for (i = 0; i < ii; ++i)
        for (j = 0; j < jj; ++j)
            for (k = 0; k < kk; ++k)
                tm[i][k] += A[i][j]*xo[j][k];

    // Calculation: ro = b - A*xo
    for (j = 0; j < jj; ++j)
        for (k = 0; k < kk; ++k)
            ro[j][k] = b[j][k] - tm[j][k];

    // Calculation: po = ro for first iteration
    for (j = 0; j < jj; ++j)
        for (k = 0; k < kk; ++k)
            po[j][k] = ro[j][k];

    // Calculation: alphao_num = roT*ro (numerator of alpha_0)
    double alphao_num = 0.0;
    for (j = 0; j < jj; ++j)
        for (k = 0; k < kk; ++k)
            alphao_num += ro[k][j]*ro[j][k];

    // Calculation: tn = poT*A (for denominator of alpha_0)
    for (i = 0; i < ii; ++i)
        for (j = 0; j < jj; ++j)
            for (k = 0; k < kk; ++k)
                tn[k][i] += po[j][k]*A[i][j];

    // Calculation: alphao_den = tn*po
    // (denominator of alpha_0 = poT*A*po)
    double alphao_den = 0.0;
    for (j = 0; j < jj; ++j)
        for (k = 0; k < kk; ++k)
            alphao_den += tn[k][j]*po[j][k];

    // Calculation: alphao
    double alphao = alphao_num/alphao_den;

    for (int n = 1; n <= 1; ++n)
    {
        // Calculation: xn = xo + alphao*po
        for (j = 0; j < jj; ++j)
            for (k = 0; k < kk; ++k)
                xn[j][k] = xo[j][k] + alphao*po[j][k];

        // Calculation: rn = ro - alphao*A*po
        for (i = 0; i < ii; ++i)
            for (j = 0; j < jj; ++j)
                for (k = 0; k < kk; ++k)
                    to[i][k] += alphao*A[i][j]*po[j][k];
        for (j = 0; j < jj; ++j)
            for (k = 0; k < kk; ++k)
                rn[j][k] = ro[j][k] - to[j][k];

        double beta = 0.0;
        double rnTrn = 0.0;
        double roTro = 0.0;
        double pnTApn = 0.0;
        double alphan = 0.0;


        // Calculation: beta = (rnT*rn)/(roT*ro)
        for (j = 0; j < jj; ++j)
            for (k = 0; k < kk; ++k)
            {
                rnTrn += rn[k][j]*rn[j][k];
                roTro += ro[k][j]*ro[j][k];
            }
        beta = rnTrn/roTro;

        // Calculation: pn = rn + beta*po
        for (j = 0; j < jj; ++j)
            for (k = 0; k < kk; ++k)
                pn[j][k] = rn[j][k] + beta*po[j][k];

        // Calculation: alphan = (rnT*rn)/(pnT*A*pn)

        for (i = 0; i < ii; ++i)
            for (j = 0; j < jj; ++j)
                for (k = 0; k < kk; ++k)
                    tp[k][i] += pn[j][k]*A[i][j];

        for (j = 0; j < jj; ++j)
            for (k = 0; k < kk; ++k)
                pnTApn += tp[k][j]*pn[j][k];

        alphan = rnTrn/pnTApn;

        // Update value of xn
        for (j = 0; j < jj; ++j)
            for (k = 0; k < kk; ++k)
                xo[j][k] = xn[j][k];

        // Calculation: xn = xo + alphan*pn
        for (j = 0; j < jj; ++j)
            for (k = 0; k < kk; ++k)
                xn[j][k] = xo[j][k] + alphan*pn[j][k];
    }

    // Output
    cout << "xn =" << endl;
    for (j = 0; j < jj; ++j)
    {
        for (k = 0; k < kk; ++k)
            cout << xn[j][k] << '\t';
        cout << endl;
    }
    cout << endl;


    return 0;
}


//// Output: A*xo
//cout << "tm = " << endl;
//for (i = 0; i < ii; ++i)
//{
//    for (k = 0; k < kk; ++k)
//        cout << tm[i][k] << '\t';
//    cout << endl;
//}
//cout << endl;
//
//// Output: ro = b - A*xo
//cout << "ro = " << endl;
//for (i = 0; i < ii; ++i)
//{
//    for (k = 0; k < kk; ++k)
//        cout << ro[i][k] << '\t';
//    cout << endl;
//}
//cout << endl;
//
//// Output: po = ro for first iteration
//cout << "po = " << endl;
//for (i = 0; i < ii; ++i)
//{
//    for (k = 0; k < kk; ++k)
//        cout << po[i][k] << '\t';
//    cout << endl;
//}
//cout << endl;
//
//// Output: alphao_num = roT*ro (numerator of alpha_0)
//cout << "alphao_num =" << endl << alphao_num << endl;
//cout << endl;
//
//// Output: tn = poT*A (for denominator of alpha_0)
//cout << "tn =" << endl;
//for (k = 0; k < kk; ++k)
//{
//    for (i = 0; i < ii; ++i)
//        cout << tn[k][i] << '\t';
//    cout << endl;
//}
//cout << endl;
//
//// Output: alphao_num = roT*ro (numerator of alpha_0)
//cout << "alphao_den =" << endl << alphao_den << endl;
//cout << endl;
//
//// Output: alphao
//cout << "alphao =" << endl << alphao << endl;
//cout << endl;
//
//// Output: xn = xo + alphao*po
//cout << "xn = " << endl;
//for (j = 0; j < jj; ++j)
//{
//    for (k = 0; k < kk; ++k)
//        cout << xn[j][k] << '\t';
//    cout << endl;
//}
//cout << endl;
