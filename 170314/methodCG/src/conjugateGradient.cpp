/*
 * conjugateGradient.cpp
 *
 *  Created on: Mar 13, 2017
 *      Author: Syed Ahmad Raza
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

using namespace std;

const int ii = 2;
const int jj = ii;
const int kk = 1;

double C[ii][kk] = {{0.0}};

double A[ii][jj] = {{4.0, 1.0},
                    {1.0, 3.0}};
double b[ii][kk] = {{1.0},
                    {2.0}};

double xo[jj][kk] = {{0.0},
                     {0.0}};
double xn[jj][kk] = {{0.0}};

double ro[ii][kk] = {{0.0}};
double rn[ii][kk] = {{0.0}};
double po[ii][kk] = {{0.0}};
double pn[ii][kk] = {{0.0}};

double tm[ii][kk] = {{0.0}};
double tn[ii][kk] = {{0.0}};
double to[ii][kk] = {{0.0}};
double tp[ii][kk] = {{0.0}};

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

int main()
{
//    mProduct(A, b);

    // Calculation: A*xo
    for (i = 0; i < ii; ++i)
        for (j = 0; j < jj; ++j)
            for (k = 0; k < kk; ++k)
                tm[i][k] += A[i][j]*xo[j][k];

    // Calculation: ro = b - A*xo
    for (i = 0; i < ii; ++i)
        for (k = 0; k < kk; ++k)
            ro[i][k] = b[i][k] - tm[i][k];

    // Calculation: po = ro for first iteration
    for (i = 0; i < ii; ++i)
        for (k = 0; k < kk; ++k)
            po[i][k] = ro[i][k];

    // Calculation: alphao_num = roT*ro (numerator of alpha_0)
    double alphao_num = 0.0;
    for (i = 0; i < ii; ++i)
        for (k = 0; k < kk; ++k)
            alphao_num += ro[k][i]*ro[i][k];

    // Calculation: tn = poT*A (for denominator of alpha_0)
    for (i = 0; i < ii; ++i)
        for (j = 0; j < jj; ++j)
            for (k = 0; k < kk; ++k)
                tn[k][i] += po[j][k]*A[i][j];

    // Calculation: alphao_den = tn*po
    // (denominator of alpha_0 = poT*A*po)
    double alphao_den = 0.0;
    for (i = 0; i < ii; ++i)
        for (k = 0; k < kk; ++k)
            alphao_den += tn[k][i]*po[i][k];

    // Calculation: alphao
    double alphao = alphao_num/alphao_den;

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

    double betao = 0.0;
    double rnTrn = 0.0;
    double roTro = 0.0;
    double pnTApn = 0.0;
    double alphan = 0.0;
//    double betan = 0.0;

    // Calculation: betao = (rnT*rn)/(roT*ro)
    for (j = 0; j < jj; ++j)
        for (k = 0; k < kk; ++k)
        {
            rnTrn += rn[k][j]*rn[j][k];
            roTro += ro[k][j]*ro[j][k];
        }
    betao = rnTrn/roTro;

    // Calculation: pn = rn + betao*po
    for (j = 0; j < jj; ++j)
        for (k = 0; k < kk; ++k)
            pn[j][k] = rn[j][k] + betao*po[j][k];

    // Calculation: alphan = (rnT*rn)/(pnT*A*pn)

    for (i = 0; i < ii; ++i)
        for (j = 0; j < jj; ++j)
            for (k = 0; k < kk; ++k)
                tp[k][i] += pn[j][k]*A[i][j];

    for (i = 0; i < ii; ++i)
        for (k = 0; k < kk; ++k)
            pnTApn += tp[k][i]*pn[i][k];

    alphan = rnTrn/pnTApn;

    // Update value of xn
    for (j = 0; j < jj; ++j)
        for (k = 0; k < kk; ++k)
            xo[j][k] = xn[j][k];

    // Calculation: xn = xo + alpha1*p1
    for (j = 0; j < jj; ++j)
        for (k = 0; k < kk; ++k)
            xn[j][k] = xo[j][k] + alphan*pn[j][k];

    // Output
    cout << "xn =" << endl;
    for (i = 0; i < ii; ++i)
    {
        for (k = 0; k < kk; ++k)
            cout << xn[i][k] << '\t';
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
