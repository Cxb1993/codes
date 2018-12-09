/*
 * sorSolver.cpp
 *
 *  Created on: Mar 6, 2017
 *      Author: Syed Ahmad Raza
 */

#include <iostream>
#include <fstream>
#include <sstream>  // for int to string conversion
#include <cmath>

using namespace std;

double a11 = 3.0;
double a12 = -0.1;
double a13 = -0.2;
double a21 = 0.1;
double a22 = 7.0;
double a23 = -0.3;
double a31 = 0.3;
double a32 = -0.2;
double a33 = 10.0;
double b1 = 7.85;
double b2 = -19.3;
double b3 = 71.4;

double omega = 1.5;

double x_1(double x2v, double x3v)
{
    return (b1 - a12 * x2v - a13 * x3v) / a11;
}

double x_2(double x1v, double x3v)
{
    return (b2 - a21 * x1v - a23 * x3v) / a22;
}

double x_3(double x1v, double x2v)
{
    return (b3 - a31 * x1v - a32 * x2v) / a33;
}

double error(double xCurrent, double xLast)
{
    return abs((xCurrent - xLast) / xCurrent) * 100;
}

//int main()
//{
//    double x1 = x_1(0.0, 0.0);
//    double x2 = x_2(x1, 0.0);
//    double x3 = x_3(x1, x2);
//
//    int i = 1;
//    cout << i << '\t' << x1 << '\t' << x2 << '\t' << x3 << '\n';
//
//    double x1o = x1;
//    double x2o = x2;
//    double x3o = x3;
//
//    x1 = x_1(x2o, x3o) * omega + (1 - omega) * x1o;
//    x2 = x_2(x1, x3o) * omega + (1 - omega) * x2o;
//    x3 = x_3(x1, x2) * omega + (1 - omega) * x3o;
//
//    ++i;
//    cout << i << '\t' << x1 << '\t' << x2 << '\t' << x3 << '\n';
//
//    double errorValue1 = 100.0;
//    double errorValue2 = 100.0;
//    double errorValue3 = 100.0;
//
//    while (errorValue1 >= 0.01 || errorValue2 >= 0.01 || errorValue3 >= 0.01)
//    {
//        x1 = x_1(x2o, x3o) * omega + (1 - omega) * x1o;
//        x2 = x_2(x1, x3o) * omega + (1 - omega) * x2o;
//        x3 = x_3(x1, x2) * omega + (1 - omega) * x3o;
//
//        errorValue1 = error(x1, x1o);
//        errorValue2 = error(x2, x2o);
//        errorValue3 = error(x3, x3o);
//
//        x1o = x1;
//        x2o = x2;
//        x3o = x3;
//
//        ++i;
//        cout << i << '\t' << x1 << '\t' << x2 << '\t' << x3 << '\n';
//    }
//
//    return 0;
//}
