// Example program
#include <iostream>
#include <iomanip>
// #include <string>
#include <limits>
#include "Matrix11.h"
#include "MatrixIO11.h"

using namespace std;
using namespace Numeric_lib;


void f(int n1, int n2, int n3)
{
    // Matrix<int,0> ai0; // error: no 0D matrices

    Matrix<double,1> ad1(5);
    Matrix<int,1> ai(5);
    Matrix<double,1> ad11(7);

    Matrix<double,3> ad3(n1,n2,n3);
    Matrix<double,3> ad33(n1,n2,n3);

    ad3 = ad33;

    cout << ad3.size() << endl;
    cout << ad3.dim1() << endl;
    cout << ad3.dim2() << endl;
    cout << ad3.dim3() << endl;

    // Matrix<double,1>   ad1(n1);
    // Matrix<int,1>      ai1(n1);
    //
    // ad1(7) = 0;
    // ad1[7] = 8;
    //
    // Matrix<double,2>    ad2(n1, n2);
    // Matrix<double,3>    ad3(n1, n2, n3);
    //
    // ad2(3, 4) = 7.5;
    // ad3(3, 4, 5) = 9.2;
    // ad3[6][0][9] = 9.6;
    //
    // cout << ad2[3][4]       << endl;
    // cout << ad3(3, 4, 5)    << endl;
    // cout << ad3[3][4][5]    << endl;
    // cout << ad3[6][9][9]    << endl;
    // cout << ad3[6][0][9]    << endl;
}


double scale(double d, double s)
{
    return d*s;
}


void scaleInPlace(double& d, double s)
{
    d *= s;
}


void init2(Matrix<double,2>& a)
{
    for (int i = 0; i < a.dim1(); ++i)
    {
        for (int j = 0; j < a.dim2(); ++j)
        {
            a[i][j] = 10*i+j;
        }
    }
}


void print2(const Matrix<double,2>& a)
{
    for (int i = 0; i < a.dim1(); ++i)
    {
        for (int j = 0; j < a.dim2(); ++j)
        {
            cout << a(i,j) << '\t';
        }
        cout << '\n';
    }
}


void init1(Matrix<double,1>& a)
{
    for (int i = 0; i < a.dim1(); ++i)
    {
        a[i] = i;
    }
}


void print1(const Matrix<double,1>& a)
{
    for (int i = 0; i < a.dim1(); ++i)
    {
        cout << a(i) << '\t';
    }
    cout << '\n';
}


void copyBuiltIn(double* p, int n)
{
    double val[] = {1.2, 2.3, 3.4, 4.5};
    Matrix<double> data(p, n);
    Matrix<double> constants(val);

    print1(data);
    print1(constants);
}

// int main()
// {
//     Matrix<double,2> a(9, 10);
//
//     init2(a);
//
//     Matrix<double,2> b = apply(scale,a,7);
//
//     a.apply(scaleInPlace, 7);
//
//     double c[5] = {1.1, 2.2, 3.3, 4.4, 5.5};
//     copyBuiltIn(c, 5);
//
//     Matrix<double,1> d(10);
//     init1(d);
//     print1(d);
//
//     print1(a[2].slice(4));
//
//     Matrix<double,1> e = a[2].slice(5,2);
//     print1(e);
//
//     Matrix<double,3> f(3, 3, 3);
//     cout << f(2,2,2);
//
//     // cout << d3[4][3];
//     // cout << d3[4].slice(4);
//     // f(8, 9, 10);
//
//     return 0;
// }

int main()
{
    Matrix<int,2> m(2,2);
    cin >> m;
    cout << m;
    return 0;
}
