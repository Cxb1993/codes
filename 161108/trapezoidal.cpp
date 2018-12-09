// Trapezoidal rule for integration with limits

#include <iostream>
#include <cmath>

long double trapezoidal(long double (*functocall)(long double),
    int n, long double a, long double b)
{
    long double segment = (b - a) / n;                  //width of one segment
    long double summation = functocall(a);              //first point
    long double j = a;
    
    for (int i = 1; i < n; ++i)
    {
        j += segment;
        summation += 2 * functocall(j);                 //middle points
    }
    
    summation += functocall(b);                         //last point
    
    long double answer = (b - a) * summation / (2 * n); //final formula
    
    return answer;
}