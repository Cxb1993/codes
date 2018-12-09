// Simpson's 1/3 rule for integration with limits

#include <iostream>
#include <cmath>

long double simpsons13(long double (*functocall)(long double),
    int n, long double a, long double b)
{
    long double segment = (b - a) / n;                  //width of one segment
    long double summation = functocall(a);              //first point
    long double j = a;
    
    for (int i = 2; i <= n - 1; i += 2)
    {
        j += segment;
        summation += 4 * functocall(j);                 //odd points
        j += segment;
        summation += 2 * functocall(j);                 //even points
    }
    
    j += segment;
    summation += 4 * functocall(j) + functocall(b);     //last two points
    
    long double answer = (b - a) * summation / (3 * n); //final forumla
    
    return answer;
}