#ifndef TRAPEZOIDAL_H
#define TRAPEZOIDAL_H

long double trapezoidal(long double (*functocall)(long double),
    int n, long double a, long double b);

#endif