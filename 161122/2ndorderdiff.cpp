#include <iostream>
#include <fstream>          //for file manipulation
#include <iomanip>          //for std::setprecision()
#include <cmath>
#include <string>

const long double pi(3.141592653589793238462643383279);
long double sine(long double x);
long double cosine(long double x);
void fwd2(long double (*f)(long double), long double a,
        long double b, int intervals, std::string file);
void bwd2(long double (*f)(long double), long double a,
        long double b, int intervals, std::string file);
void cnt2(long double (*f)(long double), long double a,
        long double b, int intervals, std::string file);

int main()
{
    using namespace std;

    // Generating files for differentiation of a sine function
    for (int m = 10; m <= 100; m *= 10)
    { 
        fwd2(sine, 0, pi, m, "2ndsinfwd" + to_string(m) + ".dat");
        bwd2(sine, 0, pi, m, "2ndsinbwd" + to_string(m) + ".dat");
        cnt2(sine, 0, pi, m, "2ndsincnt" + to_string(m) + ".dat");
    }

    // Generating files for differentiation of a cosine function
    for (int m = 10; m <= 100; m *= 10)
    { 
        fwd2(cosine, 0, pi, m, "2ndcosfwd" + to_string(m) + ".dat");
        bwd2(cosine, 0, pi, m, "2ndcosbwd" + to_string(m) + ".dat");
        cnt2(cosine, 0, pi, m, "2ndcoscnt" + to_string(m) + ".dat");
    }

    // Generating error files for sine function differentiation

    return 0;
}

long double sine(long double x)
{
    return sin(3*x);          //function to be differentiated
}

long double cosine(long double x)
{
    return cos(3*x);          //function to be differentiated
}
    

// 2nd order forward difference method

void fwd2(long double (*f)(long double), long double a,
        long double b, int intervals, std::string file)
{
    using namespace std;

    ofstream fwd2_file;
    fwd2_file << setprecision(100);
    fwd2_file.open(file);

    long double h = (b - a) / intervals;
    long double i = a;

    for (int n = 0; n < (intervals - 1); ++n)
    {
        long double An = (f(i + 2 * h) - 2 * f(i + h) + f(i)) / (h * h);
        fwd2_file << i << "  " << An << endl;
        i += h;
    }
    fwd2_file.close();
}

// 2nd order backward difference method

void bwd2(long double (*f)(long double), long double a,
        long double b, int intervals, std::string file)
{
    using namespace std;

    ofstream bwd2_file;
    bwd2_file << setprecision(100);
    bwd2_file.open(file);

    long double h = (b - a) / intervals;
    long double i = a + h;

    for (int n = 2; n <= intervals; ++n)
    {
        long double An = (f(i) - 2 * f(i - h) + f(i - 2 * h)) / (h * h);
        bwd2_file << i << "  " << An << endl;
        i += h;
    }
    bwd2_file.close();
}

// 2nd order central difference method

void cnt2(long double (*f)(long double), long double a,
        long double b, int intervals, std::string file)
{
    using namespace std;

    ofstream cnt2_file;
    cnt2_file << setprecision(100);
    cnt2_file.open(file);

    long double h = (b - a) / intervals;
    long double i = a + h;

    for (int n = 2; n < (intervals - 1); ++n)
    {
        long double An = (f(i + h) - 2 * f(i) + f(i - h)) / (h * h);
        cnt2_file << i << "  " << An << endl;
        i += h;
    }
    cnt2_file.close();
}


