#include <iostream>
#include <fstream>          //for file manipulation
#include <iomanip>          //for std::setprecision()
#include <cmath>
#include <string>

const long double pi(3.141592653589793238462643383279);
long double sine(long double x);
long double cosine(long double x);
void fwd(long double (*f)(long double), long double a,
        long double b, int intervals, std::string file);
void bwd(long double (*f)(long double), long double a,
        long double b, int intervals, std::string file);
void cnt(long double (*f)(long double), long double a,
        long double b, int intervals, std::string file);

int main()
{
    using namespace std;

    // Generating files for differentiation of a sine function
    for (int m = 10; m <= 100; m *= 10)
    { 
        fwd(sine, 0, pi, m, "1stsinfwd" + to_string(m) + ".dat");
        bwd(sine, 0, pi, m, "1stsinbwd" + to_string(m) + ".dat");
        cnt(sine, 0, pi, m, "1stsincnt" + to_string(m) + ".dat");
    }

    // Generating files for differentiation of a cosine function
    for (int m = 10; m <= 100; m *= 10)
    { 
        fwd(cosine, 0, pi, m, "1stcosfwd" + to_string(m) + ".dat");
        bwd(cosine, 0, pi, m, "1stcosbwd" + to_string(m) + ".dat");
        cnt(cosine, 0, pi, m, "1stcoscnt" + to_string(m) + ".dat");
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
    
// Forward difference method

void fwd(long double (*f)(long double), long double a,
        long double b, int intervals, std::string file)
{
    using namespace std;

    ofstream fwd_file;
    fwd_file << setprecision(100);
    fwd_file.open(file);

    long double h = (b - a) / intervals;
    long double i = a;

    for (int n = 0; n < intervals; ++n)
    {
        long double An = (f(i + h) - f(i)) / h;
        fwd_file << i << "  " << An << endl;
        i += h;
    }
    fwd_file.close();
}

// Backward difference method

void bwd(long double (*f)(long double), long double a,
        long double b, int intervals, std::string file)
{
    using namespace std;

    ofstream bwd_file;
    bwd_file << setprecision(100);
    bwd_file.open(file);

    long double h = (b - a) / intervals;
    long double i = a + h;

    for (int n = 1; n <= intervals; ++n)
    {
        long double An = (f(i) - f(i - h)) / h;
        bwd_file << i << "  " << An << endl;
        i += h;
    }
    bwd_file.close();
}

// Central difference method

void cnt(long double (*f)(long double), long double a,
        long double b, int intervals, std::string file)
{
    using namespace std;

    ofstream cnt_file;
    cnt_file << setprecision(100);
    cnt_file.open(file);

    long double h = (b - a) / intervals;
    long double i = a + h;

    for (int n = 1; n < intervals; ++n)
    {
        long double An = (f(i + h) - f(i - h)) / (2 * h);
        cnt_file << i << "  " << An << endl;
        i += h;
    }
    cnt_file.close();
}

