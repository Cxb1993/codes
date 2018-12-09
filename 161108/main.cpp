#include <iostream>
#include <fstream>
#include <iomanip>          //for std::setprecision()
#include <cmath>
#include "trapezoidal.h"
#include "simpsons13.h"

const long double pi(3.141592653589793238462643383279);

long double f(long double x)
{
    return sin(x);          //function to be integrated
}

long double integration(long double (*integFunc)(
    long double (*)(long double), int, long double, long double), int n)
    {
        return integFunc(f, n, 0.0, pi);
    }

long double diff(long double intgrResult)
{
    //For the case of the integral of sin(x) between the limits
    //0 and pi, the analytical solution is 2
    return (abs(2 - intgrResult))/2;
}

int main()
{
    using namespace std;
    
    cout << setprecision(100);
    
    //Introduction
    cout << "This program will perform integration between the limits of "
            << "0 and pi for sin(x)" << endl << endl;
    
    //Question 1
    cout << "Which method would you like to use?" << endl <<
            " ('t' for trapezoidal and 's' for simpsons13)" << endl <<
            "Input: ";
    string userInput1;
    cin >> userInput1;
    cout << endl;
    
    //Question 2
    cout << "What type of result do you want? (1/2/3)" << endl <<
        "1. Dispaly numerical results" << endl <<
        "2. Display difference between the analytical and numerical " <<
        "results" << endl << "3. Write the data (log of no. of segments" <<
        "  VS log of difference between results) to a file" << endl <<
        "Input: ";
        
    int userInput2;
    cin >> userInput2;
    
    //Calculating and displaying results
    cout << endl << "The results are:" << endl << endl;
    
    if (userInput1 == "t" && userInput2 == 1)
        for (int n = 10; n <= 1000; ++n)
        {
            cout << n << "    " << integration(trapezoidal, n) << endl;
        }
    else if (userInput1 == "t" && userInput2 == 2)
        for (int n = 10; n <= 1000; ++n)
        {
            long double intgr = integration(trapezoidal, n);
            //For the case of the integral of sin(x) between the limits
            //0 and pi, the analytical solution is 2
            cout << n << "    " << (abs(2 - intgr))/2 << endl;
        }
    else if (userInput1 == "s" && userInput2 == 1)
        for (int n = 10; n <= 1000; n += 2)
        {
            cout << n << "    " << integration(simpsons13, n) << endl;
        }
    else if (userInput1 == "s" && userInput2 == 2)
        for (int n = 10; n <= 1000; n += 2)
        {
            long double intgr = integration(simpsons13, n);
            //For the case of the integral of sin(x) between the limits
            //0 and pi, the analytical solution is 2
            cout << n << "    " << (abs(2 - intgr))/2 << endl;
        }
    else if(userInput1 == "t" && userInput2 == 3)
    {
        ofstream out_file;
        out_file << setprecision(100);
        out_file.open("trapezoidal.txt");
        for (int n = 10; n <= 1000; ++n)
        {
            long double intgr = integration(trapezoidal, n);
            //For the case of the integral of sin(x) between the limits
            //0 and pi, the analytical solution is 2
            long double N = log10(n);
            long double D = log10((abs(2 - intgr))/2);
            out_file << N << " " << D << endl;
        }
        cout << "... written to the file.";
    }
    else if (userInput1 == "s" && userInput2 == 3)
    {
        ofstream out_file;
        out_file << setprecision(100);
        out_file.open("simpsons13.txt");
        for (int n = 10; n <= 1000; n += 2)
        {
            long double intgr = integration(simpsons13, n);
            //For the case of the integral of sin(x) between the limits
            //0 and pi, the analytical solution is 2
            long double N = log10(n);
            long double D = log10((abs(2 - intgr))/2);
            out_file << N << " " << D << endl;
        }
        cout << "... written to the file.";
    }
    else    cout << "Your input was invalid." << endl;
    
    return 0;
}