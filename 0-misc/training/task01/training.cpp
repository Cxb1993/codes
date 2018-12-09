/*
 * training.cpp
 *
 *  Created on: 2017-11-01
 *      Author: Syed Ahmad Raza
 */

#include <iostream>     // functions for input and output to console
#include <fstream>      // functions for file input and output
#include <sstream>      // functions for string to number conversion

using namespace std;

int main()
{
    // CREATING ARRAYS

    // Creating a single-dimensional array with all elements initialized to zero
    double A[5] = {0.0};

    // Creating a two-dimensional array with all elements initialized to zero
    double B[5][5] = { {0.0} };
    // Notice that we must specify and fix all the dimensions of the array
    // beforehand; if we are specifying all the elements of the array, then
    // we can omit the first dimension on the left, e.g. B[][5]. But this is
    // only possible when all elements are specified, unlike here.

    // Creating a single-dimensional array with some arbitrary values
    double C[] = {2, 2, 3, 4, 4};

    // Check the output to see if it is as expected
    cout << "A[2]\t= " << A[2] << "\n";
    cout << "B[0][1]\t= " << B[0][1] << "\n";
    cout << "B[2][1]\t= " << B[2][1] << " => This value was not defined by us.\
    \n\tTherefore, the computer outputs a random value stored in that location."
    << "\n";
    cout << "C[2]\t= " << C[2] << "\n";

    // INPUT FILE OPERATIONS

    // Define an input file variable
    ifstream inputFileA;
    // Open the file for reading data and assign it to input file variable
    inputFileA.open("inputA.txt");
    // Define a string variable to store data from a single row
    string rowA;

    int i = 0;  // for indexing the loop

    // Define a loop that gets one line of data from "inputFile" variable and
    // places it in string variable "row", and loops until there is any row of
    // data left in the "inputFile" variable
    while (getline(inputFileA, rowA))
    {
        // Copy the string variable "row" of data to a stringstream variable
        // "ssA"
        stringstream ssA(rowA);
        // Convert and copy the stringstream variable "ssA" to an element in the
        // array variable "A[]"
        ssA >> A[i]; // note that it only transfers the first item of the row
        ++i;    // increase the index for the next iteration
    }

    // Use similar procedures for two-dimensional array "B[][5]"; observe the
    // differences
    ifstream inputFileB;
    inputFileB.open("inputB.txt");
    string rowB;
    int j = 0;
    while (getline(inputFileB, rowB))
    {
        stringstream ssB(rowB);
        ssB >> B[0][j]; // first element is copied
        ssB >> B[1][j]; // second element is copied
        ++j;
    }

    // ARRAY OPERATIONS
    for (int k = 0; k < 5; ++k)
    {
        C[k] = A[k] * B[0][k];  // C array values are overwritten
    }

    // OUTPUT FILE OPERATIONS

    // Define an output file variables
    ofstream outputFileC;
    // Open the file for writing data and assign it the output file variables
    outputFileC.open("outputC.txt");

    // Define a loop that writes every line of file using the array variable C[]
    for (int k = 0; k < 5; ++k)
    {
        // Write a single line of the fileQ
        outputFileC << C[k] << "\n";
    }
}
