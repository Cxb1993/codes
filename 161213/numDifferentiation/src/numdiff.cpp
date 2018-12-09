#include <iostream>
#include <fstream>  // for file manipulation
#include <iomanip>  // for std::precision
#include <cmath>
//#include <cstdlib>  // for abs function
#include <sstream>  // for integer to string function

using namespace std;

const double pi = 3.141592653589793;
const double fromLimit = 0.0;
const double toLimit = pi;

/*
 * COLUMNS FOR solFileGen
 * Define a variable for number of colsumns in the array with the following
 * values in each colsumn:
 * Column 0: Values of the nodes
 * Column 1: ANALYTICAL solutions at the respective node
 * Column 2: NUMERICAL solution at the node using FORWARD difference method
 * Column 3: Values of square of ratio of DIFFERENCE between FORWARD
 *  difference numerical and analytical solution
 * Column 4: NUMERICAL solution at the node using BACKWARD difference method
 * Column 5: Values of square of ratio of DIFFERENCE between BACKWARD
 *  difference numerical and analytical solution
 * Column 6: NUMERICAL solution at the node using CENTRAL difference method
 * Column 7: Values of square of ratio of DIFFERENCE between CENTRAL
 *  difference numerical and analytical solution
 */
const int cols = 8;

/*
 * COLUMNS FOR errFileGen
 * Define a variable for number of colsumns in the error array with the
 * the following values in each colsumn:
 * Column 0: Number of intervals
 * Column 1: Log of number of intervals
 * Column 2: Summation of colsumn 3 of the first array (sum of the errors
 *  of FORWARD difference method) for different intervals
 * Column 3: Error E of FORWARD difference method
 * Column 4: Log(E) of FORWARD difference method
 * Column 5: Summation of colsumn 5 of the first array (sum of the errors
 *  of BACKWARD difference method) for different intervals
 * Column 6: Error E of BACKWARD difference method
 * Column 7: Log(E) of BACKWARD difference method
 * Column 8: Summation of colsumn 7 of the first array (sum of the errors
 *  of CENTRAL difference method) for different intervals
 * Column 9: Error E of CENTRAL difference method
 * Column 10: Log(E) of CENTRAL difference method
 */
const int cols2 = 11;

double funcToDiff(double x) { return sin(3 * x); }
double anFirstSol(double x) { return (3 * cos(3 * x)); }
double anSecondSol(double x) { return (-9 * sin(3 * x)); }
//double funcToDiff(double x) { return cos(3 * x); }
//double anFirstSol(double x) { return -(3 * sin(3 * x)); }

/*
 * Following functions calculate the forward, backward and central difference
 * approximations of the first derivative where parameter 'i' is the value of
 * node and 'h' is the width of the interval.
 */
double firstFwd(double i, double h)
{
    return (funcToDiff(i + h) - funcToDiff(i)) / h;
}
double firstBwd(double i, double h)
{
    return (funcToDiff(i) - funcToDiff(i - h)) / h;
}
double firstCnt(double i, double h)
{
    return (funcToDiff(i + h) - funcToDiff(i - h)) / (2 * h);
}

/*
 * Following functions calculate the forward, backward and central difference
 * approximations of the second derivative where parameter 'i' is the value of
 * node and 'h' is the width of the interval.
 */
double secondFwd(double i, double h)
{
    return ((funcToDiff(i + 2 * h) - 2 * funcToDiff(i + h) + funcToDiff(i))
            / (h * h));
}
double secondBwd(double i, double h)
{
    return (funcToDiff(i) - 2 * funcToDiff(i - h) + funcToDiff(i - 2 * h))
            / (h * h);
}
double secondCnt(double i, double h)
{
    return (funcToDiff(i + h) - 2 * funcToDiff(i) + funcToDiff(i - h))
            / (h * h);
}

/*
 * When provided with the analytical solution, three numerical methods, the
 * array that is to be filled with values and the number of intervals, it
 * populates the array with the all the values.
 */
void arrayGen(double (*sol)(double),
        double (*fwd)(double, double),
        double (*bwd)(double, double),
        double (*cnt)(double, double),
        double table[][cols], int rows)
{
    double h = (toLimit - fromLimit) / static_cast<double>(rows);
                                        // width of interval
    double i = fromLimit + h;           // second node, ignoring the first

    for (int r = 1; r < rows; ++r)  // first row corresponds to first node
    {
        table[r][0] = i;            // index
        table[r][1] = sol(i);       // Aa
        table[r][2] = fwd(i, h);    // An for forward difference method
        table[r][3] = pow(abs((table[r][2] - table[r][1]) / table[r][1]), 2);
        table[r][4] = bwd(i, h);    // An for backward difference method
        table[r][5] = pow(abs((table[r][3] - table[r][1]) / table[r][1]), 2);
        table[r][6] = cnt(i, h);    // An for central difference method
        table[r][7] = pow(abs((table[r][4] - table[r][1]) / table[r][1]), 2);
        //            pow(abs((     An     -      Aa    ) /      Aa    ), 2)

        i += h;     // assigning the value of next node
    }
}

/* Generates file for first derivative using the given file name (and path),
 * and number of intervals
 */
void solFileGen(double (*sol)(double),
        double (*fwd)(double, double),
        double (*bwd)(double, double),
        double (*cnt)(double, double),
        string fileName, int rows)
{
    ofstream file; file << setprecision(50); file.open(fileName);

    double solValues[rows][cols] = { 0 };   // create the table array
    arrayGen(sol, fwd, bwd, cnt, solValues, rows);

    // File writing loop, ignoring the first and last rows of the array
    for (int r = 1; r < rows; ++r)
    {
        for (int c = 0; c < cols; ++c)
            file << solValues[r][c] << "\t";
        file << endl;
    }
    file.close();
}

void errFileGen(double (*sol)(double),
        double (*fwd)(double, double),
        double (*bwd)(double, double),
        double (*cnt)(double, double),
        string fileName, int intervalsStart, int intervalsFinal)
{
    ofstream file; file << setprecision(50); file.open(fileName);

    int rows2 = (intervalsFinal - intervalsStart) + 1; //rows of error array
    double errValues[rows2][cols2] = { 0 };
    int r = 0;  // index for rows of errValues

    // rows == intervals == nodes - 1; because starting index from 0, there
    // will always be one extra row, representing one extra node.
    for (int rows = intervalsStart; rows <= intervalsFinal; ++rows)
    {
        double solValues[rows][cols] = { 0 }; // table array
        arrayGen(sol, fwd, bwd, cnt, solValues, rows);

        double summationFwd = 0.0;
        double summationBwd = 0.0;
        double summationCnt = 0.0;

        // Indexing each row of solValues to find summation using the index
        // variable 'rr'. This block could have been included in arrayGen.
        for (int rr = 1; rr < rows; ++rr)
        {
            summationFwd += solValues[rr][3];
            summationBwd += solValues[rr][5];
            summationCnt += solValues[rr][7];
        }
        errValues[r][0] = r;
        errValues[r][1] = log10(errValues[r][0]);
        errValues[r][2] = summationFwd;
        errValues[r][3] = pow(errValues[r][2] / errValues[r][0], 0.5);
        errValues[r][4] = log10(errValues[r][3]);
        errValues[r][5] = summationBwd;
        errValues[r][6] = pow(errValues[r][5] / errValues[r][0], 0.5);
        errValues[r][7] = log10(errValues[r][6]);
        errValues[r][8] = summationCnt;
        errValues[r][9] = pow(errValues[r][5] / errValues[r][0], 0.5);
        errValues[r][10] = log10(errValues[r][9]);

        for (int c = 0; c < cols2; ++c)
            file << errValues[r][c] << "\t";
        file << endl;
        ++r;    // increasing the index for error array
        // Note that variable 'rows' corresponds to the number of intervals
        // for the respective pass whereas, variable 'r' corresponds to the
        // particular row of errValues being filled.
    }
    file.close();
}

int main()
{
    cout << setprecision(50);

    // PART 1
    // Generate files for plotting numerical solution vs analytical solution
    int intervalsStart = 10;   // number of intervals, initally
    int intervalsFinal = 100;  // number of intervals, finally

    // rows == intervals == nodes - 1; because starting index from 0, there
    // will always be one extra rows, representing one extra node.
    for (int rows = intervalsStart; rows <= intervalsFinal; rows *= 10)
    {
        ostringstream convert; convert << rows;  // for int-to-string conv.
        solFileGen(anFirstSol, firstFwd, firstBwd, firstCnt,
                "./data/firstsin" + convert.str() + ".dat", rows);
    }

    for (int rows = intervalsStart; rows <= intervalsFinal; rows *= 10)
        {
            ostringstream convert; convert << rows;  // for int-to-string conv.
            solFileGen(anSecondSol, secondFwd, secondBwd, secondCnt,
                    "./data/secondsin" + convert.str() + ".dat", rows);
        }

    // PART 2
    // Generate file with the error
    intervalsFinal = 10000;  // number of intervals, finally
    errFileGen(anFirstSol, firstFwd, firstBwd, firstCnt,
            "./data/firstsinerror.dat", intervalsStart, intervalsFinal);
    errFileGen(anSecondSol, secondFwd, secondBwd, secondCnt,
                "./data/secondsinerror.dat", intervalsStart, intervalsFinal);

    return 0;
}

/*
 * FILE GENERATOR (EXTRA FUNCTION)
 * Generates a file with the specified filename, using the given 2D array
 * with the specified number of colsumns and rowss, except for the first and the
 * last rows. All the colsumns of the array are separated using tab character
 * and the rowss are separated using newline character. The first rows is omitted
 * because in this code they generally contain zeros.
 */
//void fileGen(string fileName, double table[][cols], int rows)
//{
//    ofstream file;
//    file << setprecision(50);
//    file.open(fileName);
//
//    // File writing loop, ignoring the first and last rowss
//    for (int r = 1; r < rows; ++r)
//    {
//        for (int c = 0; c < cols; ++c)
//        {
//            file << table[r][c] << "\t";
//        }
//        file << endl;
//    }
//
//    file.close();
//}
