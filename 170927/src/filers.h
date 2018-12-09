/*
 * filer.h
 *
 *  Created on: 2017-09-27
 *      Author: Syed Ahmad Raza
 */

#ifndef FILERS_H
#define FILERS_H

#include <string>       // for file name

using namespace std;

void filer1(double Z[], int nk, string fileName);
void filer2(double Q[][ny], int nx, int ny, string fileName);
void filerAllSol(string fileNamePrefix, int timeStep = -1);

#endif
