/*
 * filer.h
 *
 *  Created on: 2017-04-20
 *      Author: Syed Ahmad Raza
 */

#ifndef FILERS_H
#define FILERS_H

#include <string>       // for file name

using namespace std;

void gridFiler(double Z[], int nk, string fileName);
void solFiler(double Q[][ny], int nx, int ny, string fileName);

#endif
