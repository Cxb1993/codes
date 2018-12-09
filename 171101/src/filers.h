/*
 * filers.h
 *
 *  Created on: 2017-09-27
 *      Author: Syed Ahmad Raza
 */

#ifndef FILERS_H
#define FILERS_H

#include <string>       // for fileName

void filer1(double Z[], int nk, std::string fileName);
void filer2(double Q[][ny], int nx, int ny, std::string fileName);
void filerAllSol(int timeStep = -1);
void progressPrinter();
void infoWriter();

#endif
