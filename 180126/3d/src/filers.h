/*
 * filers.h
 *
 *  Created on: 2018-01-26
 *      Author: Syed Ahmad Raza
 */

#ifndef FILERS_H
#define FILERS_H

void filer1(double Z[], int nk, std::string fileName);
void filer2(double Q[][ny], int nx, int ny, std::string fileName);
void filer3(double S[][my][mz], int mx, int my, int mz, std::string fileName);
void filerAllSol(int timeStep = -1);
void progressFileCreator();
void progressFiler();
void infoFiler();
void temp();

#endif
