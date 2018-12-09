/*
 * filers.h
 *
 *  Created on: 2017-09-27
 *      Author: Syed Ahmad Raza
 */

#ifndef FILERS_H
#define FILERS_H

void filer1(const Numeric_lib::Matrix<double,1>& Z,
            const int nk, const std::string fileName);
void filer2(const Numeric_lib::Matrix<double,2>& Q,
            const int mx, const int my, const std::string fileName);
void filerAllSol(const int timeStep = -1);
void progressFileCreator();
void progressFiler();
void infoFilerStart();
void infoFilerEnd();

#endif
