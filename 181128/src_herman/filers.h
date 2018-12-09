/*
 * filers.h
 *
 *     Project: Finite Volume Navier-Stokes Solver
 *      Author: Syed Ahmad Raza
 */

#ifndef FILERS_H
#define FILERS_H

void filerAllSol(const int timeStep = -1);
void filerCoordinates();
void filerInfoStart();
void filerCreateCD();
void filerCD();
void filerCreateCL();
void filerCL();
void filerInfoEnd();

#endif
