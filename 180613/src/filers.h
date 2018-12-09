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
void filerCreateProgress();
void filerProgress();
void filerInfoStart();
void filerInfoEnd();

#endif
