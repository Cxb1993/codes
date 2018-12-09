/*
 * navierFVD.h
 *
 *  Created on: 2017-09-27
 *      Author: Syed Ahmad Raza
 */

#ifndef NAVIERFVD_H
#define NAVIERFVD_H

// Solver
void initial();
void velBoundary();
void pressBoundary();
void momentum(string velScheme);
void momentumTimeStep();
void pressure();
void velUpdater();
void navierFVD();

#endif
