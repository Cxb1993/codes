/*
 * navierFVD.h
 *
 *  Created on: 2017-06-13
 *      Author: Syed Ahmad Raza
 */

#ifndef NAVIERFVD_H
#define NAVIERFVD_H

// Solver
void initial();
void velBoundary();
void pressBoundary();
void momentum();
void momentumTimeStep();
void pressure();
void velUpdater();
void navierFVD();

#endif
