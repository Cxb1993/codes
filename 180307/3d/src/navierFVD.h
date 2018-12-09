/*
 * navierFVD.h
 *
 *  Created on: 2018-01-26
 *      Author: Syed Ahmad Raza
 */

#ifndef NAVIERFVD_H
#define NAVIERFVD_H

// Solver
void initial();
void velBoundary();
void pressBoundary();
void momentum(std::string velScheme);
void momentumTimeStep();
void pressure();
void velUpdater();
void navierFVD();

#endif
