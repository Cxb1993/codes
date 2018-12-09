/*
 * navierFVD.h
 *
 *  Created on: 2017-04-24
 *      Author: Syed Ahmad Raza
 */

#ifndef NAVIERFVD_H
#define NAVIERFVD_H

// Solver
void initial();
void boundary();
void momentum();
void momentumTimeStep();
void pressure();
void velocityCorrector();
void temp();

#endif
