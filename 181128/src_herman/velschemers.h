/*
 * velschemers.h
 *
 *     Project: Finite Volume Navier-Stokes Solver
 *      Author: Syed Ahmad Raza
 */

#ifndef VELSCHEMERS_H
#define VELSCHEMERS_H

void quick(int i, int j,\
    double& ue, double& uw, double& un, double& us,\
    double& vnu, double& vsu,\
    double& ve, double& vw, double& vn, double& vs,\
    double& uev, double& uwv);
// void upwind(int i, int j, int k);

#endif
