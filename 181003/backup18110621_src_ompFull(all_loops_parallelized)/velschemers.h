/*
 * velschemers.h
 *
 *     Project: Finite Volume Navier-Stokes Solver
 *      Author: Syed Ahmad Raza
 */

#ifndef VELSCHEMERS_H
#define VELSCHEMERS_H

void quick(int i, int j, int k,\
    double& ue, double& uw, double& un, double& us, double& uf, double& ub,\
    double& vnu, double& vsu, double& wfu, double& wbu,\
    double& ve, double& vw, double& vn, double& vs, double& vf, double& vb,\
    double& uev, double& uwv, double& wfv, double& wbv,\
    double& we, double& ww, double& wn, double& ws, double& wf, double& wb,\
    double& uew, double& uww, double& vnw, double& vsw);
// void upwind(int i, int j, int k);

#endif
