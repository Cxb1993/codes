"""Step 3 of the online CFD course of Lorena Barbara.

Practice the step 3 of online CFD course of Lorena Barbara called 12 steps to
Navier-Stokes equations:
1D Diffusion
step3.py

Created on: 2017-12-23
Author: Syed Ahmad Raza
"""

import numpy as np
import matplotlib.pyplot as plt
# import time
# import sys

nx = 41             # number of grid points
dx = 2 / (nx-1)     # distance between any pair of adjacent grid points
nt = 20             # nt is the number of timesteps we want to calculate
nu = 0.3            # the value of viscosity
sigma = 0.2         # maximum Courant number
dt = sigma * dx**2 / nu     # timestep defined using Courant number

u = np.ones(nx)     # numpy function ones()
# set u = 2 between 0.5 and 1 as per our initial conditions
u[int(.5 / dx):int(1 / dx + 1)] = 2

plt.plot(np.linspace(0, 2, nx), u)

un = np.ones(nx)    # initialize a temporary array

for n in range(nt):     # loop for values of n from 0 to nt; run nt times
    un = u.copy()       # copy the existing values of u into un
    for i in range(1, nx-1):
        u[i] = un[i] + nu * dt / dx**2 * (un[i+1] - 2 * un[i] + un[i-1])

plt.plot(np.linspace(0, 2, nx), u)
plt.show()
# plt.savefig("step3.eps")
