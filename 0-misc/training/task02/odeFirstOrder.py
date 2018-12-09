"""Numerical solution of ordinary differential equations.

Numerical solution of two simple ordinary differential equations. It was
written to test the task for training practice.
odeFirstOrder.py

Created on: 2017-11-10
Author: Syed Ahmad Raza
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams['image.cmap'] = 'gist_rainbow'

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# x = np.loadtxt('x.dat')
# y = np.loadtxt('y.dat')

maxTime = 2   # seconds

# CASE 1
tSteps1 = 21

t = np.linspace(0, maxTime, tSteps1)
y = np.exp(-t)

dydt = np.zeros(tSteps1 - 1)
for i in range(tSteps1 - 1):
    dydt[i] = (y[i+1] - y[i]) / (t[i+1] - t[i])

# Figure 1
plt.figure()
plt.plot(t, y, "-", label="Function")
plt.plot(t, -y, "-", label="Analytical derivative")
plt.plot(t[:tSteps1-1], dydt, "-", label="Numerical derivative")

# For more accurate graph
# t2 = np.linspace(0.05, 1.95, tSteps1-1)
# plt.plot(t2, dydt, label="Numerical derivative")

plt.title("Plot of the function \(y=e^{-t}\)\n\
          and its numerical and analytical derivatives")
plt.xlabel("Time")
plt.ylabel("\(y\)")
plt.legend()
plt.tight_layout()
# plt.show()
plt.savefig("figure1.eps")
plt.close()

# CASE 2
deltaT = 0.5
tSteps2 = int(maxTime / deltaT)

tn = np.linspace(0, maxTime, tSteps2 + 1)
yn = np.ones(tSteps2 + 1)
n = 1
for n in range(1, tSteps2 + 1):
    yn[n] = yn[n-1] - yn[n-1]*deltaT

# Figure 2
plt.figure()
plt.plot(t, y, "-", label="Analytical solution")
plt.plot(tn, yn, "-*", label="Numerical solution")
plt.title("Plot of the analytical and numerical\n\
          solution of \( \\frac{dy}{dt}=-y \)")
plt.xlabel("Time")
plt.ylabel("\(y\)")
plt.legend()
plt.tight_layout()
# plt.show()
plt.savefig("figure2.eps")
plt.close()
