"""Rewriting a simple plotter for 2D cavity problem.

A simple plotter for 2D cavity problem with no custom functions
filename.py

Created on: 2018-05-11
Author: Syed Ahmad Raza
"""

import math
import numpy as np
# import matplotlib as mpl
import matplotlib.pyplot as plt

# mpl.rcParams['image.cmap'] = 'jet'

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# The code only works for ODD number of grids
nxy = 41
fileUniqueName = "invertedBottom"

# Cavity flow reference data input

xexact = [0.0, 0.0625, 0.0703, 0.0781, 0.0938, 0.1563, 0.2266, 0.2344, 0.5,
          0.8047, 0.8594, 0.9063, 0.9453, 0.9531, 0.9609, 0.9688, 1.0]

yexact = [0.0, 0.0547, 0.0625, 0.0703, 0.1016, 0.1719, 0.2813, 0.4531,
          0.5, 0.6172, 0.7344, 0.8516, 0.9531, 0.9609, 0.9688, 0.9766, 1.0]

uexact = [0.0, -0.03717, -0.04192, -0.04775, -0.06434, -0.10150, -0.15662,
          -0.21090, -0.20581, -0.13641, 0.00332, 0.23151, 0.68717, 0.73722,
          0.78871, 0.84123, 1.0]

vexact = [0.0, 0.09233, 0.10091, 0.10890, 0.12317, 0.16077, 0.17507, 0.17527,
          0.05454, -0.24533, -0.22445, -0.16914, -0.10313, -0.08864, -0.07391,
          -0.05906, 0.0]

x = np.loadtxt("../data/" + fileUniqueName + "_" + str(nxy)
               + "_coordinateX.dat")
y = np.loadtxt("../data/" + fileUniqueName + "_" + str(nxy)
               + "_coordinateY.dat")
xi = np.zeros(nxy)
for i in range(nxy):
    xi[i] = (x[i] + x[i+1]) / 2.0
yi = np.zeros(nxy)
for i in range(nxy):
    yi[i] = (y[i] + y[i+1]) / 2.0

# print(x)
# print(xi)
# print(y)
# print(yi)

u = np.loadtxt("../data/" + fileUniqueName + "_" + str(nxy)
               + "_qk_U-final.dat")
u = u[2:-2, 1:-2]
# print(u)
um = np.zeros(nxy)
a = math.floor(nxy/2)
b = math.ceil(nxy/2)
for i in range(nxy):
    um[i] = (u[i, a] + u[i, b]) / 2.0
umn = np.zeros(nxy)
m = nxy
for i in range(nxy):
    umn[i] = um[m-1]
    m = m - 1

v = np.loadtxt("../data/" + fileUniqueName + "_" + str(nxy)
               + "_qk_V-final.dat")
v = v[1:-2, 2:-2]
# print(v)
vm = np.zeros(nxy)
a = math.floor(nxy/2)
b = math.ceil(nxy/2)
for j in range(nxy):
    vm[j] = (v[a, j] + v[b, j]) / 2.0
vmn = np.zeros(nxy)
m = nxy
for i in range(nxy):
    vmn[i] = vm[m-1]
    m = m - 1

fig, ax1 = plt.subplots(figsize=[5.39749, 5.39749])

ax1.plot(xexact, vexact, '.r', label="Ghia et. al")

ax2 = ax1.twinx()
ax2.invert_yaxis()
ax2.plot(yi, vmn, '-b', label=fileUniqueName)

plt.legend()
plt.tight_layout()
# plt.xlabel("\(u\)-velocity")
# plt.ylabel("\(y\)")
# plt.legend()
plt.tight_layout()
# plt.show()
plt.savefig("../figures/" + fileUniqueName + "_" + str(nxy)
            + "_cavityFlowV.eps")
plt.close()

fig, ax1 = plt.subplots(figsize=[5.39749, 5.39749])

ax1.plot(uexact, yexact, '.r', label="Ghia et. al")

ax2 = ax1.twiny()
ax2.invert_xaxis()
ax2.plot(umn, yi, '-b', label=fileUniqueName)
# plt.xlabel("\(x\)")
# plt.ylabel("\(v\)-velocity")
# plt.legend()
plt.tight_layout()
# plt.show()
plt.savefig("../figures/" + fileUniqueName + "_" + str(nxy)
            + "_cavityFlowU.eps")
plt.close()
