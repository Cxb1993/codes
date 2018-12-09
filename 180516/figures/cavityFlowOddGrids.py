"""Rewriting a simple plotter for 2D cavity problem.

A simple plotter for 2D cavity problem with no custom functions.
The code only works for ODD number of grids for now.

cavityFlowOddGrids.py

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


def gridPlotter(nm):
    """Plot the grid."""
    # x, xi, y, yi = coordinateLoader(nm)
    # xv, yv = np.meshgrid(x, y)
    # linesx = plt.plot(xv, yv)
    # xv, yv = np.meshgrid(y[1:-1], x[1:-1])
    # linesy = plt.plot(yv, xv)
    # for lines in [linesy]:
    #     for line in lines:
    #         line.set_lw(0.2)
    #         line.set_color('k')


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


plt.figure(1, figsize=(5.39749, 5.39749))
plt.plot(uexact, yexact, '.r', label="Ghia et. al", markersize=7)


# nxy = 41
fileUniqueName = "unequal"


def uRepeater(nxy, fileUniqueName):
    """Function to allow plotting of several u profiles."""
    y = np.loadtxt("../data/" + fileUniqueName + "_" + str(nxy)
                   + "_coordinateY.dat")
    # yi = np.loadtxt("../data/" + fileUniqueName + "_" + str(nxy)
    #                 + "_coordinateYs.dat")
    yi = np.zeros(nxy)
    for i in range(nxy):
        yi[i] = (y[i] + y[i + 1]) / 2.0

    u = np.loadtxt("../data/" + fileUniqueName + "_" + str(nxy)
                   + "_qk_U-final.dat")
    u = u[2:-2, 1:-2]

    um = np.zeros(nxy)
    a = math.floor(nxy / 2)
    b = math.ceil(nxy / 2)
    for i in range(nxy):
        um[i] = (u[i, a] + u[i, b]) / 2.0

    plt.plot(um, yi, '-', label=str(nxy) + "\(\\times\)" + str(nxy)
             + " grids", linewidth=0.5)

    plt.xlabel("\(u\)-velocity")
    plt.ylabel("\(y\)")
    plt.legend()

    xv, yv = np.meshgrid(yi, um)
    linesy = plt.plot(yv, xv)
    for lines in [linesy]:
        for line in lines:
            line.set_lw(0.2)
            line.set_color('0.5')


# uRepeater(11, fileUniqueName)
# uRepeater(21, fileUniqueName)
uRepeater(41, fileUniqueName)
# uRepeater(61, fileUniqueName)
# uRepeater(81, fileUniqueName)
# uRepeater(101, fileUniqueName)
# uRepeater(121, fileUniqueName)
# uRepeater(141, fileUniqueName)

plt.legend()
plt.tight_layout()
# plt.show()

plt.tight_layout()
plt.savefig("../figures/" + fileUniqueName + "_cavityFlowU.eps")
plt.close(1)


plt.figure(2, figsize=(5.39749, 5.39749))
plt.plot(xexact, vexact, '.r', label="Ghia et. al", markersize=7)


def vRepeater(nxy, fileUniqueName):
    """Function to allow plotting of several v profiles."""
    x = np.loadtxt("../data/" + fileUniqueName + "_" + str(nxy)
                   + "_coordinateX.dat")
    # xi = np.loadtxt("../data/" + fileUniqueName + "_" + str(nxy)
    #                 + "_coordinateXs.dat")
    xi = np.zeros(nxy)
    for i in range(nxy):
        xi[i] = (x[i] + x[i + 1]) / 2.0

    v = np.loadtxt("../data/" + fileUniqueName + "_" + str(nxy)
                   + "_qk_V-final.dat")
    v = v[1:-2, 2:-2]

    vm = np.zeros(nxy)
    a = math.floor(nxy / 2)
    b = math.ceil(nxy / 2)
    for j in range(nxy):
        vm[j] = (v[a, j] + v[b, j]) / 2.0

    plt.plot(xi, vm, '-', label=str(nxy) + "\(\\times\)" + str(nxy)
             + " grids", linewidth=0.5)

    plt.xlabel("\(x\)")
    plt.ylabel("\(v\)-velocity")
    plt.legend()

    xv, yv = np.meshgrid(vm, xi)
    linesy = plt.plot(yv, xv)
    for lines in [linesy]:
        for line in lines:
            line.set_lw(0.2)
            line.set_color('0.5')


# vRepeater(11, fileUniqueName)
# vRepeater(21, fileUniqueName)
vRepeater(41, fileUniqueName)
# vRepeater(61, fileUniqueName)
# vRepeater(81, fileUniqueName)
# vRepeater(101, fileUniqueName)
# vRepeater(121, fileUniqueName)
# vRepeater(141, fileUniqueName)
# plt.show()
plt.tight_layout()
plt.savefig("../figures/" + fileUniqueName + "_cavityFlowV.eps")
plt.close(2)
