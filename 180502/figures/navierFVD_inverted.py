"""Generate figures for solvers coded in C++.

Generate figures for Navier-Stokes equations solved using Finite Volume Method,
coded in C++ and data saved to .dat files.
numerical.py

Created on: 2017-11-28
Author: Syed Ahmad Raza
"""

import math
import numpy as np
# import matplotlib as mpl
import matplotlib.pyplot as plt

# plt.style.use('seaborn-paper')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


def analyticalCalculator(deltaP):
    """Calculate and return analytical solution."""
    y = np.linspace(0.0, 0.01, 101)
    D = 0.01
    x_y_ratio = 20
    L = D * x_y_ratio
    nu = 1.0e-4
    rho = 1000.0
    mu = nu * rho
    Uavg = ((D**2) * deltaP) / (12 * mu * L)

    u = [None for j in range(len(y))]
    for j in range(int(np.ceil(len(y)/2))):
        u[j+int(len(y)/2)] = 1.5 * Uavg * (1.0 - (y[j] / (D/2.0))**2)
    i = 0
    for j in reversed(range(int(np.ceil(len(y)/2)))):
        u[i] = 1.5 * Uavg * (1.0 - (-y[j] / (D / 2.0))**2)
        i += 1

    return u, y


def analyticalCalculatedPlotter(deltaP):
    """Plot calculated analytical results."""
    ua, ya = analyticalCalculator(deltaP)
    plt.plot(ua, ya, label="Calculated exact")


def numericalPlotter(nm, velScheme, color):
    """Plot loaded numerical results."""
    x, xi, y, yi = coordinateLoader(nm)
    if velScheme == "up":
        plotLabel = "Upwind"
    elif velScheme == "ct":
        plotLabel = "Central"
    elif velScheme == "qk":
        plotLabel = "QUICK"

    u = numericalSolutionLoader(nm, velScheme, "u")
    plt.plot(u[1:-1, int(0.5 * u.shape[1])], yi[1:-1], label=plotLabel,
             color=color)

    v = numericalSolutionLoader(nm, velScheme, "v")
    plt.plot(xi[1:-1], v[1:-1, int(0.5 * v.shape[1])], label=plotLabel,
             color=color)


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


def compareVelProfilesChannelFlow(nm, deltaP, output=""):
    """Compare the velocity profiles of analytical and numerical solutions."""
    plt.figure(7, figsize=(5.39749, 5.39749))
    # analyticalCalculatedPlotter(deltaP)  # 2D channel flow
    # numericalPlotter(ny, "up", "b")
    # numericalPlotter(ny, "qk", "r")
    # numericalPlotter(ny, "ct", "g")
    # gridPlotter(ny)

    # Plot loaded analytical results
    # u_a = np.loadtxt(iFPath + "up_ny=" + str(ny) + "_exact_U" + iFExt)
    # u_a_up = np.loadtxt(iFPath + "up_ny=" + str(ny) + "_exact_U" + iFExt)
    # u_a_qk = np.loadtxt(iFPath + "qk_ny=" + str(ny) + "_exact_U" + iFExt)
    # u_a_ct = np.loadtxt(iFPath + "ct_ny=" + str(ny) + "_exact_U" + iFExt)

    # plt.plot(u_a, y, label="Analytical")
    # plt.plot(u_a_up, y, label="Analytical_{up}")
    # plt.plot(u_a_qk, y, label="Analytical_{qk}")
    # plt.plot(u_a_ct, y, label="Analytical_{ct}")

    plotLabeller("Comparision of velocity profiles of analytical solution\
                 \nversus various numerical schemes for \(ny\) = " + str(nm),
                 "\(u\)-velocity (m/s)", "\(D\) (m)")

    if output == "show":
        plt.show()
    else:
        plt.savefig(oFPath + "ny=" + str(nm) + "_profilesComparison" + oFExt)

    plt.close(7)


def gridPlotter(nm):
    """Plot the grid."""
    x, xi, y, yi = coordinateLoader(nm)
    # xv, yv = np.meshgrid(x, y)
    # linesx = plt.plot(xv, yv)
    xv, yv = np.meshgrid(y[1:-1], x[1:-1])
    linesy = plt.plot(yv, xv)
    for lines in [linesy]:
        for line in lines:
            line.set_lw(0.2)
            line.set_color('k')


def coordinateLoader(nm, name):
    """Load and return coordinates from file."""
    x = np.loadtxt(iFPath + name + "_" + str(nm) + "_coordinateX"
                   + iFExt)
    y = np.loadtxt(iFPath + name + "_" + str(nm) + "_coordinateY"
                   + iFExt)

    xi = np.zeros(nm)
    for i in range(nm):
        xi[i] = (x[i] + x[i+1]) / 2.0
    yi = np.zeros(nm)
    for i in range(nm):
        yi[i] = (y[i] + y[i+1]) / 2.0

    return x, xi, y, yi


def numericalSolutionLoader(nm, name, s):
    """Load and return required numerical solution from file."""
    if s == "u":
        S = "U"
    elif s == "v":
        S = "V"
    elif s == "p":
        S = "P"
    else:
            print("ERROR: Wrong input in numericalSolutionLoader")
    numSol = np.loadtxt(iFPath + name + "_" + str(nm) + "_"
                        + velScheme + "_" + S + "-final" + iFExt)
    return numSol


def dataExtractorU(nm, name):
    """Extract u velocity data at a y cross-section."""
    u = numericalSolutionLoader(nm, name, "u")
    nt = nm + 4
    um = np.zeros(nt)

    if nt % 2 == 0:
        # Calculate the middle velocities by averaging for even number of grids
        a = int(nt/2 - 1)
        b = int(nt/2)
        for i in range(nt):
            um[i] = (u[i, a] + u[i, b]) / 2.0

    else:
        c = math.ceil(nt/2)
        for i in range(nt):
            um[i] = u[i, c]

    # Ghost cell values are also returned
    return um


def dataExtractorV(nm, name):
    """Extract v velocity data at an x cross-section."""
    v = numericalSolutionLoader(nm, name, "v")
    nt = nm + 4
    vm = np.zeros(nt)

    if nt % 2 == 0:
        # Calculate the middle velocities by averaging for even number of grids
        a = int(nt/2 - 1)
        b = int(nt/2)
        for j in range(nt):
            vm[j] = (v[a, j] + v[b, j]) / 2.0

    else:
        c = math.ceil(nt/2)
        for j in range(nt):
            vm[j] = v[c, j]

    # Ghost cell values are also returned
    return vm


def plotLabeller(xLabel, yLabel, title=""):
    """Add title and labels on x- and y-axes in the legend."""
    if title:
        plt.title(title)
    if xLabel:
        plt.xlabel(xLabel)
    if yLabel:
        plt.ylabel(yLabel)

    plt.legend()
    plt.tight_layout()


def cavityPlotterU(nm, name, lineStyle):
    """Plot the u velocity cavity flow comparison against y coordinates."""
    x, xi, y, yi = coordinateLoader(nm, name)
    uy = dataExtractorU(nm, name)

    plotLabel = name

    plt.plot(yi, uy[2:-2], label=plotLabel, linestyle=lineStyle)
    plt.gca().invert_yaxis()


def cavityPlotterV(nm, name, lineStyle):
    """Plot the v velocity cavity flow comparison against x coordinates."""
    x, xi, y, yi = coordinateLoader(nm, name)
    vx = dataExtractorV(nm, name)

    plotLabel = name

    plt.plot(vx[2:-2], xi, label=plotLabel, linestyle=lineStyle)
    plt.gca().invert_yaxis()


def uContours(nm, name, output=""):
    """Plot u velocity contours."""
    plt.figure(1, figsize=(5.39749, 5.39749))

    x, xi, y, yi = coordinateLoader(nm)

    u = numericalSolutionLoader(nm, velScheme, "u")

    plt.contourf(xi[1:-1], yi[1:-1], u[1:-1, 1:-1])

    plt.colorbar()
    plt.title("Time step = Final")
    plt.tight_layout()

    # if output == "show":
    #     plt.show()
    # else:
    #     # figureSaver(nm, velScheme, "_uContours")

    # plt.close(1)


def vContours(nm, name, output=""):
    """Plot v velocity contours."""
    plt.figure(1, figsize=(5.39749, 5.39749))

    x, xi, y, yi = coordinateLoader(nm)

    v = numericalSolutionLoader(nm, velScheme, "v")

    plt.contourf(xi[1:-1], yi[1:-1], v[1:-1, 1:-1])

    plt.colorbar()
    plt.title("Time step = Final")
    plt.tight_layout()

    if output == "show":
        plt.show()
    # else:
    #     figureSaver(nm, velScheme, "_vContours")

    # plt.close(1)


def velVectors(nm, name, output=""):
    """Plot velocity vectors."""
    plt.figure(1, figsize=(5.39749, 5.39749))

    x, xi, y, yi = coordinateLoader(nm)

    u = numericalSolutionLoader(nm, velScheme, "u")
    v = numericalSolutionLoader(nm, velScheme, "v")

    # for n in range(0, 4):
    #     plt.subplot2grid((4, 1), sections[n])
    #     fileInput = files[n]
    #     u = np.loadtxt(iFPath + "U_" + fileInput + iFExt)
    #     v = np.loadtxt(iFPath + "V_" + fileInput + iFExt)
    u = u[1:-1, 1:-1]
    v = v[1:-1, 1:-1]
    c = np.hypot(u, v)
    plt.quiver(xi[1:-1], yi[1:-1], u[::, ::], v[::, ::], c[::, ::],
               units='height', width=0.002, angles='xy', scale=10)
    #     plt.title('Time step = ' + files[n] + '')

    # c = np.hypot(u[1:-1, 1:-1], v[1:-1, 1:-1])
    # plt.quiver(x[1:-1:], y[1:-1:], u[::, ::], v[::, ::],
    #            c[::, 1:-1:2], units='height', width=0.01, angles='xy',
    #            scale=0.3)
    # print(xi.shape[::])
    # print(yi.shape[::])
    # print(u.shape)
    # print(v.shape)

    plt.title("Time step = Final")
    plt.tight_layout()

    # if output == "show":
    #     plt.show()
    # else:
    #     figureSaver(nm, velScheme, "_velVectors")

    # plt.close(1)


def compareVelProfilesCavityFlow(nxy, name):
    """Compare velocity profiles of different cases for cavity flow."""
    # outputName = ""
    # for m in range(len(nameList)):
    #     outputName = outputName + nameList[m] + "-&-"
    outputName = name

    # Plot u-velocity against y coordinates

    fig, ax1 = plt.subplots(figsize=[5.39749, 5.39749])

    # for n in range(len(nameList)):
    #     cavityPlotterU(nxy, nameList[n], lineStyleList[n])
    x, xi, y, yi = coordinateLoader(nxy, name)
    uy = dataExtractorU(nxy, name)

    ax1.plot(xexact, vexact, '.b')

    ax2 = ax1.twinx()
    ax2.invert_yaxis()
    ax2.plot(yi, uy[2:-2], linestyle='-')

    ax1.set_xlabel("\(y\)")
    ax2.set_ylabel("\(u\)-velocity")

    plt.tight_layout()
    plt.savefig(oFPath + outputName + "_" + str(nxy) + "_" + velScheme
                + "_cavityFlowU" + oFExt)

    # Plot v-velocity against x coordinates

    fig, ax1 = plt.subplots(figsize=[5.39749, 5.39749])

    # for n in range(len(nameList)):
    #     cavityPlotterV(nxy, nameList[n], lineStyleList[n])
    x, xi, y, yi = coordinateLoader(nxy, name)
    vx = dataExtractorV(nxy, name)

    ax1.plot(uexact, yexact, '.b')

    ax2 = ax1.twinx()
    ax2.invert_yaxis()
    ax2.plot(vx[2:-2], xi, linestyle='-')

    ax1.set_xlabel("\(v\)-velocity")
    ax2.set_ylabel("\(x\)")

    plt.tight_layout()
    plt.savefig(oFPath + outputName + "_" + str(nxy) + "_" + velScheme
                + "_cavityFlowV" + oFExt)


velScheme = "qk"
iFPath = "../data/"
oFPath = "../figures/"
iFExt = ".dat"
oFExt = ".eps"

compareVelProfilesCavityFlow(41, 'invert')
# uContours(nxy, "qk")
# vContours(nxy, "qk")
# velVectors(nxy, "qk")
