"""Generate figures for solvers coded in C++.

Generate figures for Navier-Stokes equations solved using Finite Volume Method,
coded in C++ and data saved to .dat files.
numerical.py

Created on: 2017-11-28
Author: Syed Ahmad Raza
"""

import numpy as np
# import matplotlib as mpl
import matplotlib.pyplot as plt

# mpl.rcParams['image.cmap'] = 'jet'

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

iFPath = "../data/"
oFPath = "../figures/"
iFExt = ".dat"
oFExt = ".eps"


def coordinateLoader(ny):
    """Load and return coordinates from file."""
    x = np.loadtxt(iFPath + "ny=" + str(ny) + "_coordinateX" + iFExt)
    y = np.loadtxt(iFPath + "ny=" + str(ny) + "_coordinateY" + iFExt)
    return x, y


def numericalSolutionLoader(ny, velScheme, s):
    """Load and return required numerical solution from file."""
    if s == "u":
        S = "U"
    elif s == "v":
        S = "V"
    elif s == "p":
        S = "P"
    else:
            print("ERROR: Wrong input in numericalSolutionLoader")
    s_num = np.loadtxt(
        iFPath + "ny=" + str(ny) + "_" + velScheme + "_" + S + "-final"
        + iFExt)
    return s_num


def plotLabeller(title, xLabel, yLabel):
    """Add title and labels on x- and y-axes in the legend."""
    plt.title(title)
    plt.legend()
    plt.xlabel(xLabel)
    plt.ylabel(yLabel)
    plt.tight_layout()


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


def numericalPlotter(ny, velScheme, color):
    """Plot loaded numerical results."""
    nx = ny
    x, y = coordinateLoader(ny)
    if velScheme == "up":
        plotLabel = "Upwind"
    elif velScheme == "ct":
        plotLabel = "Central"
    elif velScheme == "qk":
        plotLabel = "QUICK"

    u = numericalSolutionLoader(ny, velScheme, "u")
    yi = np.zeros(ny)
    for i in range(ny):
        yi[i] = (y[i] + y[i+1]) / 2.0
    plt.plot(u[1:-1, int(0.5 * u.shape[1])], yi[1:-1], label=plotLabel,
             color=color)

    v = numericalSolutionLoader(nx, velScheme, "v")
    xi = np.zeros(nx)
    for i in range(nx):
        xi[i] = (x[i] + x[i+1]) / 2.0
    plt.plot(xi[1:-1], v[1:-1, int(0.5 * v.shape[1])], label=plotLabel,
             color=color)


def gridPlotter(ny):
    """Plot the grid."""
    x, y = coordinateLoader(ny)
    # xv, yv = np.meshgrid(x, y)
    # linesx = plt.plot(xv, yv)
    xv, yv = np.meshgrid(y[1:-1], x[1:-1])
    linesy = plt.plot(yv, xv)
    for lines in [linesy]:
        for line in lines:
            line.set_lw(0.2)
            line.set_color('k')


def profileComparison(ny, deltaP, output=""):
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
                 \nversus various numerical schemes for \(ny\) = " + str(ny),
                 "\(u\)-velocity (m/s)", "\(D\) (m)")

    if output == "show":
        plt.show()
    else:
        plt.savefig(oFPath + "ny=" + str(ny) + "_profilesComparison" + oFExt)

    plt.close(7)


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

um = [0.005261021234852, -0.005261021234852, -0.01459727566199,
      -0.023097051487705, -0.03103322946768, -0.038634671835561,
      -0.046094397454977, -0.053575069840726, -0.061211941039321,
      -0.069113519414523, -0.077360958436328, -0.086003727788147,
      -0.095051655675385, -0.104463922650381, -0.114140912563216,
      -0.123939498504926, -0.133640843395742, -0.142943584130178,
      -0.151467976558968, -0.158768838715268, -0.164359990948764,
      -0.167749980009032, -0.168476793069092, -0.166117339324564,
      -0.16031455508894, -0.15079501098321, -0.137356883474078,
      -0.119820277170244, -0.097927373281853, -0.07117987787095,
      -0.038603286233231, 0.001575832553336, 0.052353028845071,
      0.117916206681697, 0.204169607717874, 0.318431186368153,
      0.467948940931559, 0.656653437712442, 0.880082050041757,
      1.11991794995824]

vm = [-0.019810895429051, 0.019810895429051, 0.052964229277109,
      0.079962853756459, 0.101265067337296, 0.117470029766532,
      0.129240228379919, 0.137239827583481, 0.142088465710401,
      0.144330416719219, 0.144417774754493, 0.142705134694506,
      0.139452671391822, 0.13483460032218, 0.128950546899336,
      0.121838080877072, 0.11348541802632, 0.103843941907158,
      0.092840722329121, 0.080391614116642, 0.066415823039785,
      0.050853013137642, 0.033684002643871, 0.014955372681887,
      -0.005194838151024, -0.026524281163, -0.048124780605274,
      -0.069323280917563, -0.089133171268464, -0.106376144225763,
      -0.119782621934902, -0.128121725740745, -0.130359432105232,
      -0.125821636284078, -0.114329186344116, -0.096275862999233,
      -0.072635443844495, -0.04490147892252, -0.014974557489423,
      0.014974557489423]

ny = 40
nx = ny

x, y = coordinateLoader(ny)
xi = np.zeros(nx)
for i in range(nx):
    xi[i] = (x[i] + x[i+1]) / 2.0
yi = np.zeros(ny)
for i in range(ny):
    yi[i] = (y[i] + y[i+1]) / 2.0


def cavityPlotterY(ny, velScheme, color, output=""):
    """Plot the Y cavity flow comparison."""
    plt.figure(1, figsize=(5.39749, 5.39749))

    plt.plot(uexact, yexact, "*")

    if velScheme == "up":
        plotLabel = "Upwind"
    elif velScheme == "ct":
        plotLabel = "Central"
    elif velScheme == "qk":
        plotLabel = "QUICK"

    # u = numericalSolutionLoader(ny, velScheme, "u")

    # plt.plot(u[1:-1, int(0.5 * u.shape[1])], yi[1:-1], label=plotLabel,
    #          color=color)
    plt.plot(um[1:-1], yi[1:-1], label=plotLabel, color=color)

    plotLabeller("Comparision of velocity profiles of Ghia et. al solution\
                 \nversus numerical schemes", "\(u\)-velocity",
                 "\(Y\)")

    if output == "show":
        plt.show()
    else:
        plt.savefig(oFPath + "n,xy=" + str(ny) + "_" + velScheme
                    + "_cavityFlowY" + oFExt)
    plt.close(1)


def cavityPlotterX(nx, velScheme, color, output=""):
    """Plot the X cavity flow comparison."""
    plt.figure(2, figsize=(5.39749, 5.39749))

    plt.plot(xexact, vexact, ".")

    if velScheme == "up":
        plotLabel = "Upwind"
    elif velScheme == "ct":
        plotLabel = "Central"
    elif velScheme == "qk":
        plotLabel = "QUICK"

    # v = numericalSolutionLoader(nx, velScheme, "v")

    plt.plot(xi[1:-1], vm[1:-1], label=plotLabel, color=color)

    plotLabeller("Comparision of velocity profiles of Ghia et. al solution\
                 \nversus numerical schemes", "\(X\)",
                 "\(v\)-velocity")

    if output == "show":
        plt.show()
    else:
        plt.savefig(oFPath + "n,xy=" + str(nx) + "_" + velScheme
                    + "_cavityFlowX" + oFExt)
    plt.close(2)


def uContours(ny, velScheme, output=""):
    """Plot u velocity contours."""
    plt.figure(1, figsize=(5.39749, 5.39749))

    u = numericalSolutionLoader(ny, velScheme, "u")

    plt.contourf(xi[1:-1], yi[1:-1], u[1:-1, 1:-1])

    plt.colorbar()
    plt.title("Time step = Final")
    plt.tight_layout()

    if output == "show":
        plt.show()
    else:
        plt.savefig(oFPath + "n,xy=" + str(ny) + "_" + velScheme
                    + "_uContours" + oFExt)

    plt.close(1)


def vContours(ny, velScheme, output=""):
    """Plot v velocity contours."""
    plt.figure(2, figsize=(5.39749, 5.39749))

    v = numericalSolutionLoader(ny, velScheme, "v")

    plt.contourf(xi[1:-1], yi[1:-1], v[1:-1, 1:-1])

    plt.colorbar()
    plt.title("Time step = Final")
    plt.tight_layout()

    if output == "show":
        plt.show()
    else:
        plt.savefig(oFPath + "n,xy=" + str(ny) + "_" + velScheme
                    + "_vContours" + oFExt)

    plt.close(2)


def velVectors(ny, velScheme, output=""):
    """Plot velocity vectors."""
    plt.figure(4, figsize=(5.39749, 5.39749))

    u = numericalSolutionLoader(ny, velScheme, "u")
    v = numericalSolutionLoader(ny, velScheme, "v")

    # for n in range(0, 4):
    #     plt.subplot2grid((4, 1), sections[n])
    #     fileInput = files[n]
    #     u = np.loadtxt(iFPath + "U_" + fileInput + iFExt)
    #     v = np.loadtxt(iFPath + "V_" + fileInput + iFExt)
    u = u[1:-1, 1:-1]
    v = v[1:-1, 1:-1]
    c = np.hypot(u, v)
    plt.quiver(xi[1:-1], yi[1:-1], u[::, ::], v[::, ::], c[::, ::],
               units='height', width=0.003, angles='xy', scale=3)
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

    if output == "show":
        plt.show()
    else:
        plt.savefig(oFPath + "n,xy=" + str(ny) + "_" + velScheme
                    + "_velVectors" + oFExt)

    plt.close(4)


# cavityPlotterY(ny, "up", "r")
# cavityPlotterX(nx, "up", "r")
# vContours(ny, "up")
# uContours(ny, "up")
velVectors(ny, "qk")
