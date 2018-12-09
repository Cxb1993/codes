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
    x = np.loadtxt(iFPath + "coordinateX_ny-" + str(ny) + iFExt)
    y = np.loadtxt(iFPath + "coordinateY_ny-" + str(ny) + iFExt)
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
        iFPath + velScheme + "_ny-" + str(ny) + "_numerical_" + S + "-Final"
        + iFExt)
    return s_num


def uContours(ny, velScheme, output=""):
    """Plot u velocity contours."""
    plt.figure(1, figsize=(5.39749, 5.39749/3))

    x, y = coordinateLoader(ny)
    u = numericalSolutionLoader(ny, velScheme, "u")

    plt.contourf(x[1:-1:2], y[1:-1:2], u[1:-1, 1:-1])

    plt.colorbar()
    plt.title("Time step = Final")
    plt.tight_layout()

    if output == "show":
        plt.show()
    else:
        plt.savefig(oFPath + "ny-" + str(ny) + "_numerical_uContours" + oFExt)

    plt.close(1)


def vContours(ny, velScheme, output=""):
    """Plot v velocity contours."""
    plt.figure(2, figsize=(5.39749, 5.39749/3))

    x, y = coordinateLoader(ny)
    v = numericalSolutionLoader(ny, velScheme, "v")

    plt.contourf(x[1:-1:2], y[1:-1:2], v[1:-1, 1:-1])

    plt.colorbar()
    plt.title("Time step = Final")
    plt.tight_layout()

    if output == "show":
        plt.show()
    else:
        plt.savefig(oFPath + "ny-" + str(ny) + "_numerical_vContours" + oFExt)

    plt.close(2)


def pContours(ny, velScheme, output=""):
    """Plot pressure contours."""
    plt.figure(3, figsize=(5.39749, 5.39749/3))

    x, y = coordinateLoader(ny)
    p = numericalSolutionLoader(ny, velScheme, "p")

    plt.contourf(x[1:-1:2], y[1:-1:2], p[1:-1, 1:-1])

    plt.colorbar()
    plt.title("Time step = Final")
    plt.tight_layout()

    if output == "show":
        plt.show()
    else:
        plt.savefig(oFPath + "ny-" + str(ny) + "_numerical_pContours" + oFExt)

    plt.close(3)


def velVectors(ny, velScheme, output=""):
    """Plot velocity vectors."""
    plt.figure(4, figsize=(5.39749, 5.39749))

    x, y = coordinateLoader(ny)
    u = numericalSolutionLoader(ny, velScheme, "u")
    v = numericalSolutionLoader(ny, velScheme, "v")

    # for n in range(0, 4):
    #     plt.subplot2grid((4, 1), sections[n])
    #     fileInput = files[n]
    #     u = np.loadtxt(iFPath + "U_" + fileInput + iFExt)
    #     v = np.loadtxt(iFPath + "V_" + fileInput + iFExt)
    #     u = u[1:-1, 1:-1]
    #     v = v[1:-1, 1:-1]
    #     c = np.hypot(u, v)
    #     plt.quiver(x[120::60], y[::], u[::, 120::60], v[::, 120::60],
    #                c[::, 120::60], units='height', width=0.01, angles='xy',
    #                scale=0.3)
    #     plt.title('Time step = ' + files[n] + '')

    # c = np.hypot(u[1:-1, 1:-1], v[1:-1, 1:-1])
    # plt.quiver(x[1:-1:], y[1:-1:], u[::, ::], v[::, ::],
    #            c[::, 1:-1:2], units='height', width=0.01, angles='xy',
    #            scale=0.3)
    print(x.shape[::])
    print(y.shape[::])
    print(u.shape)
    print(v.shape)

    plt.title("Time step = Final")
    plt.tight_layout()

    if output == "show":
        plt.show()
    else:
        plt.savefig(oFPath + "ny-" + str(ny) + "_numerical_velVectors" + oFExt)

    plt.close(4)


def analytical(ny):
    """Plot u contours."""
    plt.figure(6, figsize=(5.39749, 5.39749))

    x, y = coordinateLoader(ny)
    u = np.loadtxt(iFPath + "exact_U.dat")

    # x = np.linspace(0, 100, 101)
    # y = np.linspace(-10, 10, 21)
    # D = 10
    # L = 100
    # mu = 8.9e-4
    # deltaP = 1000
    # Uavg = D**2*deltaP/(12*mu*L)
    # u = [[None for i in range(len(x))] for j in range(len(y))]
    # for i in range(len(x)):
    #     for j in range(len(y)):
    #         u[j][i] = (3.0/2.0)*(1.0 - (y[j]/(D/2.0))**2)*Uavg
    plt.plot(u, y)
    # plt.contourf(x, y, u)
    plt.tight_layout
    plt.show()


def plotLabeller(title, xLabel, yLabel):
    """Add title and labels on x- and y-axes in the legend."""
    plt.title(title)
    plt.legend()
    plt.xlabel(xLabel)
    plt.ylabel(yLabel)
    plt.tight_layout()


def analyticalCalculator():
    """Calculate and return analytical solution."""
    y = np.linspace(0.0, 0.01, 101)
    D = 0.01
    x_y_ratio = 20
    L = D * x_y_ratio
    nu = 1.0e-4
    rho = 1000.0
    mu = nu * rho
    deltaP = 1226800 - 1226640
    Uavg = ((D**2) * deltaP) / (12 * mu * L)

    u = [None for j in range(len(y))]
    for j in range(int(np.ceil(len(y)/2))):
        u[j+int(len(y)/2)] = 1.5 * Uavg * (1.0 - (y[j] / (D/2.0))**2)
    i = 0
    for j in reversed(range(int(np.ceil(len(y)/2)))):
        u[i] = 1.5 * Uavg * (1.0 - (-y[j] / (D / 2.0))**2)
        i += 1

    return u, y


def analyticalCalculatedPlotter():
    """Plot calculated analytical results."""
    ua, ya = analyticalCalculator()
    plt.plot(ua, ya, label="Calculated exact")


def numericalPlotter(ny, velScheme):
    """Plot loaded numerical results."""
    u = numericalSolutionLoader(ny, velScheme, "u")
    x, y = coordinateLoader(ny)
    if velScheme == "up":
        plotLabel = "Upwind"
    elif velScheme == "ct":
        plotLabel = "Central"
    elif velScheme == "qk":
        plotLabel = "QUICK"
    plt.plot(u[1:-1, int(0.9 * u.shape[1])], y[1:-1:2], label=plotLabel)


def profileComparison(ny, output=""):
    """Compare the velocity profiles of analytical and numerical solutions."""
    plt.figure(7, figsize=(5.39749, 5.39749))
    analyticalCalculatedPlotter()
    numericalPlotter(ny, "up")
    numericalPlotter(ny, "qk")
    numericalPlotter(ny, "ct")

    # Plot loaded analytical results
    # u_a = np.loadtxt(iFPath + "up_ny-" + str(ny) + "_exact_U" + iFExt)
    # u_a_up = np.loadtxt(iFPath + "up_ny-" + str(ny) + "_exact_U" + iFExt)
    # u_a_qk = np.loadtxt(iFPath + "qk_ny-" + str(ny) + "_exact_U" + iFExt)
    # u_a_ct = np.loadtxt(iFPath + "ct_ny-" + str(ny) + "_exact_U" + iFExt)

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
        plt.savefig(oFPath + "ny-" + str(ny) + "_profilesComparison" + oFExt)

    plt.close(7)


def plotLogFiles(ny):
    """Plot log files."""
    log = np.loadtxt(iFPath + "log.dat")

    plt.figure(8)
    plt.semilogy(log[:, 2])
    plt.xlabel("Timesteps")
    plt.ylabel("Log of maximum change in pressure")
    plt.title("Log plot of maximum change in pressure versus timesteps")
    plt.tight_layout()
    plt.show()
    # plt.savefig(oFPath + "ny-" + str(ny) + "_maxP_logy" + oFExt)
    plt.close()

    plt.figure(9)
    plt.loglog(log[:, 2])
    plt.xlabel("Log of timesteps")
    plt.ylabel("Log of maximum change in pressure")
    plt.title("Log plot of maximum change in pressure versus timesteps")
    plt.tight_layout()
    plt.show()
    # plt.savefig(oFPath + "ny-" + str(ny) + "_maxP_loglog" + oFExt)
    plt.close()

    plt.figure(10)
    plt.semilogy(log[:, 3])
    plt.xlabel("Timesteps")
    plt.ylabel("Log of maximum change in velocity")
    plt.title("Log plot of maximum change in velocity versus timesteps")
    plt.tight_layout()
    plt.show()
    # plt.savefig(oFPath + "ny-" + str(ny) + "_maxP_logy" + oFExt)
    plt.close()

    plt.figure(11)
    plt.loglog(log[:, 3])
    plt.xlabel("Log of timesteps")
    plt.ylabel("Log of maximum change in velocity")
    plt.title("Log plot of maximum change in pressure versus timesteps")
    plt.tight_layout()
    plt.show()
    # plt.savefig(oFPath + "ny-" + str(ny) + "_maxP_loglog" + oFExt)
    plt.close()


# plotLogFiles(10)

ny = 10
profileComparison(ny, "show")
# uContours(ny, "up", "show")
# vContours(ny, "up", "show")
# pContours(ny, "up", "show")
# velVectors(ny, "up", "show")
