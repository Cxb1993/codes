##
## navierFVD.py
##
##  Created on: 2017-04-29
##      Author: Syed Ahmad Raza
##

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# mpl.rcParams['image.cmap'] = 'jet'

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

x = np.loadtxt('../data/coordinateX.dat')
x = x[1:-1]
y = np.loadtxt('../data/coordinateY.dat')
y = y[1:-1]

files = ['1', '500', '3500', '7292']
# titles = ['$0$', '$1$', '$100$', '$200$']
sections = [(0,0), (1,0), (2,0), (3,0)]

plt.figure(1, figsize=(5.39749, 5.39749))

for n in range(0,4):
    plt.subplot2grid((4,1), sections[n])
    fileInput = files[n];
    u = np.loadtxt('../data/navierFVD_U_' + fileInput + '.dat')
    u = u[1:-1,1:-1]
    plt.contourf(x, y, u)
    plt.colorbar()
    # plt.title('Time step = ' + titles[n] + '')
    plt.title('Time step = ' + files[n] + '')

plt.tight_layout()
plt.savefig("../figures/navierFVD-UvelocityContours.eps")
plt.close(1)

plt.figure(2, figsize=(5.39749, 5.39749))

for n in range(0,4):
    plt.subplot2grid((4,1), sections[n])
    fileInput = files[n];
    v = np.loadtxt('../data/navierFVD_V_' + fileInput + '.dat')
    v = v[1:-1,1:-1]
    plt.contourf(x, y, v)
    plt.colorbar()
    # plt.title('Time step = ' + titles[n] + '')
    plt.title('Time step = ' + files[n] + '')

plt.tight_layout()
plt.savefig("../figures/navierFVD-VvelocityContours.eps")
plt.close(2)

plt.figure(3, figsize=(5.39749, 5.39749))

for n in range(0,4):
    plt.subplot2grid((4,1), sections[n])
    fileInput = files[n];
    p = np.loadtxt('../data/navierFVD_P_' + fileInput + '.dat')
    p = p[1:-1,1:-1]
    plt.contourf(x, y, p)
    plt.colorbar()
    # plt.title('Time step = ' + titles[n] + '')
    plt.title('Time step = ' + files[n] + '')

plt.tight_layout()
plt.savefig("../figures/navierFVD-pressureContours.eps")
plt.close(3)

plt.figure(4, figsize=(5.39749, 5.39749))

for n in range(0,4):
    plt.subplot2grid((4,1), sections[n])
    fileInput = files[n];
    u = np.loadtxt('../data/navierFVD_U_' + fileInput + '.dat')
    u = u[1:-1,1:-1]
    v = np.loadtxt('../data/navierFVD_V_' + fileInput + '.dat')
    v = v[1:-1,1:-1]
    c = np.hypot(u, v)
    plt.quiver(x[::60], y, u[::,::60], v[::,::60], c[::,::60], units='height', width=0.01, angles='xy', scale=3.0)
    plt.title('Time step = ' + files[n] + '')

plt.tight_layout()
plt.savefig("../figures/navierFVD-velocityVectors.eps")
plt.close(4)

plt.figure(5, figsize=(5.39749, 5.39749))

for n in range(0,4):
    plt.subplot2grid((4,1), sections[n])
    fileInput = files[n];
    u = np.loadtxt('../data/navierFVD_U_' + fileInput + '.dat')
    u = u[1:-1,1:-1]
    v = np.loadtxt('../data/navierFVD_V_' + fileInput + '.dat')
    v = v[1:-1,1:-1]
    c = np.hypot(u, v)
    plt.streamplot(x, y, u, v, density=(20,1), linewidth=0.2, arrowsize=0.1, minlength=0.01)
    # plt.title('Time step = ' + titles[n] + '')
    plt.title('Time step = ' + files[n] + '')

plt.tight_layout()
plt.savefig("../figures/navierFVD-velocityStream.eps")
plt.close(5)

# # Analytical solution
#
# plt.figure(6, figsize=(5.39749, 5.39749))
#
# u = np.loadtxt('../data/navierAnalytical_U.dat')
# u = u[1:-1,1:-1]
#
#
# x = np.linspace(0,100,101)
# y = np.linspace(-10,10,21)
# D = 10
# L = 100
# mu = 8.9e-4
# deltaP = 1000
# Uavg = D**2*deltaP/(12*mu*L)
# u = [[None for i in range(len(x))] for j in range(len(y))]
# for i in range(len(x)):
#     for j in range(len(y)):
#         u[j][i] = (3.0/2.0)*(1.0 - (y[j]/(D/2.0))**2)*Uavg
# # plt.plot(u, y)
# plt.contourf(x, y, u)
# plt.tight_layout
# plt.show()
