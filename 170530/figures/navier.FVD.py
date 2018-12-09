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

files = ['1', '100', '1000', '10000']
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
    plt.quiver(x[::60], y, u[::,::60], v[::,::60], c[::,::60], units='height', width=0.01, angles='xy', scale=0.175)
    plt.title('Time step = ' + files[n] + '')

plt.tight_layout()
plt.savefig("../figures/navierFVD-velocityVectors.eps")
plt.close(4)

# plt.figure(5, figsize=(5.39749, 5.39749))
#
# for n in range(0,4):
#     plt.subplot2grid((4,1), sections[n])
#     fileInput = files[n];
#     u = np.loadtxt('../data/navierFVD_U_' + fileInput + '.dat')
#     u = u[1:-1,1:-1]
#     v = np.loadtxt('../data/navierFVD_V_' + fileInput + '.dat')
#     v = v[1:-1,1:-1]
#     c = np.hypot(u, v)
#     plt.streamplot(x, y, u, v, density=(20,1), linewidth=0.2, arrowsize=0.1, minlength=0.01)
#     # plt.title('Time step = ' + titles[n] + '')
#     plt.title('Time step = ' + files[n] + '')
#
# plt.tight_layout()
# plt.savefig("../figures/navierFVD-velocityStream.eps")
# plt.close(5)
