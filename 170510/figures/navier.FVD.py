##
## navierFVD.py
##
##  Created on: 2017-04-29
##      Author: Syed Ahmad Raza
##

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams['image.cmap'] = 'jet'

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.figure(figsize=(5.39749, 5.39749))

x = np.loadtxt('../data/coordinateX.dat')
x = x[1:-1]
y = np.loadtxt('../data/coordinateY.dat')
y = y[1:-1]

files = ['0', '1', '2', '100']
titles = [r'$0$', r'$1$', r'$2$', r'$100$']
sections = [(0,0), (1,0), (2,0), (3,0)]

for n in range(0,4):
    plt.subplot2grid((4,1), sections[n])
    fileInput = files[n];
    u = np.loadtxt('../data/navierFVD_U_' + fileInput + '.dat')
    u = u[1:-1,1:-1]
    # v = np.loadtxt('../data/navierFVD_V_' + fileInput + '.dat')
    # v = v[1:-1,1:-1]
    # p = np.loadtxt('../data/navierFVD_P_' + fileInput + '.dat')
    # p = p[1:-1,1:-1]
    # plt.quiver(x, y, u, v, angles='xy', scale=100, color='b')
    plt.contourf(x, y, u)
    plt.colorbar()
    # plt.streamplot(x, y, u, v, density=(20,1), linewidth=0.2, arrowsize=0.1, minlength=0.01)
    plt.title(r'Time step = ' + titles[n] + r'')

plt.tight_layout()
plt.savefig("../figures/navierFVD-Uvelocity.eps")
