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

cases = ['0', '0.3', '0.6', '1']
titles = [r'$0.0$', r'$0.3$', r'$0.6$', r'$1.0$']
sections = [(0,0), (1,0), (2,0), (3,0)]

data = np.loadtxt('../data/navierFVD_data_1e-05.dat')
data = data[1:-1,1:-1]
plt.contourf(x, y, data)
plt.colorbar()
plt.title(r'$t= 1e-05$ s')

# for n in range(0,4):
#     plt.subplot2grid((4,1), sections[n])
#     fileInput = cases[n];
#     u = np.loadtxt('../data/navierFVD_U_' + fileInput + '.dat')
#     u = u[1:-1,1:-1]
#     v = np.loadtxt('../data/navierFVD_V_' + fileInput + '.dat')
#     v = v[1:-1,1:-1]
#     data = np.loadtxt('../data/navierFVD_data_' + fileInput + '.dat')
#     data = data[1:-1,1:-1]
#     plt.contourf(x, y, data)
#     plt.colorbar()
#     # plt.quiver(x, y, u, v, angles='xy', scale=20, color='b')
#     # plt.streamplot(x, y, u, v, density=(20,1), linewidth=0.2, arrowsize=0.1, minlength=0.01)
#     plt.title(r'$t= $' + titles[n] + r' s')

plt.tight_layout()
plt.savefig("../figures/navierFVD_P.eps")
