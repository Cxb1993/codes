"""Testing plotting.

Simple file to plot an array.
testing.py

Created on: 2018-03-23
Author: Syed Ahmad Raza
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams['image.cmap'] = 'jet'

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

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
plt.plot(uexact, yexact, "*")

u = np.loadtxt('../data/old/ny=100_qk_U-final.dat')
y = np.loadtxt('../data/old/ny=100_coordinateY.dat')
# print(u.shape)
plt.plot(u[1:-1, 50], y[2:-1], 'b', label='Label Me')

# plt.title(r'Plot Title (raw string for $Latex math input$)')
# plt.xlabel('x-axis label')
# plt.ylabel('y-axis label')
# plt.legend()
plt.tight_layout()
# plt.savefig("../figures/figureName.eps")
plt.show()
plt.close(1)

plt.figure(2, figsize=(5.39749, 5.39749))
plt.plot(xexact, vexact, "*")

v = np.loadtxt('../data/old/ny=100_qk_V-final.dat')
x = np.loadtxt('../data/old/ny=100_coordinateX.dat')
# print(u.shape)
plt.plot(x[2:-1], v[50, 1:-1], 'b', label='Label Me')

# plt.title(r'Plot Title (raw string for $Latex math input$)')
# plt.xlabel('x-axis label')
# plt.ylabel('y-axis label')
# plt.legend()
plt.tight_layout()
# plt.savefig("../figures/figureName.eps")
plt.show()
plt.close(2)
