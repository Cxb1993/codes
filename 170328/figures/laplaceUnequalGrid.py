##
## laplaceUnequalGrid.py
##
##  Created on: 2017-04-09
##      Author: Syed Ahmad Raza
##

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams['image.cmap'] = 'jet'

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.figure(figsize=(5.39749, 5.39749))

y = np.loadtxt('../data/coordinatesY.dat')

analytical = np.loadtxt('../data/laplaceAnalyticalError.dat')
plt.plot(y, analytical, '.r', label='Analytical')

numerical = np.loadtxt('../data/laplaceNumericalError.dat')
plt.plot(y, numerical, label=r'Numerical: tolerance = $1 \times 10^{-9}$')

plt.title('Profile of $\phi$ along $y$-axis at $x = 0.507933$')
plt.xlabel('y-axis values')
plt.ylabel('$\phi$')
plt.legend()
plt.tight_layout()
plt.savefig('../figures/laplaceUnequalConvergence.eps')

# x = np.loadtxt('../data/coordinatesX.dat')
# y = np.loadtxt('../data/coordinatesY.dat')
# v = np.linspace(1, 2, 21)
# t = np.loadtxt('../data/laplaceUnequalGridFinal.dat')
#
# plt.contourf(x, y, t, v)
# plt.colorbar(ticks=v)
# plt.title('Numerical solution for Laplace Equation')
# plt.tight_layout()
# plt.savefig('laplaceUnequalContours.eps')
