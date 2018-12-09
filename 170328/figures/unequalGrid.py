##
## filename.py
##
##  Created on: 2017-04-10
##      Author: Syed Ahmad Raza
##

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams['image.cmap'] = 'jet'

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.figure(figsize=(5.39749, 5.39749 * 2 / 3))
gridSize = (2,3)    # size of grid in the figure

x = np.loadtxt("../data/coordinatesX.dat")
y = np.loadtxt("../data/coordinatesY.dat")

xIntervalSize = 1 / (x.shape[0] - 1)
yIntervalSize = 1 / (y.shape[0] - 1)

# Main grid
plt.subplot2grid(gridSize, (0, 0), rowspan = 2, colspan = 2)
xv, yv = np.meshgrid(x, y)
linesx = plt.plot(xv, yv)
xv, yv = np.meshgrid(y, x)
linesy = plt.plot(yv, xv)
plt.title("Grid")
plt.xlabel("$x$")
plt.ylabel("$y$")
for lines in [linesx, linesy]:
    for line in lines:
        line.set_lw(0.2)
        line.set_color('k')

# x-axis graph
xValuesForX = np.arange(0.0, 1.0 + xIntervalSize, xIntervalSize)
plt.subplot2grid(gridSize, (0, 2))
plt.plot(xValuesForX, x)
plt.xticks(np.arange(2), ['$0$', r'$\frac{i_{nx}}{n_x}$'])
plt.yticks(np.arange(2), ['$0$', '$L$'])
plt.title("Function for x coordinates")
plt.xlabel(r"$f(x)=\frac{L}{2}\left[ 1 + \sin\left(\pi\left(\frac{i}{n_x} - \frac{1}{2}\right)\right)\right]$", fontsize = 10)

# y-axis graph
plt.subplot2grid(gridSize, (1, 2))
xValues = np.arange(0.0, 1.0 + yIntervalSize, yIntervalSize)
plt.plot(xValues, y)
plt.xticks(np.arange(2), ['$0$', r'$\frac{j_{ny}}{n_y}$'])
plt.yticks(np.arange(2), ['$0$', '$W$'])
plt.title("Function for y coordinates")
plt.xlabel(r"$f(y)=W\sin\left(\frac{\pi}{2}\times\frac{j}{n_y}\right)$", fontsize = 10)

plt.tight_layout()
plt.savefig("gridAndGraph.eps")
