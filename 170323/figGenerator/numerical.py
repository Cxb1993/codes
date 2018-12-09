import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams['image.cmap'] = 'jet'

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.figure(figsize=(5.39749, 5.39749))

x = np.linspace(0, 1, 50, True)
y = np.linspace(0, 1, 50, True)
v = np.linspace(1, 2, 21)
t = np.loadtxt("../data/laplaceSORFinal.dat")

plt.contourf(x, y, t, v)
plt.colorbar(ticks=v)
plt.title("Numerical solution for Laplace Equation using SOR method")
plt.tight_layout()
plt.savefig("../170323_report/sorContours.eps")
