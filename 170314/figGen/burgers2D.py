import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams['image.cmap'] = 'jet'

plt.rc('text', usetex=True)
plt.rc('font', family='sertif')

x = np.loadtxt("../170314_burgers2D/data/coordinatesX.dat")
y = np.loadtxt("../170314_burgers2D/data/coordinatesY.dat")

u = np.loadtxt("../170314_burgers2D/data/analyticalU.dat")
un = np.loadtxt("../170314_burgers2D/data/numericalU_0.1.dat")
v = np.loadtxt("../170314_burgers2D/data/analyticalV.dat", unpack=True)
vn = np.loadtxt("../170314_burgers2D/data/numericalV_0.1.dat", unpack=True)

plt.figure(figsize=(5.39749, 5.39749))
plt.plot(x, u[:][25], label='Analytical result')
plt.plot(x, un[:][25], '.r', label='Numerical result')
plt.title("Profile of $u$ along $x$-axis at $y= $" + str(y[25]))
plt.xlabel("x-axis")
plt.ylabel("$u$")
plt.legend()
plt.tight_layout()
plt.savefig("../170314_report/comparisonU.eps")

plt.figure(figsize=(5.39749, 5.39749))
plt.plot(y, v[25][:], label='Analytical result')
plt.plot(y, vn[:][25], '.r', label='Numerical result')
plt.title("Profile of $v$ along $y$-axis at $x= $" + str(x[25]))
plt.xlabel("y-axis")
plt.ylabel("$v$")
plt.legend()
plt.tight_layout()
plt.savefig("../170314_report/comparisonV.eps")

plt.figure(figsize=(5.39749, 5.39749))
q = plt.streamplot(x, y, u, v, color='b')
plt.tight_layout()
plt.savefig("../170314_report/streamplot.eps")