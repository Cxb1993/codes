import numpy as np
import matplotlib as mpl
mpl.use('PS')
mpl.rcParams['image.cmap'] = 'jet'
import matplotlib.pyplot as plt

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# Data loading

x = np.loadtxt("../161220_heatTransfer2D/data/coordinatesX.dat")
y = np.loadtxt("../161220_heatTransfer2D/data/coordinatesY.dat")

plt.figure(figsize=(5.39749, 5.39749))

plt.subplot2grid((2,2), (0,0))
fileInput = "ht2dCase02T100"
t = np.loadtxt("../161220_heatTransfer2D/data/" + fileInput + ".dat")
plt.contourf(x, y, t)
plt.colorbar()
plt.title(r"$t = 100$")

plt.subplot2grid((2,2), (0,1))
fileInput = "ht2dCase02T1000"
t = np.loadtxt("../161220_heatTransfer2D/data/" + fileInput + ".dat")
plt.contourf(x, y, t)
plt.colorbar()
plt.title(r"$t = 1000$")

plt.subplot2grid((2,2), (1,0))
fileInput = "ht2dCase02T2000"
t = np.loadtxt("../161220_heatTransfer2D/data/" + fileInput + ".dat")
plt.contourf(x, y, t)
plt.colorbar()
plt.title(r"$t = 2000$")

plt.subplot2grid((2,2), (1,1))
fileInput = "ht2dCase02T3000"
t = np.loadtxt("../161220_heatTransfer2D/data/" + fileInput + ".dat")
plt.contourf(x, y, t)
plt.colorbar()
plt.title(r"$t = 3000$")

plt.tight_layout()
plt.savefig("../161220_report/ht2dCase02.eps")