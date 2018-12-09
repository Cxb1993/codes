import numpy as np
import matplotlib as mpl
#mpl.use('PS')
mpl.rcParams['image.cmap'] = 'jet'
import matplotlib.pyplot as plt

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# Loading the coordinates
coordinatesX = ("../161227_heatTransfer2D/data/coordinatesXCase01.dat",
              "../161227_heatTransfer2D/data/coordinatesXCase02.dat")
coordinatesY = ("../161227_heatTransfer2D/data/coordinatesYCase01.dat",
              "../161227_heatTransfer2D/data/coordinatesYCase02.dat")

# Plotting the contours of numerical solution
plt.figure(figsize=(5.39749, 5.39749))
  
cases = [("1T0100", "1T1000", "1T2000", "1T3741.219"),
         ("2T0100", "2T1000", "2T2000", "2T3883.242")] 
titles = [(r"$100$", r"$1000$", r"$2000$", r"$3741.219$"),
          (r"$100$", r"$1000$", r"$2000$", r"$3883.242$")]
sections = [(0,0), (0,1), (1,0), (1,1)]
  
for i in range(0,2):    
    for j in range(0,4):
        plt.subplot2grid((2,2), sections[j])
        x = np.loadtxt(coordinatesX[i])
        y = np.loadtxt(coordinatesY[i])
        fileInput = cases[i][j];
        t = np.loadtxt("../161227_heatTransfer2D/data/ht2dCase0"
                       + fileInput + ".dat")
        plt.contourf(x, y, t, np.linspace(25.0, 50.0, 26))
        plt.colorbar(ticks = [25, 30, 35, 40, 45, 50])
        plt.title(r"$t = $" + titles[i][j] + r" s")
  
    plt.tight_layout()
    plt.savefig("../161227_report/ht2dCase0" + str(i+1) + ".eps")

# Plotting the contours of analytical solution
plt.figure(figsize=(5.39749, 5.39749))
x = np.loadtxt(coordinatesX[0])
y = np.loadtxt(coordinatesY[0])
t = np.loadtxt("../161227_heatTransfer2D/data/ht2dCase01Analytical.dat")
v = np.linspace(25.0, 50.0, 26)
plt.contourf(x, y, t, v)
plt.colorbar(ticks=v)
plt.title(r"Analytical solution")
plt.savefig("../161227_report/ht2dCase01Analytical.eps")