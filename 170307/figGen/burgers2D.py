import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import math

mpl.rcParams['image.cmap'] = 'jet'
 
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

lm = 25.0

x = np.linspace(-1.0, 1.0, 51)
y = np.linspace(0.0, 2.0, 51)

#q = plt.quiver(X,Y,u,v,angles='xy',scale=1000,color='r')

# Plotting analytical solution
u = np.loadtxt("../170307_burgers2D/data/analyticalU.dat")
v = np.loadtxt("../170307_burgers2D/data/analyticalV.dat")

plt.figure(figsize=(5.39749, 5.39749))

q = plt.quiver(x, y, u, v, angles='xy', scale=7, color='b')

plt.tight_layout()
plt.savefig("../170307_report/analytical.eps")

# # Plotting numerical solution 0.2
# u = np.loadtxt("../170307_burgers2D/data/numericalU_0.2.dat")
# v = np.loadtxt("../170307_burgers2D/data/numericalV_0.2.dat")
# 
# plt.figure(figsize=(5.39749, 5.39749))
# 
# q = plt.quiver(x, y, u, v, angles='xy', scale=7, color='b')
# 
# plt.tight_layout()
# plt.savefig("../170307_report/numerical_0.2.eps")

# Plotting the numerical solution
plt.figure(figsize=(5.39749, 5.39749))
  
cases = ["1000", "2000", "3000", "5000"] 
titles = [r"$1000$", r"$2000$", r"$3000$", r"$5000$"]
sections = [(0,0), (0,1), (1,0), (1,1)]
  
    
for j in range(0, 4):
    plt.subplot2grid((2, 2), sections[j])
    fileInput = cases[j];
    u = np.loadtxt("../170307_burgers2D/data/numericalU_"
                   + fileInput + ".dat")
    v = np.loadtxt("../170307_burgers2D/data/numericalV_"
                   + fileInput + ".dat")
    q = plt.quiver(x, y, u, v, angles='xy', scale=7, color='b')
    plt.title(r"$t = $" + titles[j] + r" s")
  
    plt.tight_layout()
    plt.savefig("../170307_report/numerical.eps")
