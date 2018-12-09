import math
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams['image.cmap'] = 'jet'

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

U = 1

# liq
Pr = [5.49, 0.0248, 155]
v = [8.02e-7, 1.12e-7, 1.1e-5]
liquid = [r'Water', r'Hg', r'Oil']

# var
a1 = [3.464, 4.646, 4.800]
a2 = [0.577, 0.646, 0.654]

titles = [r'$\delta$ versus $x$ with $m=$',
          r'$\delta_T$ versus $x$ with $m=$',
          r'$C_f$ versus $x$ with $m=$']
m = [r'$n$', r'$(\frac{n}{2})(3-n^2)$', r'sin$(\frac{n\pi}{2})$']

ylabels = [r'$\delta$',
           r'$\delta_T$',
           r'$C_f$']
# fileNames = ['delta',
#              'Cfx',
#              'deltaT']

x = np.linspace(0.1, 1, 20)
deltas = []
Cfxs = []
deltaTs = []

gridSize = (5,4)
plt.figure(figsize=(7.26913, 10.695))

for liq in range(0,3):
    delta = []
    deltaT = []
    Cfx = []
    for j in range(20):
        delta.append(x[j]*a1[0]/math.sqrt((U*x[j])/v[liq]))
        deltaT.append(delta[j]/pow(Pr[liq],(1/3)))
        Cfx.append(a2[0]/math.sqrt((U*x[j])/v[liq]))
    deltas.append(delta)
    deltaTs.append(deltaT)
    Cfxs.append(Cfx)
        
plt.subplot2grid(gridSize, (0,0), rowspan = 1, colspan = 2)
for liq in range(0,3):
    plt.plot(x, deltas[liq], label=liquid[liq])
    plt.title(titles[0] + m[0])
    plt.xlabel(r'$x$')
    plt.ylabel(ylabels[0])
    plt.legend()
    
plt.subplot2grid(gridSize, (0,2), rowspan = 1, colspan = 2)
for liq in range(0,3):
    plt.plot(x, deltaTs[liq], label=liquid[liq])
    plt.title(titles[1] + m[0])
    plt.xlabel(r'$x$')
    plt.ylabel(ylabels[1])
    plt.legend()
    
plt.subplot2grid(gridSize, (3,0), rowspan = 1, colspan = 2)
for liq in range(0,3):
    plt.plot(x, Cfxs[liq], label=liquid[liq])
    plt.title(titles[2] + m[0])
    plt.xlabel(r'$x$')
    plt.ylabel(ylabels[2])
    plt.legend()

for liq in range(0,3):
    delta = []
    deltaT = []
    Cfx = []
    for j in range(20):
        delta.append(x[j]*a1[1]/math.sqrt((U*x[j])/v[liq]))
        deltaT.append(delta[j]/pow(Pr[liq],(1/3)))
        Cfx.append(a2[1]/math.sqrt((U*x[j])/v[liq]))
    deltas.append(delta)
    deltaTs.append(deltaT)
    Cfxs.append(Cfx)

plt.subplot2grid(gridSize, (1,0), rowspan = 1, colspan = 2)
for liq in range(0,3):
    plt.plot(x, deltas[liq], label=liquid[liq])
    plt.title(titles[0] + m[1])
    plt.xlabel(r'$x$')
    plt.ylabel(ylabels[0])
    plt.legend()
    
plt.subplot2grid(gridSize, (1,2), rowspan = 1, colspan = 2)
for liq in range(0,3):
    plt.plot(x, deltaTs[liq], label=liquid[liq])
    plt.title(titles[1] + m[1])
    plt.xlabel(r'$x$')
    plt.ylabel(ylabels[1])
    plt.legend()
    
plt.subplot2grid(gridSize, (3,2), rowspan = 1, colspan = 2)
for liq in range(0,3):
    plt.plot(x, Cfxs[liq], label=liquid[liq])
    plt.title(titles[2] + m[1])
    plt.xlabel(r'$x$')
    plt.ylabel(ylabels[2])
    plt.legend()

for liq in range(0,3):
    delta = []
    deltaT = []
    Cfx = []
    for j in range(20):
        delta.append(x[j]*a1[2]/math.sqrt((U*x[j])/v[liq]))
        deltaT.append(delta[j]/pow(Pr[liq],(1/3)))
        Cfx.append(a2[2]/math.sqrt((U*x[j])/v[liq]))
    deltas.append(delta)
    deltaTs.append(deltaT)
    Cfxs.append(Cfx)

plt.subplot2grid(gridSize, (2,0), rowspan = 1, colspan = 2)
for liq in range(0,3):
    plt.plot(x, deltas[liq], label=liquid[liq])
    plt.title(titles[0] + m[2])
    plt.xlabel(r'$x$')
    plt.ylabel(ylabels[0])
    plt.legend()
    
plt.subplot2grid(gridSize, (2,2), rowspan = 1, colspan = 2)
for liq in range(0,3):
    plt.plot(x, deltaTs[liq], label=liquid[liq])
    plt.title(titles[1] + m[2])
    plt.xlabel(r'$x$')
    plt.ylabel(ylabels[1])
    plt.legend()
    
plt.subplot2grid(gridSize, (4,1), rowspan = 1, colspan = 2)
for liq in range(0,3):
    plt.plot(x, Cfxs[liq], label=liquid[liq])
    plt.title(titles[2] + m[2])
    plt.xlabel(r'$x$')
    plt.ylabel(ylabels[2])
    plt.legend()

plt.tight_layout()
plt.savefig('./170329Report/all.eps')