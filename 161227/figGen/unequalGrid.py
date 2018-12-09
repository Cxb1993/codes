import numpy as np
import matplotlib as mpl
#mpl.use('PS')
mpl.rcParams['image.cmap'] = 'jet'
import matplotlib.pyplot as plt

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# Grid operations
indexer = 0
for file in [("../161227_heatTransfer2D/data/coordinatesXSimple.dat",
              "../161227_heatTransfer2D/data/coordinatesYSimple.dat",
              "../161227_report/gridSimple.eps"),
             ("../161227_heatTransfer2D/data/coordinatesXCase01.dat",
              "../161227_heatTransfer2D/data/coordinatesYCase01.dat",
              "../161227_report/gridCase01.eps"),
             ("../161227_heatTransfer2D/data/coordinatesXCase02.dat",
              "../161227_heatTransfer2D/data/coordinatesYCase02.dat",
              "../161227_report/gridCase02.eps")]:
    x = np.loadtxt(file[0])
    y = np.loadtxt(file[1])
    xIntervalSize = 1 / (x.shape[0] - 1)
    yIntervalSize = 1 / (y.shape[0] - 1)
    gridSize = (2, 3)
    plt.figure(figsize=(5.39749, 5.39749 * 2 / 3))
    
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
    if indexer == 0:
        plt.xlabel(r"$f(x)=L\sin\left(\frac{\pi}{L}\times\frac{i}{n_x}\right)$", fontsize = 10)
    elif indexer == 1:
        plt.xlabel(r"$f(x)=\frac{L}{2}\left[ 1 + \sin\left(\pi\left(\frac{i}{n_x} - \frac{1}{2}\right)\right)\right]$", fontsize = 10)
    elif indexer == 2:
        plt.xlabel(r"$f(x)=L\sin\left(\frac{\pi}{2}\times\frac{i}{n_x}\right)$", fontsize = 10)

    # y-axis graph
    plt.subplot2grid(gridSize, (1, 2))
    xValues = np.arange(0.0, 1.0 + yIntervalSize, yIntervalSize)
    plt.plot(xValues, y)
    plt.xticks(np.arange(2), ['$0$', r'$\frac{j_{ny}}{n_y}$'])
    plt.yticks(np.arange(2), ['$0$', '$W$'])
    plt.title("Function for y coordinates")
    if indexer == 0:
        plt.xlabel(r"$f(y)=W\sin\left(\frac{\pi}{W}\times\frac{j}{n_y}\right)$", fontsize = 10)
    elif indexer == 1:
        plt.xlabel(r"$f(y)=W\sin\left(\frac{\pi}{2}\times\frac{j}{n_y}\right)$", fontsize = 10)
    elif indexer == 2:
        plt.xlabel(r"$f(y)=W\sin\left(\frac{\pi}{2}\times\frac{j}{n_y}\right)$", fontsize = 10)
    
    plt.tight_layout()
    plt.savefig(file[2])
    indexer += 1
