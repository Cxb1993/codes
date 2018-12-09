import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams['image.cmap'] = 'jet'

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.figure(figsize=(5.39749, 5.39749))

analytical = np.loadtxt("../data/laplaceAnalyticalError.dat")
plt.plot(analytical[1][:], analytical[0][:], '.k', label='Analytical')

files = ['0.01', '0.001', '0.0001', '1e-005']
diffs = [r'$1 \times 10^{-2}$', r'$1 \times 10^{-3}$', r'$1 \times 10^{-4}$',
          r'$1 \times 10^{-5}$']
for s in range(0,4):
    data = np.loadtxt("../data/laplaceSORError" + files[s] + ".dat")
    plt.plot(data[1][:], data[0][:], label='Numerical: difference = '
             + diffs[s])

plt.title("Profile of $\phi$ along $y$-axis at $x = 0.5102$")
plt.xlabel("y-axis values")
plt.ylabel("$\phi$")
plt.legend()
plt.tight_layout()
plt.savefig("../170323_report/sorError.eps")