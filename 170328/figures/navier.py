import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from PIL import Image

mpl.rcParams['image.cmap'] = 'jet'

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

lm = 25.0

x = np.linspace(0.0, 0.010, 48)
y = np.linspace(0.0, 0.005, 48)

#q = plt.quiver(X,Y,u,v,angles='xy',scale=1000,color='r')

# Plotting the numerical solution
plt.figure(figsize=(5.39749, 5.39749))

cases = ["0", "0.0005", "0.01001", "0.01999"]
titles = [r"$0$", r"$0.0005$", r"$0.01001$", r"$0.01999$"]
sections = [(0,0), (0,1), (1,0), (1,1)]

for j in range(0, 4):
    plt.subplot2grid((2, 2), sections[j])
    fileInput = cases[j];
    u = np.loadtxt("../data/Navier/U_"
                   + fileInput + ".dat")
    v = np.loadtxt("../data/Navier/V_"
                   + fileInput + ".dat")
    q = plt.quiver(x, y, u, v, angles='xy', scale=50, color='b')
    plt.title(r"$t = $" + titles[j] + r" s")

plt.tight_layout()
plt.savefig("flow.eps")

# ims = []
# for j in range(0,4):
#     fileInput = cases[j];
#     u = np.loadtxt("../TurbulentCPP/data/U_"
#                    + fileInput + ".dat")
#     v = np.loadtxt("../TurbulentCPP/data/V_"
#                    + fileInput + ".dat")
#     q = plt.quiver(x, y, u, v, angles='xy', scale=1000, color='b')
#     plt.title(r"$t = $" + titles[j] + r" s")
#     plt.tight_layout()
#     plt.savefig("flow" + fileInput + ".jpg")
#     im = plt.imshow(Image.open("flow" + fileInput + ".jpg"), animated=True)
#     ims.append([im])
#
# fig = plt.figure()
# ani = animation.ArtistAnimation(fig, ims, interval=50, blit=True,
#                                 repeat_delay=1000)
# plt.show()
# ani.save('dynamic_images.mp4')
