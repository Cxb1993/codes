import numpy as np
import matplotlib as mpl
mpl.use('PS')
import matplotlib.pyplot as plt


# Data loading
t_01 = np.loadtxt('D:\\PhD\\Tasks\\005\\heatTransfer1D\\data\\ht1dn158t0.1.dat')
t_1 = np.loadtxt('D:\\PhD\\Tasks\\005\\heatTransfer1D\\data\\ht1dn158t1.dat')
t_10 = np.loadtxt('D:\\PhD\\Tasks\\005\\heatTransfer1D\\data\\ht1dn158t10.dat')
t_100 = np.loadtxt('D:\\PhD\\Tasks\\005\\heatTransfer1D\\data\\ht1dn158t100.dat')
t_1000 = np.loadtxt('D:\\PhD\\Tasks\\005\\heatTransfer1D\\data\\ht1dn158t1000.dat')
t_2100 = np.loadtxt('D:\\PhD\\Tasks\\005\\heatTransfer1D\\data\\ht1dn158t2100.dat')
t_3000 = np.loadtxt('D:\\PhD\\Tasks\\005\\heatTransfer1D\\data\\ht1dn158t3000.dat')
t_3500 = np.loadtxt('D:\\PhD\\Tasks\\005\\heatTransfer1D\\data\\ht1dn158t3500.dat')

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.figure(figsize=(5.39749, 5.39749))

plt.plot(t_01[:,0], t_01[:,1], label='t = 0.1')
plt.plot(t_1[:,0], t_1[:,1], label='t = 1')
plt.plot(t_10[:,0], t_10[:,1], label='t = 10')
plt.plot(t_100[:,0], t_100[:,1], label='t = 100')
plt.plot(t_1000[:,0], t_1000[:,1], label='t = 1000')
plt.plot(t_2100[:,0], t_2100[:,1], label='t = 2100')
plt.plot(t_3000[:,0], t_3000[:,1], label='t = 3000')
plt.plot(t_3500[:,0], t_3500[:,1], label='t = 3500')

plt.legend(loc=0)
plt.xlabel(r'Distance, $x$', fontsize=12)
plt.ylabel(r'Temperature ($^{\circ}$C)',fontsize=12)

plt.savefig('D:\\PhD\\Tasks\\005\\report161215\\ht1dn158.eps')