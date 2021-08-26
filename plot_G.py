import matplotlib.pyplot as plt
import numpy as np

#y, t = np.loadtxt ('tmp_mean.txt', unpack=True)

Year, FF, LU, ATM, OS, LS, CEM, IMB = np.loadtxt ('GCB.txt', unpack=True)
y, b, s, nbp = np.loadtxt ('transient.txt', unpack=True)
yr_inv, inv_flux = np.loadtxt ('inv_global.txt', delimiter=',', unpack=True)
kyr10, hyr10, NPP10, Rh10, NEE10, B10, SOM10 = np.loadtxt ('/home/jhj34/HYBRID10/RUN3/SAVE4/global_means06420.txt', unpack=True, skiprows=1)

plt.xlim(1958, 2025)
#s = 1961-1901+1-1
#e = 1990-1901+1-1
#base = sum (t[s:e]) / 30.0
#print (s, e, base)

#plt.plot(y,t)
a = 0.0
b = 0.5
y = y + 1900
L = FF - ATM - OS
plt.plot(Year, L, 'g')
plt.plot (y, a+ b*nbp, 'r')
plt.plot(yr_inv, inv_flux, 'b')
#plt.plot(hyr10, 0.5*NEE10,'y')

plt.show( )