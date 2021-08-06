import matplotlib.pyplot as plt
import numpy as np

#y, t = np.loadtxt ('tmp_mean.txt', unpack=True)
y, b, s, nbp = np.loadtxt ('transient.txt', unpack=True)
yr_inv, inv_flux = np.loadtxt ('inv_global.txt', delimiter=',', unpack=True)
kyr10, hyr10, NPP10, Rh10, NEE10, B10, SOM10 = np.loadtxt ('/home/adf10/HYBRID10/RUN3/SAVE4/global_means06420.txt', unpack=True, skiprows=1)

plt.xlim(1958, 2025)
#s = 1961-1901+1-1
#e = 1990-1901+1-1
#base = sum (t[s:e]) / 30.0
#print (s, e, base)

#plt.plot(y,t)
a = 0.0
b = 1.0
y = y + 1900
plt.plot (y, b*nbp)
plt.plot(yr_inv, inv_flux)
plt.plot(hyr10, b*NEE10,'y')

plt.show( )