import matplotlib.pyplot as plt
import numpy as np

#y, t = np.loadtxt ('tmp_mean.txt', unpack=True)
y, b, s, nbp = np.loadtxt ('transient.txt', unpack=True)
yr_inv, inv_flux = np.loadtxt ('inv_global.txt', delimiter=',', unpack=True)

#s = 1961-1901+1-1
#e = 1990-1901+1-1
#base = sum (t[s:e]) / 30.0
#print (s, e, base)

#plt.plot(y,t)
b = 0.5
y = y + 1900
plt.plot (y, b*nbp)
plt.plot(yr_inv, inv_flux)

plt.show( )