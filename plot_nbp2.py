import matplotlib.pyplot as plt
import numpy as np

# (width, height) inches
#plt.figure(figsize=(12,12))

#xright = 1750
#ytop   =    1
#plt.xlim(0,xright)
#plt.ylim(0,ytop)

Year, FF, LU, ATM, OS, LS, CEM, IMB = np.loadtxt ('GCB.txt', unpack=True)
yr_inv, inv_flux = np.loadtxt ('inv_global.txt', delimiter=',', unpack=True)
kyr3, hyr3, NPP3, Rh3, NEE3, B3, SOM3 = np.loadtxt ('/home/adf10/HYBRID10/RUN1/SAVE3/global_means00600.txt', unpack=True, skiprows=1)
kyr4, hyr4, NPP4, Rh4, NEE4, B4, SOM4 = np.loadtxt ('/home/adf10/HYBRID10/RUN1/SAVE4/global_means00640.txt', unpack=True, skiprows=1)
kyr5, hyr5, NPP5, Rh5, NEE5, B5, SOM5 = np.loadtxt ('/home/adf10/HYBRID10/RUN1/SAVE5/global_means00670.txt', unpack=True, skiprows=1)
kyr6, hyr6, NPP6, Rh6, NEE6, B6, SOM6 = np.loadtxt ('/home/adf10/HYBRID10/RUN2/SAVE12/global_means04340.txt', unpack=True, skiprows=1)
kyr7, hyr7, NPP7, Rh7, NEE7, B7, SOM7 = np.loadtxt ('/home/adf10/HYBRID10/RUN2/SAVE13/global_means04380.txt', unpack=True, skiprows=1)
kyr8, hyr8, NPP8, Rh8, NEE8, B8, SOM8 = np.loadtxt ('/home/adf10/HYBRID10/RUN2/SAVE14/global_means04420.txt', unpack=True, skiprows=1)
kyr9, hyr9, NPP9, Rh9, NEE9, B9, SOM9 = np.loadtxt ('/home/adf10/HYBRID10/RUN3/SAVE3/global_means06380.txt', unpack=True, skiprows=1)
kyr10, hyr10, NPP10, Rh10, NEE10, B10, SOM10 = np.loadtxt ('/home/adf10/HYBRID10/RUN3/SAVE4/global_means06420.txt', unpack=True, skiprows=1)

vyr, vNBP = np.loadtxt ('output.txt', unpack=True)

plt.xlim(1958, 2026)
plt.ylim(-4, 4)

L = FF - ATM - OS
plt.plot(Year, L)

a = 0.5
# b = 0.04 # * b * (hyr3-1950)
#plt.plot(hyr3, a*NEE3, 'r')
#plt.plot(hyr4, a*NEE4, 'r')
#plt.plot(hyr5, a*NEE5, 'r')
#plt.plot(hyr6, a*NEE6, '--')
#plt.plot(hyr7, a*NEE7, '--')
#plt.plot(hyr8, a*NEE8, '--')
plt.plot(hyr9, a*NEE9, 'y') #done before
plt.plot(hyr10, a*NEE10,'y') #done before
#b = 0.03
#plt.plot(hyr9, a*NEE9 * b * (hyr9-1950), 'r')
#plt.plot(hyr10, a*NEE10 * b * (hyr10-1950),'r')

#plt.plot(hyr4, a*(NPP4-110.0), 'm')
#plt.plot(hyr5, a*(NPP5-110.0), 'm')

plt.plot(yr_inv, inv_flux)

offset = vyr * 7 / (2020 - 1975) - 307.0
plt.plot(vyr, vNBP + offset)

#b = 0.03
#plt.plot(hyr3, a*NEE3 * b * (hyr3-1950), 'r')
#plt.plot(hyr4, a*NEE4 * b * (hyr4-1950), 'r')
#plt.plot(hyr5, a*NEE5 * b * (hyr5-1950), 'r')

plt.show( )