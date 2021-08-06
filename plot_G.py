import matplotlib.pyplot as plt
import numpy as np

y, t = np.loadtxt ('tmp_mean.txt', unpack=True)

s = 1961-1901+1-1
e = 1990-1901+1-1
base = sum (t[s:e]) / 30.0

plt.plot(y,t-base)

plt.show( )