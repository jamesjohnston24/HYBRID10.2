import matplotlib.pyplot as plt
import numpy as np

y, t = np.loadtxt ('tmp_mean.txt', unpack=True)

plt.plot(y,t)
what
plt.show( )