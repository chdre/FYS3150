import matplotlib.pyplot as plt
import numpy as np
import os

E = np.loadtxt("data/mean_vals.txt", usecols=0)
M = np.loadtxt("data/mean_vals.txt", usecols=1)
C_V = np.loadtxt("data/mean_vals.txt", usecols=2)
chi = np.loadtxt("data/mean_vals.txt", usecols=3)
# temp = np.loadtxt("data/mean_vals.txt", usecols=4)

mcs = np.linspace(0,len(M),len(M))

plt.plot(mcs,E)
plt.show()

os.remove("data/mean_vals.txt")
