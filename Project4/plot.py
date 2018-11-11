import matplotlib.pyplot as plt
import numpy as np
import os

file = "data/test.txt"

E = np.loadtxt(file, usecols=0)
M = np.loadtxt(file, usecols=1)
C_V = np.loadtxt(file, usecols=2)
chi = np.loadtxt(file, usecols=3)
Mabs = np.loadtxt(file, usecols=4)
# temp = np.loadtxt("data/test.txt", usecols=4)

mcs = np.linspace(0,len(M),len(M))

plt.plot(mcs,M)
plt.show()

os.remove("data/test.txt")
