import matplotlib.pyplot as plt
import numpy as np
import os

file = "data/L20-T1-mcc1e4.txt"

E = np.loadtxt(file, usecols=0)
M = np.loadtxt(file, usecols=1)
# C_V = np.loadtxt(file, usecols=2, skiprows=2)
# chi = np.loadtxt(file, usecols=3, skiprows=2)
# Mabs = np.loadtxt(file, usecols=4, skiprows=2)
# temp = np.loadtxt("data/test.txt", usecols=4)

mcs = np.linspace(0,len(M),len(M))

plt.plot(mcs,E)
plt.legend(["Energy"], prop={'size': 15})
plt.title("Energy with T = 1 for 20x20 lattice", size=15)
plt.xlabel("Monte Carlo cycles", size=15); plt.ylabel("E", size=15)
plt.show()

#os.remove("data/test.txt")
