import matplotlib.pyplot as plt
import numpy as np
import os

file = "data/test.txt"

E = np.loadtxt(file, usecols=0, skiprows=2)
M = np.loadtxt(file, usecols=1, skiprows=2)
# C_V = np.loadtxt(file, usecols=2, skiprows=2)
# chi = np.loadtxt(file, usecols=3, skiprows=2)
# Mabs = np.loadtxt(file, usecols=4, skiprows=2)
# temp = np.loadtxt("data/test.txt", usecols=4)

mcs = np.linspace(0,len(M),len(M))

plt.plot(mcs,E)
plt.legend(["Energy"])
plt.xlabel("Monte Carlo cycles"); plt.ylabel("$\\langle E \\rangle$")
plt.show()

os.remove("data/test.txt")
