import matplotlib.pyplot as plt
import numpy as np

energy = np.loadtxt("data/mean_vals.txt", usecols=0)
energy2 = np.loadtxt("data/mean_vals.txt", usecols=1)
magMoment = np.loadtxt("data/mean_vals.txt", usecols=2)
magMoment2 = np.loadtxt("data/mean_vals.txt", usecols=3)
# temp = np.loadtxt("data/mean_vals.txt", usecols=4)

print(np.size(energy))

C_V = (energy2 - energy**2)
chi = (magMoment2 - magMoment**2)

mcs = np.linspace(0,len(energy),len(energy))
print(len(mcs))

plt.plot(mcs,energy)
plt.show()
