import matplotlib.pyplot as plt
import numpy as np

kinetic = np.loadtxt("data/test.txt", usecols=0)
potential = np.loadtxt("data/test.txt", usecols=1)
time = np.loadtxt("data/test.txt", usecols=2)

plt.figure()
plt.plot(time,kinetic)
plt.plot(time,potential)

plt.title("Energy of Earth in orbit around sun",size=15)
plt.xlabel("time [years]",size=15); plt.ylabel("Energy",size=15)
plt.legend(['Kinetic', 'Potential'],prop={'size': 15})
plt.show()
