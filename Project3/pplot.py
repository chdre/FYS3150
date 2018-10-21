import matplotlib.pyplot as plt
import sys


x = np.loadtxt("test.txt", usecols=0)
y = np.loadtxt("test.txt", usecols=1)


plt.figure()
plt.plot(x,y)
plt.title("Earth orbit around sun using Velocity Verlet",size=15)
plt.xlabel("x",size=15); plt.ylabel("y",size=15)
plt.legend(['Eearth'],prop={'size': 15})
plt.show()
