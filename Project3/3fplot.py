import matplotlib.pyplot as plt
import numpy as np

x = np.loadtxt("data/testdata3f.txt")
t = np.linspace(0,100,415)

plt.plot(t,x)
plt.show()
