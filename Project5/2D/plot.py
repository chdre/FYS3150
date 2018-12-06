import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

Euler = True

j = 10
j2 = 30

if Euler == True:
    u = np.loadtxt("FWEuler2D.txt")
    n = len(u[0,:])

    vals = []

    t_steps = int(len(u[:,0])/n)

    index = 0 # dummy index to add matrices to vals array

    for t in range(t_steps):
        start = n*t
        end = start+n
        #print(start, end, t, t_steps)
        #print(u[start:end])

        vals.append(u[start:end])

        index += 1

    vals = np.array(vals)

    print(vals[-1][0], vals[-1][2])

    x = np.linspace(0,1,n)

    plt.plot(vals[-1])
    plt.show()
