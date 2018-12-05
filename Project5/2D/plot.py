import numpy as np
import matplotlib.pyplot as plt

Euler = True
Back = True
CN = True

j = 103
j2 = 30

if Euler == True:
    u = np.loadtxt("FWEuler.txt")
    n = len(u[0,:])
    y_length = len(u[:,0]) # length of y
    timesteps = y_length/n

    vals = np.zeros(timesteps)

    max = n
    for i in range(n):
        vals = u[:max,:max]
        max += n

    x = np.linspace(0,1,n)

    plt.plot(x,u[j,:])
    plt.plot(x,u[j2,:])
    plt.show()
