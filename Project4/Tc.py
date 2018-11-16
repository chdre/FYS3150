import numpy as np
import matplotlib.pyplot as plt


file = "data/old2/L40-1e6.txt"
file2 = "data/old2/L60-1e6.txt"
file3 = "data/old2/L80-1e6.txt"
file4 = "data/old2/L100-1e6.txt"


chi, T = np.loadtxt(file, usecols=(3,4), unpack=True)
chi2 = np.loadtxt(file2, usecols=3)
chi3 = np.loadtxt(file3, usecols=3)
chi4 = np.loadtxt(file4, usecols=3)

files = [file,file2,file3,file4]

L = np.array([40,60,80,100])
Lx = np.linspace(0,120,1)

def Tcinf(Tc,L):
    return Tc - 1.0/L

ind = np.argmax(chi)
Tc = T[ind]

ind2 = np.argmax(chi2)
Tc2 = T[ind2]

ind3 = np.argmax(chi3)
Tc3 = T[ind3]

ind4 = np.argmax(chi4)
Tc4 = T[ind4]

plt.plot(Tcinf(Tc,L[0]),Lx, 'o')
plt.plot(Tcinf(Tc2,L[1]),Lx, 'o')
plt.plot(Tcinf(Tc3,L[2]),Lx, 'o')
plt.plot(Tcinf(Tc4,L[3]),Lx, 'o')

plt.show()
