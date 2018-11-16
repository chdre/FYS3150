import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg


file = "data/L40-1e6.txt"
file2 = "data/L60-1e6.txt"
file3 = "data/L80-1e6.txt"
file4 = "data/L100-1e6.txt"


chi, T = np.loadtxt(file, usecols=(3,4), unpack=True)
chi2 = np.loadtxt(file2, usecols=3)
chi3 = np.loadtxt(file3, usecols=3)
chi4 = np.loadtxt(file4, usecols=3)

files = [file,file2,file3,file4]

L = np.array([1.0/100, 1.0/80, 1.0/60, 1.0/40])
T = np.array([T[np.argmax(chi4)] - L[3], T[np.argmax(chi3)]- L[2], T[np.argmax(chi2)]- L[1], T[np.argmax(chi)]- L[0]])

z = np.polyfit(L,T,1)
zz = np.poly1d(z)

print(zz[0], zz[1])

plt.plot(L, T,'o')
plt.plot(L, zz(L))
plt.legend(['Data', 'Fitted'], prop={'size':15})
plt.title('Linear regression of critical temperature as a function of L',size=15)
plt.xlabel('1/L', size=15); plt.ylabel('$T_C$ [k/J]')


plt.show()
