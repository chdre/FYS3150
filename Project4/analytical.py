import numpy as np

n = 2
T = 1.0

E = -8.0*np.sinh(8/T)/(3.0 + np.cosh(8/T))
CV = 64/(T*T) * (np.cosh(8)*(3 + np.cosh(8)) - (np.sinh(8)**2))/(3 + np.cosh(8))**2
Mabs = (2*np.exp(8/T) + 4.0)/(3.0 + np.cosh(8/T))
M2 = 8.0*(np.exp(8/T) + 1.0)/(3.0 + np.cosh(8/T))
chi = M2/T

print(E, Mabs, CV, chi)
