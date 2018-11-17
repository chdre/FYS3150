import numpy as np

n = 2
T = 1.0

E = -8.0*np.sinh(8/T)/(3.0 + np.cosh(8/T))
#CV = 1/T**2*(64.0*np.cosh(8/T)*(3 + np.cosh(8/T)) - 64*np.sinh(8/T)**2)/(3+np.cosh(8/T))**2
CV = (64*np.cosh(8.0/T))/(T**2*(3 + np.cosh(8.0/T)))
Mabs = (2*np.exp(8/T) + 4.0)/(3.0 + np.cosh(8/T))
M2 = 8.0*(np.exp(8/T) + 1.0)/(3.0 + np.cosh(8/T))
chi = M2/T

print(E, Mabs, CV, chi)
