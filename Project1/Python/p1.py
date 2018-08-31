import numpy as np

n = 1000    #Number of points
h = 1./(n+1)

a_diag = 2.*np.ones(n+1)
a_lowd = -1.*np.ones(n+1)
a_upd = -1.*np.ones(n+1)

v = np.ones(n+1)

btde = h**2*np.zeros(n+1)
