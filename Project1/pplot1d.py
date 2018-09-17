import matplotlib.pyplot as plt
import sys
import numpy as np

infile = sys.argv[1]

def y(h):
    h = np.array(h)
    return 2*h

log10h = []
eps = []

for lines in open(infile, "r"):
    col = lines.split()
    eps.append(float(col[0]))
    log10h.append(float(col[1]))

plt.figure()
plt.title('Logarithmic plot of error as a function of step size h',size=18)
plt.plot(log10h,eps)
plt.plot(log10h,y(log10h))
plt.xlabel('$\\log_{10}(h)$', size=18);plt.ylabel('$\\epsilon_i$',size=18)
plt.legend(['$log_{10}(\\epsilon_i)$', '$2log_{10}(h)$'],prop={'size': 18})
plt.show()

eps = 10**(np.array(eps))
print(eps)
