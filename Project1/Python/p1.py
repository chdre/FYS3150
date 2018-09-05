import numpy as np
import sys
import matplotlib.pyplot as plt

n = int(sys.argv[1])    #Number of points
h = 1./(n+1)

#Diagonal vectors forming matrix A
a = np.zeros(n-1)     #Diagonal
b = np.zeros(n)   #Upper diagonal
c = np.zeros(n-1)   #Lower diagonal

u = np.zeros(n)     #Vector with v_i's.
bt = np.zeros(n)    #Vector b~
ft = np.zeros(n)

def matrix_fill(value, i, matrix):
    matrix[i] = value

#Boundary conditions
u[0] = 0
u[-1] = u[0]

x = np.linspace(0,1,n)      #x_i's
f = 100.*np.exp(-10.*x)     #f_i's
fprime = h**2*f             #RHS of equation, f~

#Filling lower and upper diagonals (vectors)
a.fill(-1); c.fill(-1); b.fill(2)      #Value of coefficient of v_{i+1} and v_{i-1}.
bt[0] = b[0]
ft[0] = fprime[0]

for i in range(n-2):
    bt[i+1] = b[i+1] - a[i]*c[i]/bt[i]
    ft[i+1] = fprime[i+1] - ft[i]*a[i]/bt[i]

#Equations
for i in range(n-2,0,-1):
    u[i] = (ft[i] - c[i]*u[i+1])/bt[i]

cc_sol = 1. - (1. - np.exp(-10.))*x - np.exp(-10.*x)
cc_sol = cc_sol

plt.plot(x,cc_sol)
plt.plot(x,u)
plt.legend(['CC', 'u'])
plt.show()

print(cc_sol)












#
