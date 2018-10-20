import matplotlib.pyplot as plt
import sys

infile = sys.argv[1]

x = []
y = []
x2 = []
y2 = []

for lines in open(infile, "r"):
    col = lines.split()
    x.append(float(col[0]))
    y.append(float(col[1]))
    x2.append(float(col[2]))
    y2.append(float(col[3]))

plt.figure()
plt.plot(x,y)
plt.plot(x2,y2)
plt.legend(['Earth','SuNN'],prop={'size': 18})
plt.show()
