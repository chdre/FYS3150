import sys

N = int(sys.argv[1])

Sup = 0

for n in range(1,N+1):
    Sup += 1.0/n

Sdown = 0

for n in range(N,1+1,-1):
    Sdown += 1.0/n

print Sup, Sdown
