CPPflags = g++ -O3

LIB = -DARMA_DONT_USE_WRAPPER -lblas -llapack

3a: main.o euler.o planet.o
		${CPPflags} main.cpp planet.cpp euler.cpp  -o test.x  ${LIB}
