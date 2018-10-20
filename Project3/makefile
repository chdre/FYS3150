CPPflags = g++ -O3

LIB = -DARMA_DONT_USE_WRAPPER -lblas -llapack

main: main.o planet.o euler.o
		${CPPflags} main.cpp planet.cpp euler.cpp  -o test.x  ${LIB}
