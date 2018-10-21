CPPflags = g++ -O3

LIB = -DARMA_DONT_USE_WRAPPER -lblas -llapack

main: main.o planet.o solver.o
		${CPPflags} main.cpp planet.cpp solver.cpp  -o test.x  ${LIB}
