CPPflags = g++ -Wall -std=c++11 -O3


LIB = -DARMA_DONT_USE_WRAPPER -lblas -llapack


main: main.o ForwEuler.o
			${CPPflags} main.cpp ForwEuler.cpp -o test.x  ${LIB}
