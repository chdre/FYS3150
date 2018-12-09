CPPflags = g++ -Wall -std=c++11 -O3


LIB = -DARMA_DONT_USE_WRAPPER -lblas -llapack


forw: main.o ForwEuler.o
			${CPPflags} main.cpp ForwEuler.cpp -o test.x  ${LIB}

back: main.o BackEuler.o tridiag.o
			${CPPflags} main.cpp BackEuler.cpp tridiag.cpp -o test.x  ${LIB}

cn: main.o CrankNic.o tridiag.o
			${CPPflags} main.cpp CrankNic.cpp tridiag.cpp -o test.x  ${LIB}

an1d: main.o analytical.o
			${CPPflags} main.cpp analytical.cpp -o test.x  ${LIB}

all: main.o ForwEuler.o BackEuler.o CrankNic.o tridiag.o analytical.o
			${CPPflags} main.cpp ForwEuler.cpp BackEuler.cpp  analytical.cpp CrankNic.cpp tridiag.cpp -o test.x  ${LIB}
