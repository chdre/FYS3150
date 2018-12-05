CPPflags = g++ -Wall -std=c++11


LIB = -DARMA_DONT_USE_WRAPPER -lblas -llapack


forw: main.o ForwEuler.o
			${CPPflags} main.cpp ForwEuler.cpp -o test.x  ${LIB}

back: main.o BackEuler.o tridiag.o
			${CPPflags} main.cpp BackEuler.cpp tridiag.cpp -o test.x  ${LIB}

cn: main.o CrankNic.o tridiag.o
			${CPPflags} main.cpp CrankNic.cpp tridiag.cpp -o test.x  ${LIB}
