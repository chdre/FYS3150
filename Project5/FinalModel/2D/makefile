CPPflags = g++ -std=c++11 -O3


LIB = -DARMA_DONT_USE_WRAPPER -lblas -llapack


main: main.o ForwEuler.o BackEuler.o CrankNic.o tridiag.o analytical.o
			${CPPflags} main.cpp ForwEuler.cpp BackEuler.cpp CrankNic.cpp tridiag.cpp analytical.cpp -o test.x  ${LIB}

#If need to run methods below, please remove calls in main.cpp to functions not
#relevant for the specific method
forw: main.o ForwEuler.o
			${CPPflags} main.cpp ForwEuler.cpp -o test.x  ${LIB}

back: main.o BackEuler.o Jacobi.o
			${CPPflags} main.cpp BackEuler.cpp Jacobi.cpp -o test.x  ${LIB}

an2d: main.o analytical.o
			${CPPflags} main.cpp analytical.cpp -o test.x  ${LIB}
