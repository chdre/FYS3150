CPPflags = g++ -O3

LIB = -DARMA_DONT_USE_WRAPPER -lblas -llapack

main: p4.o
		${CPPflags} p4.cpp -o test.x  ${LIB}
