CPPflags = g++ -O3

LIB = -DARMA_DONT_USE_WRAPPER -lblas -llapack

#2b: 2b.o eigvals.o
#		${CPPflags} -o final_2b.x 2b.cpp eigvals.cpp ${LIB}

#test: test-2.x test-2.cpp
#		${CPPflags} -o test-2.x test-2.cpp eigvals.cpp ${LIB}

2d: 2d.o eigvals.o
		${CPPflags} -o test2d.x 2d.cpp eigvals.cpp ${LIB}
