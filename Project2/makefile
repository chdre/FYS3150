CPPflags = g++ -O3

LIB = -DARMA_DONT_USE_WRAPPER -lblas -llapack

#main: final_2b.x
#			${CPPflags} -o final_2b.x 2.cpp eigvals.cpp ${LIB}

test: test-2.x test-2.cpp
			${CPPflags} -o test-2.x test-2.cpp eigvals.cpp ${LIB}
