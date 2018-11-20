CC = gcc
CXX = mpic++
CPPflags = -g -Wall -std=c++11 -O3

LIB = -DARMA_DONT_USE_WRAPPER -lblas -llapack

main: main.o
		${CXX} ${CPPflags} main.cpp -o main.x  ${LIB}
