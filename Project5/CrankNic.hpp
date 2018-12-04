#include <iostream>
#include <cmath>
#include <armadillo>
#include "tridiag.hpp"

using namespace std;
using namespace arma;
ofstream outfile;

void CNSolver(int n, double alpha, int tmax);
void tridiagonal(double d, double e, vec &u, vec rhs, int n);
