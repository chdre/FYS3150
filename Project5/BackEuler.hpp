#include <iostream>
#include <cmath>
#include <armadillo>
#include "tridiag.hpp"

using namespace  std;
using namespace  arma;

void FWSolver(int n, double alpha, int tmax);
void tridiagonal(double d, double e, vec &u, vec rhs, int n);
