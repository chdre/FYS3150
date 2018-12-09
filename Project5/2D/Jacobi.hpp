#include <iostream>
#include <cmath>
#include <armadillo>

using namespace  std;
using namespace  arma;

void JSolver(mat &u, mat &rhs, double dx, double dt, double tol, int n, double alpha);
