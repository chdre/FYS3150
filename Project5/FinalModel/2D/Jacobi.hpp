#include <iostream>
#include <cmath>
#include <armadillo>

using namespace  std;
using namespace  arma;

//void JSolver(mat A, mat &u, double dx, double dt, int n, double alpha);
//void JSolver(double e, double d, int n, mat &u);
void JSolver(double d, int n, double alpha, mat &u, mat rhs);
