#include <iostream>
#include <cmath>
#include <armadillo>

using namespace  std;
using namespace  arma;

double max_offdiag(mat &A, int n, int *l, int *k);
void rotate(mat &A, mat &R, int &k, int &l, int n);
vec eigvals_arma(int n, double h, double a, double d);
vec eigvals_analytical(int n, double a, double d);
vec jacobi(mat &A, mat &R, int n, double h);
