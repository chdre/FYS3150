#include <armadillo>
#include <cmath>
#include <iostream>
#include <fstream>
#include "eigvals.h"

using namespace std;
using namespace arma;

main(int argc, char *argv[]) {
        int n = atoi(argv[1]); // size of matrix

        double rho_0 = 0.0;
        double rho_n = 1.0;
        double h = (rho_n - rho_0)/double(n);               // step length, preserving u(L) = 1
        double a = -1.0/pow(h,2);
        double d = 2.0/pow(h,2);

        vec ivec = linspace<vec>(1, n, n);
        vec rho = rho_0 + ivec*h;
        cout << rho << endl;

        // Creating tridiagonal matrix
        mat A = zeros<mat>(n, n);        // matrix for A
        A.diag() += d; // center diagonal
        A.diag(1) += a; // upper diagonal
        A.diag(-1) += a; // lower diagonal

        mat R;  // matrix of eigenvectors

        return 0;
}
