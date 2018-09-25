#include <armadillo>
#include <cmath>
#include <iostream>
#include <fstream>
#include "eigvals.h"


using namespace std;
using namespace arma;

//-DARMA_DONT_USE_WRAPPER -lblas -llapack

void write_to_file(mat &A, mat &R, double a, double d, double h, int n){
        vec eigvals_arm = eigvals_arma(n, h, a, d); // eigenvalues from Armadillo
        vec eigvals_ana = eigvals_analytical(n, a, d); // analytical eigenvalues
        vec eigvals_jac = jacobi(A, R, n, h);

        // writing to file
        ofstream outfile;
        outfile.open("2b_results.txt");
        for (int i = 0; i < n; i++) {
                outfile << eigvals_jac[i];
                outfile << " " << eigvals_arm[i];
                outfile << " " << eigvals_ana[i] << endl;
        }
        outfile.close();
}

main(int argc, char *argv[]) {
        int n = atoi(argv[1]); // size of matrix

        double h = 1.0/double(n); // step length, preserving u(L) = 1
        double a = -1.0/pow(h,2);
        double d = 2.0/pow(h,2);

        // Creating tridiagonal matrix
        mat A = zeros<mat>(n, n);        // matrix for A
        A.diag() += d; // center diagonal
        A.diag(1) += a; // upper diagonal
        A.diag(-1) += a; // lower diagonal

        mat R;  // matrix of eigenvectors

        return 0;
}
