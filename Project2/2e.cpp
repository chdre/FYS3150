#include <armadillo>
#include <cmath>
#include <iostream>
#include <fstream>
#include "eigvals.h"

using namespace std;
using namespace arma;

void write_to_file(mat &A, mat &R, double a, double d, double h, int n, double w_r){
        /* Writing values to file. Takes the adress of matrix A and R, and the
           values along the diagonal, upper diagonal and lower diagonal of A, theta
           step size h and matrix dimension n.*/
        vec eigvals_jac = jacobi(A, R, n, h);

        // writing to file
        ofstream outfile;
        outfile.open("2e_results.txt");
        for (int i = 0; i < n; i++) {
                // writing eigenvalues of the jacobis method
                outfile << eigvals_jac[i] << endl;
        }
        outfile.close();
}

main(int argc, char *argv[]) {
        string filename;
        int n;
        double w_r, rho_n;

        if(argc <= 1) {
                cout << "Missing arguments:" << argv[0] << " specify output filename, value of n and omega_r" << endl;
                exit(1);
        }
        else{
                // reading filename, value of n, omega_r and rho_n from command line
                filename = argv[1];
                n = atoi(argv[1]); // size of matrix
                w_r = atof(argv[2]); // angular frequency
                rho_n = atof(argv[3]);  // end value of rho
        }

        double rho_0 = 0.0; // start value of dimensionless length variable
        double h = (rho_n - rho_0)/double(n);  // step length
        double a = -1.0/pow(h,2);
        double d = 2.0/pow(h,2);

        vec rho = rho_0 + linspace<vec>(1,n,n)*h;

        // Creating tridiagonal matrix
        mat A = zeros<mat>(n, n);        // matrix for A
        A.diag() += d + pow(w_r,2)*pow(rho,2) + 1.0/rho; // center diagonal
        A.diag(1) += a; // upper diagonal
        A.diag(-1) += a; // lower diagonal

        mat R;  // matrix of eigenvectors

        write_to_file(A, R, a, d, h, n, w_r);

        return 0;
}
