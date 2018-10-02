#include <armadillo>
#include <cmath>
#include <iostream>
#include <fstream>
#include "eigvals.h"

using namespace std;
using namespace arma;

void write_to_file(mat &A, mat &R, double a, double d, double h, int n){
        /* Writing values to file. Takes the adress of matrix A and R, and the
           values along the diagonal, upper diagonal and lower diagonal of A, theta
           step size h and matrix dimension n.*/
        vec eigvals_arm = eigvals_arma(n, h, a, d); // eigenvalues from Armadillo
        vec eigvals_ana = eigvals_analytical(n, a, d); // analytical eigenvalues
        vec eigvals_jac = jacobi(A, R, n, h);

        // writing to file
        ofstream outfile;
        outfile.open(filename);
        for (int i = 0; i < n; i++) {
                // writing all eigenvalues to file
                outfile << eigvals_jac[i];
                outfile << " " << eigvals_arm[i];
                outfile << " " << eigvals_ana[i] << endl;
        }
        outfile.close();
}

main(int argc, char *argv[]) {
        string filename;
        int n;
        if(argc <= 1) {
                cout << "Missing arguments:" << argv[0] << " specify output filename and value of n" << endl;
                exit(1);
        }
        else{
                // reading filename and value of n from command line
                filename = argv[1];
                n = atoi(argv[2]);      // setting n from ascii to integer
        }

        double rho_0 = 0.0;   // start value of dimensionless length variable
        double rho_n = 1.0;   // end value of dimensionless length variable
        double h = (rho_n - rho_0)/double(n);  // step length
        double a = -1.0/pow(h,2);
        double d = 2.0/pow(h,2);

        // Creating tridiagonal matrix
        mat A = zeros<mat>(n, n);        // matrix for A
        A.diag() += d; // center diagonal
        A.diag(1) += a; // upper diagonal
        A.diag(-1) += a; // lower diagonal

        mat R;  // matrix of eigenvectors

        write_to_file(A,R,a,d,h,n);

        return 0;
}
