#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using namespace std;
using namespace arma;

main(int argc, char *argv[]) {
  int n; // size of matrix

  if (argc <= 1) {
    cout << "Missing arguments:" << argv[0] << " specify size of matrix N"
         << endl;
    exit(1);
  } else {
    // reading value of N
    n = atoi(argv[1]); // N from asci
  }

  double h = 1.0 / (n + 1); // step length
  double eps = 1.0e-8;      // tolerance

  // Creating tridiagonal matrix
  mat A = zeros<mat>(n + 2, n + 2);      // matrix for A
  A.diag() += double(2.0) / pow(h, 2);   // center diagonal
  A.diag(1) -= double(1.0) / pow(h, 2);  // upper diagonal
  A.diag(-1) -= double(1.0) / pow(h, 2); // lower diagonal

  // creating a while loop that checks whether the off diagonal elements are
  // larger than eps
  while (maxelm > eps) {
  }

  double max_offdiag_test(double *maxelm, int *k, int *l) {
    // function to find the maximum value of
    mat A_temp = A.diag() *= 0.0;         // diagonal = 0
    uword max_index = A_temp.index_max(); // finding index of max element
    *maxelm = A_temp(max_index);          // max element
  }

  double max_offdiag(double *maxelm, ) {
    *maxelm = A[1, 1];
    for (int i = 1; i < n; i++) {
      for (int j = 1; j < n; j++) {
        if (abs(A[i, j]) > maxelm) { // checking whether
          *maxelm = abs(A[i, j]);
        }
      }
    }
  }
  return 0;
}
