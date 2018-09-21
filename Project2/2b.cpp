#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using namespace std;
using namespace arma;

double max_offdiag_test(mat &A, int n, double h) {
  /* function to find the maximum value of the array A. The diagonal is set to
   0. The index of the maximum value is found, since the function returns a
   single number we must transform this to indices for the matrix.*/
  mat A_temp = abs(A);
  A_temp -= eye(size(A)) * (double(2.0) / pow(h, 2)); // diagonal = 0
  A_temp(0, 1) = 0;
  uword max_index = A_temp.index_max(); // finding index of max element
  int l = int(max_index) / n;           // find correct index for column
  int k = int(max_index) - (n * l);     // find correct index for row
  double maxelm = A_temp(k, l);         // max element
  // cout << A_temp << endl;    // printing array
  // cout << k << " " << l << endl; // printing indices
  return maxelm;
}

main(int argc, char *argv[]) {
  int n; // size of matrix

  if (argc <= 1) {
    cout << "Missing arguments:" << argv[0] << " specify size of matrix N"
         << endl;
    exit(1);
  } else {
    n = atoi(argv[1]); // N from as ci
  }

  double h = 1.0 / (n + 1); // step length
  double eps = 1.0e-8;      // tolerance

  // Creating tridiagonal matrix
  mat A = zeros<mat>(n, n);             // matrix for A
  A.diag() += double(2.0) / pow(h, 2);  // center diagonal
  A.diag(1) -= double(1.0) / pow(h, 2); // upper diagonal
  // A.diag(-1) -= double(1.0) / pow(h, 2); // lower diagonal

  // creating a while loop that checks whether the off diagonal elements are
  // larger than eps
  // while (maxelm > eps) {
  //}
  return 0;
}
