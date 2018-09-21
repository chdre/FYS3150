#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using namespace std;
using namespace arma;

main(int argc, char *argv[]) {
  int N; // size of matrix

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
  while (max(pow(a[i, j]), 2) > eps) {
  }

  double; //?
  for (int i = 1; i < n; i++) {
    for (int j = 1; j < n; j++) {
      if (pow(A[i, j], 2) > eps) {
      }
    }
  }
  return 0;
}
