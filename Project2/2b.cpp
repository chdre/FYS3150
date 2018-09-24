#include <algorithm>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <typeinfo>
//-DARMA_DONT_USE_WRAPPER -lblas -llapack

using namespace std;
using namespace arma;

ofstream outfile;

double max_offdiag(mat &A, int n, int *l, int *k) {
  // function to find the maximum value of the array A
  double maxelm = 0.0;

  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      if (abs(A(i, j)) > maxelm) {
        maxelm = abs(A(i, j));
        *l = i; // saving indices of max element i,j
        *k = j;
      }
    }
  }
  return maxelm;
}

void rotate(mat &A, mat &R, int &k, int &l, int n) {
  // rotating the matrix A, saving eigenvalues in vector R
  double c, s, t, tau;

  if (A(k, l) != 0.0) { // making sure we do not divide by 0
    tau = (A(l, l) - A(k, k)) / (2 * A(k, l)); // tau = cot2\theta
    if (tau >= 0) {
      t = 1.0 / (tau + sqrt(1 + pow(tau, 2)));
    } else {
      t = -1.0 / (-tau + sqrt(1 + pow(tau, 2))); // tau + sqrt(1 + pow(tau, 2));
    }
    c = 1.0 / sqrt(1 + pow(t, 2));
    s = t * c;
  } else {
    c = 1.0;
    s = 0.0;
  }
  double A_ll = A(l, l);
  double A_kk = A(k, k);
  double A_kl = A(k, l);

  // changing matrix elements
  A(k, k) = pow(c, 2) * A_kk - 2.0 * c * s * A_kl + pow(s, 2) * A_ll;
  A(l, l) = pow(s, 2) * A_kk + 2.0 * c * s * A_kl + pow(c, 2) * A_ll;
  A(k, l) = 0.0;
  A(l, k) = 0.0;

  // change remaining elements
  for (int i = 0; i < n; i++) {
    if (i != k && i != l) {
      double A_ik = A(i, k);
      double A_il = A(i, l);

      A(i, k) = c * A_ik - s * A_il;
      A(k, i) = A_ik;
      A(i, l) = c * A_il + s * A_ik;
      A(l, i) = A_il;
    }
    // compute eigenvectors
    double R_il = R(i, l);
    double R_ik = R(i, k);

    R(i, k) = c * R_ik - s * R_il;
    R(i, l) = c * R_il + s * R_ik;
  }
  return;
}

void jacobi(mat &A, mat &R, int n, double h) {
  int k, l;            // indices for largest off diagonal element
  double eps = 1.0e-8; // tolerance

  R = zeros<mat>(n, n); // eigenvector matrix
  R.diag() += double(1.0);

  double max_offdiagval = max_offdiag(A, n, &l, &k); // max offdiag element
  double max_iter = pow(double(n), 3);               // max number of iterations
  int iter = 0;                                      // counter for iterations

  // creating a while loop that checks whether the off diagonal elements are
  // larger than eps. Calling the function max_offdiag to find.
  while (abs(max_offdiagval) > eps && iter < max_iter) {
    max_offdiagval = max_offdiag(A, n, &l, &k);
    rotate(A, R, k, l, n);
    iter++;
  }
  cout << "Number of iterations: " << iter << "\n";
  return;
}

main(int argc, char *argv[]) {
  int n;           // size of matrix
  string filename; // filename for output file

  if (argc <= 1) {
    cout << "Missing arguments:" << argv[0]
         << " specify filename for output and size n of matrices" << endl;
    exit(1);
  } else {
    filename = argv[1];
    n = atoi(argv[2]); // N from as ci
  }

  double h = 1.0 / (n - 1); // step length, preserving u(L) = 1

  // Creating tridiagonal matrix
  mat A = zeros<mat>(n, n);              // matrix for A
  A.diag() += double(2.0) / pow(h, 2);   // center diagonal
  A.diag(1) -= double(1.0) / pow(h, 2);  // upper diagonal
  A.diag(-1) -= double(1.0) / pow(h, 2); // lower diagonal

  mat R;

  vec eigval;
  mat eigvec;

  eig_sym(eigval, eigvec, A);

  jacobi(A, R, n, h);

  // writing to file
  outfile.open(filename);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (i == j) {
        outfile << A(i, j);
        outfile << " " << eigval[i] << endl;
      }
    }
  }
  outfile.close();

  return 0;
}
