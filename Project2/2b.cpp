#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using namespace std;
using namespace arma;


void jacobi(double ** A, int n){
  double max_offdiagval = max_offdiag(A, n, h, &l, &k, )
  // creating a while loop that checks whether the off diagonal elements are
  // larger than eps. Calling the function max_offdiag to find.
  while (max_offdiag > eps) {

 }
}

void rotate(double ** A, double ** R, int k, int l, int n) {
  // rotating the matrix A, saving eigenvalues in vector R
double c, s, tau, t;

tau = (A[l,l] - A[k,k])/(2*A[k,l]); // tau = cot2\theta
if(A[k,l] != 0.0){  // making sure we do not divide by 0 (orthognal rotation)
  if(tau > 0){
    t = -tau + sqrt(1 + pow(tau,2));
  }
  else {
    t = tau + sqrt(1 +pow(tau,2));
  }
  c = 1.0/sqrt(1 + pow(tau,2));
}
else{ // rotation is orthogonal
  c = 1.0;
  s = 0.0;
}

// changing matrix elements
A[k,k] = pow(c,2)*A[k,k] - 2.0*c*s*A[k,l] + pow(s,2)*A[l,l];
A[l,l] = pow(s,2)*A[k,k] + 2.0*c*s*A[k,l] + pow(c,2)*A[l,l];
A[k,l] = 0.0;
A[l,k] = 0.0;
// change remaining elements
for(int i = 0; i < n; i++){
  if(i != k && i != l) {
    A[i,k] = c*A[i,k] - s*A[i,l];
    A[k,i] = A[i,k];
    A[i,l] = c*A[i,l] + s*A[i,k];
    A[l,i] = A[i,l];
  }
  // compute eingenvectors
  R[i,k] = c*R[i,k] - s*R[i,l];
  R[i,l] = c*R[i,l] + s*R[i,k];
}
return;
}

double max_offdiag(mat &A, int n, double h, int *l, int *k) {
  /* function to find the maximum value of the array A. The diagonal is set to
   0. The index of the maximum value is found, since the function returns a
   single number we must transform this to indices for the matrix.*/
  mat A_temp = abs(A);  // absolute value of all elements in array
  A_temp -= eye(size(A)) * (double(2.0) / pow(h, 2)); // setting diagonal = 0
  uword max_index = A_temp.index_max(); // finding index of max element
  *l = int(max_index) / n;           // find correct index for column
  *k = int(max_index) - (n * l);     // find correct index for row
  double maxelm = A_temp(k, l);         // max element
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
  int k, l;

  // Creating tridiagonal matrix
  mat * A = zeros<mat>(n, n);             // matrix for A
  A.diag() += double(2.0) / pow(h, 2);  // center diagonal
  A.diag(1) -= double(1.0) / pow(h, 2); // upper diagonal
  A.diag(-1) -= double(1.0) / pow(h, 2); // lower diagonal

  mat * R = eye(n,n);  // eigenvector matrix

  return 0;
}
