#include "tridiag.h"
//
//   Taking the adresse of vector u, and the setting up a Toeplitz matrix with
//   diagonal of value d, and off diagonal values e (the same for upper and lower).
//   Solving the equation Au = rhs by row reduction.

void tridiagonal(double d, double e, vec &u, vec rhs, int n){
        // for solving Au = rhs. (remember: Au = f, A: b(diag) and e(offdiag))
        vec dvec = ones<vec>(n)*d;
        vec evec = ones<vec>(n)*e;

        vec dt = zeros<vec>(n); // \tilde{d}
        vec rhst = zeros<vec>(n);

        // forward substitution
        for(int i = 2; i < n+2; i++) {
                temp = e[i]/dt[i-1];
                dt[i] = d[i] - e[i]*temp;
                rhst[i] = rhs[i] - rhst[i-1]*temp;
        }
        // backward substitution
        for(int i = n; i > 0; i--) {
                u[i] = (rhst[i] - c[i]*u[i+1])/dt[i]
        }
}
