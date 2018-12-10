#include "tridiag.hpp"

//   Taking the adresse of vector u, and the setting up a Toeplitz matrix with
//   diagonal of value d, and off diagonal values e (the same for upper and lower).
//   Solving the equation Au = rhs by row reduction.

void tridiagonal(double d, double e, vec &u, vec rhs, int n){
        // for solving Au = rhs. (remember: Au = f, A: b(diag) and e(offdiag))
        vec dt = zeros<vec>(n+2); // \tilde{d}
        vec rhst = zeros<vec>(n+2);

        dt(0) = d;

        // forward substitution
        for(int i = 1; i < n+2; i++) {
                double temp = e/dt(i-1);
                dt(i) = d - e*temp;
                rhst(i) = rhs(i) - rhst(i-1)*temp;
        }
        // backward substitution
        for(int i = n; i > 0; i--) {
                u(i) = (rhst(i) - e*u(i+1))/dt(i);
        }
}
