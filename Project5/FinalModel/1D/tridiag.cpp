#include "tridiag.hpp"

//   Taking the adresse of vector u, and the setting up a Toeplitz matrix with
//   diagonal of value d, and off diagonal values e (the same for upper and lower).
//   Solving the equation Au = rhs by row reduction.

void tridiagonal(double d, double e, vec &u, vec rhs, int n){
        double Qd; // new potential, depends on debth

        // for solving Au = rhs.
        vec dt = zeros<vec>(n+2); // \tilde{d}
        vec rhst = zeros<vec>(n+2);

        dt(0) = d;

        double MYr = 3600*24*365*1E6;
        double CpRho = 1000*3.510;  // specific heat capacity times the density


        // forward substitution
        for(int i = 1; i < n+2; i++) {
                if(i <= 20) {
                        Qd = d + 1.4E-6*MYr/CpRho;
                }
                else if(i <= 40) {
                        Qd = d + 0.35E-6*MYr/CpRho;
                }
                else if(i <= 120) {
                        Qd = d + 0.05E-6*MYr/CpRho;
                }

                double temp = e/dt(i-1); // temporary value
                dt(i) = Qd - e*temp;
                rhst(i) = rhs(i) - rhst(i-1)*temp;
        }
        // backward substitution
        for(int i = n; i > 0; i--) {
                u(i) = (rhst(i) - e*u(i+1))/dt(i);
        }
}
