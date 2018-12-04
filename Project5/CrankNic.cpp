#include "CrankNic.hpp"

void CNSolver(int n, double alpha, int tmax){
        vec u = zeros<vec>(n);  // Au = r
        vec r = zeros<vec>(n);  // vector for current time step of Crank-Nicolson

        // Boundary conditions (u(0) set by zeros)
        u(n) = 1.0;

        // Matrix elements of tridiagonal matrix
        double e = -alpha;
        double d = 2.0 + 2.0*alpha;

        for(int j = 1; j < tmax; j++) {
                for(int i = 0; i < n+1; i++) {
                        r(i) = alpha*u(i-1) + (2.0 - 2.0*alpha)*u(i) + alpha*u(i+1);
                }
                r(0) = 0.0;
                r(n) = 1.0;

                tridiagonal(d, e, u, r, n);

                // preserving boundary conditions
                u(0) = 0.0;
                u(n) = 1.0;

                // writing to file
                outfile.open("CrankNic.txt");
                for(int i = 0; i < n+1; i++) {
                        outfile << u[i] << endl;
                }
                outfile.close();
        }
}
