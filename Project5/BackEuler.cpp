#include "BackEuler.hpp"

void BESolver(int n, double alpha, int tmax){
        vec u = zeros<vec>(n);  // Au = r
        vec r = zeros<vec>(n);  // vector for current time step of Crank-Nicolson

        // Boundary conditions (u(0) set by zeros)
        u(n) = 1.0;

        // Matrix elements of tridiagonal matrix
        double e = -alpha;
        double d = 1.0 + 2.0*alpha;

        for(int j = 1; j < tmax; j++) {
                tridiagonal(d, e, u, r, n);

                // Preserving boundary conditions
                u(0) = 0.0;
                u(n) = 1.0;

                r = u;  // Setting right hand side of equation to u, for all i

                // writing to file
                outfile.open("FWEuler.txt");
                for(int i = 0; i < n+1; i++) {
                        outfile << u[i] << endl;
                }
                outfile.close();

        }
}
