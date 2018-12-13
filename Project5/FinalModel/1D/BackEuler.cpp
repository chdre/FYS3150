#include "BackEuler.hpp"

// Backward Euler scheme for using a tridiagonal solver.

void BESolver(int n, double alpha, int tmax){
        vec u = zeros<vec>(n+2);  // Au = r
        vec r = zeros<vec>(n+2);  // vector for current time step

        // Boundary conditions
        u(0) = 8.0/1300;
        u(n+1) = 1.0;

        // Initial condition, setting the relative temperature from 8/1300->1 C
        // from top -> bottom
        double dT = (1.0 - 8.0/1300)/(n+1); // Temperature step
        for(int i = 0; i < n+2; i++) {
                u(i) = (i+1)*dT;
        }

        // Matrix elements of tridiagonal matrix
        double e = -alpha;
        double d = 1.0 + 2.0*alpha;

        ofstream outfile;
        outfile.open("BWEuler.txt");

        // writing initial state to file
        for(int i = 0; i < n+2; i++) outfile << u(i) << " ";
        outfile << endl;

        for(int j = 1; j < tmax; j++) {
                tridiagonal(d, e, u, r, n);

                // Preserving boundary conditions
                u(0) = 8.0/1300;
                u(n+1) = 1.0;
                r(0) = 8.0/1300;
                r(n+1) = 1.0;

                r = u;  // Setting right hand side of equation to u, for all i

                // writing to file
                for(int i = 0; i < n+2; i++) {
                        outfile << u(i) << " ";
                }
                outfile << endl;
        }
        outfile.close();
}
