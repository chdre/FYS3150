#include "BackEuler.hpp"
#include "Jacobi.hpp"

// Backward Euler scheme using the Jacobi method.

void BESolver(int n, double alpha, int tmax, double dt){
        mat u = zeros<mat>(n+2,n+2);  // Au = r

        // Initial condition, setting the relative temperature from 8/1300->1 C
        // from top -> bottom
        double dT = (1.0 - 8.0/1300)/(n+1); // Temperature step
        for(int i = 0; i < n+2; i++) {
                for(int j = 0; j < n+2; j++) {
                        u(i,j) = (i+1)*dT;
                }
        }

        ofstream outfile;
        outfile.open("BWEuler2D.txt");

        // writing initial state to file
        for(int i = 0; i < n+2; i++) {
                for(int j = 0; j < n+2; j++) {
                        outfile << u(i,j) << " ";
                }
                outfile << endl;
        }

        int counter = 1;  // counter for writing to file

        for(int j = 1; j < tmax; j++) {
                JSolver(n, dt, alpha, u, u);

                // Preserving boundary conditions
                for(int i = 0; i < n+2; i++) {
                        u(n+1,i) = 1.0;       // bottom
                        u(0,i) = 8.0/1300;    // top

                }
                // writing to file
                if(j == counter || j == tmax-1) {
                        for(int i = 0; i < n+2; i++) {
                                for(int j = 0; j < n+2; j++) {
                                        outfile << u(i,j) << " ";
                                }
                                outfile << endl;
                        }
                        counter += 100;   // stepsize of time when printing
                }
        }
        outfile.close();
}
