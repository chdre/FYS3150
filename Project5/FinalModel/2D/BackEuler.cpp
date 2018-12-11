#include "BackEuler.hpp"
#include "Jacobi.hpp"

void BESolver(int n, double alpha, int tmax){
        mat u = zeros<mat>(n+2,n+2);  // Au = r

        // Initial condition, setting the temperature from 8->1300 C
        double dT = (1.0 - 8.0/1300)/(n+1); // Temperature step
        for(int i = 0; i < n+2; i++) {
                for(int j = 0; j < n+2; j++) {
                        u(i,j) = (i+1)*dT;
                }
        }

        // Boundary conditions (u(0) set by zeros)
        for(int i = 0; i < n+2; i++) {
                u(n+1,i) = 1.0;       // bottom
                u(0,i) = 8.0/1300;    // top
                u(i,0) = (i+1)*dT;    // left side
                u(i,n+1) = (i+1)*dT;  // right side
        }

        // Matrix elements of tridiagonal matrix
        double e = -alpha;
        double d = 1.0 + 4.0*alpha;

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
                JSolver(d, n, alpha, u, u);

                // Preserving boundary conditions
                for(int i = 0; i < n+2; i++) {
                        u(n+1,i) = 1.0;
                }

                // writing to file
                //if(j == counter || j == n-1) {
                for(int i = 0; i < n+2; i++) {
                        for(int j = 0; j < n+2; j++) {
                                outfile << u(i,j) << " ";
                        }
                        outfile << endl;
                }
                counter += 100;
                //  }
        }
        outfile.close();

}
