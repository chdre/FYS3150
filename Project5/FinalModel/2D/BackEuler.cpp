#include "BackEuler.hpp"
#include "Jacobi.hpp"

void BESolver(int n, double alpha, int tmax){
        mat u = ones<mat>(n+2,n+2);  // Au = r
        //mat r = zeros<mat>(n+2,n+2);  // vector for current time step of Crank-Nicolson
        u *= 8.0/1300;  // temperature 8/1300
        // Boundary conditions (u(0) set by zeros)
        for(int i = 0; i < n+2; i++) {
                u(n+1,i) = 1.0;
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

        int counter = 1;

        for(int j = 1; j < tmax; j++) {
                JSolver(d, n, alpha, u, u);

                // Preserving boundary conditions
                for(int i = 0; i < n+2; i++) {
                        u(n+1,i) = 1.0;
                        //r(0,i) = 8.0/1300;
                        //r(n+1,i) = 1.0;
                }

                //r = u;  // Setting right hand side of equation to u, for all i

                // writing to file
                if(j == counter || j == n-1) {
                        for(int i = 0; i < n+2; i++) {
                                for(int j = 0; j < n+2; j++) {
                                        outfile << u(i,j) << " ";
                                }
                                outfile << endl;
                        }
                        counter += 1000;
                }
        }
        outfile.close();

}
