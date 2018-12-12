#include "BackEuler.hpp"
#include "Jacobi.hpp"

void BESolver(int n, double alpha, int tmax, double dx, double dt){
        mat u = zeros<mat>(n+2,n+2);  // Au = r

        // Boundary conditions (u(0) set by zeros)
        for(int i = 0; i < n+2; i++) {
                u(n+1,i) = 1.0;
        }

        // Matrix elements of tridiagonal matrix
        double e = -alpha;
        double d = 1.0 + 4.0*alpha;

        mat A = zeros<mat>(n+2,n+2);
        A.diag() += d; // center diagonal
        A.diag(1) += e; // upper diagonal
        A.diag(-1) += e; // lower diagonal

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
                JSolver(e, d, n, alpha, u, u);

                // Preserving boundary conditions
                for(int i = 0; i < n+2; i++) {
                        u(0,i) = 0;
                        u(n+1,i) = 1.0;
                }

                r = u;  // Setting right hand side of equation to u, for all i

                // writing to file
                if(j == counter || j == tmax-1) {
                        for(int i = 0; i < n+2; i++) {
                                for(int j = 0; j < n+2; j++) {
                                        outfile << u(i,j) << " ";
                                }
                                outfile << endl;
                        }
                        counter += 100;
                }
        }
        outfile.close();

}
