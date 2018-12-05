#include "ForwEuler.hpp"

void FESolver(int n, double alpha, int tmax){
        mat u = zeros<mat>(n+2,n+2);

        // Boundary condition (u(0) set by zeros)
        for(int i = 0; i < n+2; i++) {
                u(n+1,i) = 1.0;
        }

        ofstream outfile;
        outfile.open("FWEuler.txt");

        for(int l = 1; l < tmax; l++) {
                for(int i = 1; i < n+1; i++) {
                        for(int j = 1; j < n+1; j++) {
                                double delta = u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1);
                                u(i,j) = (1.0 - 4.0*alpha)*u(i,j) + alpha*delta;
                        }
                }

                // Preserving boundary conditions
                for(int i = 0; i < n+2; i++) {
                        u(0,i) = 0;
                        u(n+1,i) = 1.0;
                }

                // writing to file
                for(int i = 0; i < n+2; i++) {
                        for(int j = 0; j < n+2; j++) {
                                outfile << u(i,j) << " ";
                        }
                        outfile << endl;
                }

        }
        outfile.close();
}
