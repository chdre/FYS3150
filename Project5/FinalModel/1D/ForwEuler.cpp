#include "ForwEuler.hpp"

void FESolver(int n, double alpha, int tmax){
        vec u = zeros<vec>(n+2);

        // Boundary condition (u(0) set by zeros)
        u(n+1) = 1.0;

        ofstream outfile;
        outfile.open("FWEuler.txt");

        // writing initial state to file
        for(int i = 0; i < n+2; i++) outfile << u(i) << " ";
        outfile << endl;

        for(int j = 1; j < tmax; j++) {
                for(int i = 1; i < n+1; i++) {
                        u(i) = (1.0 - 2.0*alpha)*u(i) + alpha*u(i+1) + alpha*u(i-1);    // RHS of equation
                }

                // Preserving boundary conditions
                u(0) = 0;
                u(n+1) = 1.0;

                // writing to file
                for(int i = 0; i < n+2; i++) {
                        outfile << u(i) << " ";
                }
                outfile << endl;
        }
        outfile.close();
}
