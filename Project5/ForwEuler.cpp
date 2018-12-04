#include "ForwEuler.hpp"

void FWSolver(int n, double alpha, int tmax){
        vec u = zeros<vec>(n+1);

        // Boundary condition (u(0) set by zeros)
        u(n) = 1.0;

        for(int j = 1; j < tmax; j++) {
                for(int i = 1; i < n; i++) {
                        u(i) = (1.0 - 2.0*alpha)*u(i) + alpha*u(i+1) + alpha*u(i-1);    // RHS of equation
                }
                // writing to file
                ofstream outfile;
                outfile.open("FWEuler.txt");
                for(int i = 0; i < n+1; i++) {
                        outfile << u(i) << endl;
                }
                outfile.close();
                cout << "Closed file" << endl;
        }
}
