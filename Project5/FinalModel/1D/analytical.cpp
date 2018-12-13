#include "analytical.hpp"

// analytical solution of 1D. We have set L = 1.

void analytical1D(int n, double dx, double tmax, double dt){
        double u, x, t, mpi;

        ofstream outfile;
        outfile.open("Analytical1D.txt");
        for(double l = 0; l < tmax; l++) {
                t = l*dt;
                outfile << 0 << " ";  // BC
                for(int i = 1; i < n+1; i++) {
                        x = i*dx;
                        u = 0;
                        for(int m = 1; m < 1000; m++) {
                                mpi = m*M_PI;
                                u += 2*L*cos(mpi)/(mpi)*sin(mpi*x)*exp(-pow(mpi,2)*t);
                        }
                        outfile << u + x << " ";
                }
                outfile << 1 << endl; // BC
        }
        outfile.close();
}
