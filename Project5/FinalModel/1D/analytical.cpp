#include "analytical.hpp"

void analytical1D(int n, double dx, double tmax, double L){
        double u, x, t, piL, mpiL, mpi;

        piL = M_PI/L;

        ofstream outfile;
        outfile.open("Analytical1D.txt");
        for(double l = 0; l < tmax; l++) {
                t = l/tmax;
                outfile << 0 << " ";  // BC
                for(int i = 1; i < n+1; i++) {
                        x = i*dx;
                        u = 0;
                        for(int m = 1; m < 1000; m++) {
                                mpiL = m*piL;
                                mpi = m*M_PI;
                                u += 2*L*cos(mpi)/(mpi)*sin(mpiL*x)*exp(-pow(mpi/L,2)*t);
                        }
                        outfile << u + x << " ";
                }
                outfile << 1 << endl; // BC
        }
        outfile.close();
}
