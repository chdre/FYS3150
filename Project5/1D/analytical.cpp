#include "analytical.hpp"

void analytical1D(int n, double dx, double tmax){
        double u, x, t;

        ofstream outfile;
        outfile.open("Analytical1D.txt");
        for(double l = 0; l < tmax; l++) {
                t = (double) l/tmax;
                outfile << 0 << " ";
                for(double i = 1; i < n+1; i++) {
                        x = i*dx;
                        u = 0;
                        for(int m = 1; m < 10; m++) {
                                u += 2*cos(M_PI*m)/(M_PI*m)*sin(m*M_PI*x)*exp(-pow(m*M_PI,2)*t);
                        }
                        outfile << u + x << " ";
                }
                outfile << 1 << endl;
        }
        outfile.close();
}
