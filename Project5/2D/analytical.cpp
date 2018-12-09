#include "analytical.hpp"

void analytical2D(int n, double dx, double t){
        mat u = zeros<mat>(n+1,n+1);
        double x, y;

        ofstream outfile;
        outfile.open("Analytical2D.txt");
        for(double i = 0; i < n+1; i++) {
                x = i*dx;
                cout << x << endl;
                for(double j = 0; j < n+1; j++) {
                        y = j*dx;
                        u(i,j) = sin(M_PI*x)*sin(M_PI*y)*exp(-2*pow(M_PI,2)*t);
                }
        }
        // writing to file
        for(int i = 0; i < n+1; i++) {
                for(int j = 0; j < n+1; j++) {
                        outfile << u(i,j) << " ";
                }
                outfile << endl;
        }
        outfile.close();
}
