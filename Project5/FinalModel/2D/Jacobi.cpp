#include "Jacobi.hpp"

void JSolver(double d, int n, double alpha, mat &u, mat rhs){
        int maxIter = 1000;
        double diff = 1.0;
        int iter = 0;
        double tol = 1E-10;
        mat u_old;

        while( iter < maxIter && diff > tol) {
                diff = 0;
                u_old = u;
                for(int i = 1; i < n+1; i++) {
                        for(int j = 1; j < n+1; j++) {
                                if(i <= 20) {
                                        d += 1.4E-6/2.5;
                                }
                                else if(i <= 40) {
                                        d += 0.35E-6/2.5;
                                }
                                else if(i <= 120) {
                                        d += 0.05E-6/2.5;
                                }
                                u(i,j) = 1.0/d*(alpha*(u_old(i+1,j) + u_old(i,j+1) +
                                                       u_old(i-1,j)+ u_old(i,j-1)) + rhs(i,j));
                                diff += abs(u(i,j) - u_old(i,j));
                        }
                }
                iter++;
        }
        cout << iter << endl;
}
