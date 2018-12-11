#include "Jacobi.hpp"

void JSolver(double d, int n, double alpha, mat &u, mat rhs){
        double Qd;
        int maxIter = 1000;
        double diff = 1.0;
        int iter = 0;
        double tol = 1E-10;
        mat u_old;

        double CpRho = 1000*2.5;

        while( iter < maxIter && diff > tol) {
                diff = 0;
                u_old = u;
                for(int i = 1; i < n+1; i++) {
                        if(i <= 20) {
                                Qd = d + 1.4E-6/CpRho;
                        }
                        else if(i <= 40) {
                                Qd = d + 0.35E-6/CpRho;
                        }
                        else if(i <= 120) {
                                Qd = d + 0.05E-6/CpRho;
                        }
                        for(int j = 1; j < n+1; j++) {
                                u(i,j) = 1.0/Qd*(alpha*(u_old(i+1,j) + u_old(i,j+1) +
                                                        u_old(i-1,j)+ u_old(i,j-1)) + rhs(i,j));
                                diff += abs(u(i,j) - u_old(i,j));
                        }
                }
                iter++;
        }
        cout << iter << endl;
}
