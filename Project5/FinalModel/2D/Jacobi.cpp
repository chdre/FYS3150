#include "Jacobi.hpp"

/* Jacobi iterative method for Backward Euler scheme. The solver takes the adresse
   of matrix u, which it changes over each iteration. The solver also takes rhs,
   the right hand side of the equation Au(t+dt)=u(t), which is unchanged during
   the iterations.
 */
void JSolver(double d, int n, double alpha, mat &u, mat rhs){
        double Qd;    // new potential, depends on debth
        mat u_old;    // u(t), used to calculate u(t+dt)

        int maxIter = 1000;   // maximum number of iterations
        int iter = 0;         // counter for iterations

        double diff = 1.0;
        double tol = 1E-10;

        double MYr = 3600*24*365*1E6;

        double CpRho = 1000*3.510;  // specific heat capacity times the density

        while( iter < maxIter && diff > tol) {
                diff = 0;
                u_old = u;
                for(int i = 1; i < n+1; i++) {
                        if(i <= 20) {
                                Qd = d + 1.4E-6*MYr/CpRho;
                        }
                        else if(i <= 40) {
                                Qd = d + 0.35E-6*MYr/CpRho;
                        }
                        else if(i <= 120) {
                                Qd = d + 0.05E-6*MYr/CpRho;
                        }
                        for(int j = 1; j < n+1; j++) {
                                // Backward Euler
                                u(i,j) = 1.0/Qd*(alpha*(u_old(i+1,j) + u_old(i,j+1) +
                                                        u_old(i-1,j)+ u_old(i,j-1)) + rhs(i,j));
                                diff += abs(u(i,j) - u_old(i,j));
                        }
                }
                iter++;
        }
        cout << iter << endl;
}
