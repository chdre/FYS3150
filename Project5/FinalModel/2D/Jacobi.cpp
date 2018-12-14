#include "Jacobi.hpp"

/* Jacobi iterative method for Backward Euler scheme. The solver takes the adresse
   of matrix u, which it changes over each iteration. The solver also takes rhs,
   the right hand side of the equation Au(t+dt)=u(t), which is unchanged during
   the iterations.
 */
void JSolver(int n, double dt, double alpha, mat &u, mat rhs){
        double Q;    // new potential, depends on debth
        mat u_old;    // u(t), used to calculate u(t+dt)

        bool radactive = true; // set false if radioactive enrichment not true

        int maxIter = 1000;   // maximum number of iterations
        int iter = 0;         // counter for iterations

        double diff = 1.0;
        double tol = 1E-10;

        double GYr = 3600*24*365*1E9; // GYr conversion factor
        double lscale = 120000;     // length of grid [m]
        double rho = 3.51E3;        // density []
        double Cp = 1000;           // specific heat capacity
        double k = 2.5;             // heat conductivity
        double tempsc = 1300;       // max temperature
        double scale1 = (k*GYr)/(Cp*rho*lscale*lscale);  // scaling length and time
        double scale2 = (dt*GYr)/(Cp*rho*tempsc);     // scaling temperature to radiactive perturbation

        double factor1 = 1.0/(1+4*alpha*scale1);  // to save FLOPS
        double factor2 = alpha*scale1;  // to save FLOPS

        while( iter < maxIter && diff > tol) {
                diff = 0;
                u_old = u;
                for(int i = 1; i < n+1; i++) {
                        if(i <= 20) {
                                Q = 1.4E-6*scale2;
                        }
                        else if(i <= 40) {
                                Q = 0.35E-6*scale2;
                        }
                        else if(i <= 120) {
                                // checking if radioective enrichment has happenened
                                if(radactive == true) Q = 0.05E-6*scale2 + 0.50E-6*scale2;
                                else if(radactive == false) Q = 0.05E-6*scale2;
                        }
                        for(int j = 1; j < n+1; j++) {
                                // Backward Euler
                                u(i,j) = factor1*(factor2*(u_old(i+1,j) + u_old(i,j+1) +
                                                           u_old(i-1,j)+ u_old(i,j-1)) + Q + rhs(i,j));
                                diff += abs(u(i,j) - u_old(i,j));
                        }
                }
                iter++;
        }
        cout << iter << endl;
}
