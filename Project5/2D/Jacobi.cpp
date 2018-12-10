#include "Jacobi.hpp"

void JSolver(double e, double d, int n, double alpha, mat &u){
        int maxIter = 1000;
        double diff = 1.0;
        int iter = 0;
        double tol = 0;
        mat u_old;

        while( iter < maxIter) {
                diff = 0;
                u_old = u;
                for(int i = 1; i < n+1; i++) {
                        for(int j = 1; j < n+1; j++) {
                                u(i,j) = 1.0/d*(alpha*(u_old(i+1,j) + u_old(i,j+1) +
                                                       u_old(i-1,j)+ u_old(i,j-1)) + u_old(i,j));
                                diff += abs(u(i,j) - u_old(i,j));
                        }
                }
                iter++;
        }
        cout << iter << endl;
}

/*
   void JSolver(double e, double d, int n, double alpha, mat &u){
        double sum, diff;

        double tol = 1E-10;

        int maxIter = 100000;
        int iter = 0;

        mat u_old = u;
        diff = 1;

        mat A = zeros<mat>(n+2,n+2);
        A.diag() += d; // center diagonal
        A.diag(1) += e; // upper diagonal
        A.diag(-1) += e; // lower diagonal

        while(iter < maxIter && diff > tol) {
                diff = 0;
                for(int i = 0; i < n+2; i++) {
                        sum = 0;
                        for(int j = 0; j < n+2; j++) {
                                if(i != j) {
                                        sum +=  A(i,j)*u_old(i,j);
                                }// end if
                        } // end j1
                        for(int j = 0; j < n+2; j++) {
                                u(i,j) = (u_old(i,j) - sum)/d;
                                diff += abs(u(i,j) - u_old(i,j));
                        }// end j2
                }// end i
                u_old = u;
                iter++;
        }
        cout << iter << endl;
   }
 */
/*
   void JSolver(int n, ){
        int MaxIter = 1000000;
        abstol =

        for(int k = 0; k < MaxIterations; k++) {
                for(int i = 1; i < n+2; i++) {
                        for(int j=1; j < n+2; j++) {
                                u(i,j) = dt*u_old(i,j) + u_old(i,j) +
                                         alpha*(u_old(i+1,j) + u_old(i,j+1) - 4.0*u_old(i,j) +
                                                u_old(i-1,j) + u_old(i,j-1));
                        }
                }
                double sum = 0.0;
                for(int i = 0; i < N; i++) {
                        for(int j = 0; j < N; j++) {
                                sum += (u_old(i,j) - u_(i,j))*(u_old(i,j)-u(i,j));
                                u_old(i,j) = u_(i,j);
                        }
                }
                if(sqrt (sum) < abstol) {
                        return k;
                }
        }
   }*/
/*
   void JSolver(mat A, mat &u, double dx, double dt, int n, double alpha) {
        mat u_temp;

        double max_itr = pow(n,3);  // maximum iteration
        double itr = 0;                    // iteration counter
        double tol = 1e-10;
        double difference = 1;

        while((itr <= max_itr) && (difference > tol)) {
                u_temp = u;
                difference = 0;

                for(int i = 1; i < n+2; i++) {
                        for(int j = 1; j < n+2; j++) {
                                u(i,j) = dt*u(i,j) + u_temp(i,j) +
                                         alpha*(u_temp(i+1,j) + u_temp(i,j+1) -
                                                4*u_temp(i,j) + u_temp(i-1,j)+ u_temp(i,j-1));

                                difference += fabs(u(i,j) - u_temp(i,j));
                        }
                }
                itr++;
        }
        cout << itr << endl;
   }
 */


















//
