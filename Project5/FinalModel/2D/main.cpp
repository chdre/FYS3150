#include <armadillo>
#include <cmath>
#include <iostream>
#include "ForwEuler.hpp"
#include "BackEuler.hpp"
//#include "CrankNic.hpp"
//#include "tridiag.hpp"
#include "analytical.hpp"

using namespace std;
using namespace arma;

int main() {
        int n = 118;
        double dx = 1.0/(n+2);
        double dt = 0.0001;//pow(dx,2)/4.0;

        int tmax = 1/dt;

        double alpha = dt/pow(dx,2);

        //FESolver(n, alpha, tmax);
        BESolver(n, alpha, tmax);
        //CNSolver(n, alpha, tmax);

        //analytical2D(n, dx, 1, dt);

        mat A = zeros<mat>(10,10);
        for(int i = 0; i < 10; i++) {
                A(9,i) = 1;
        }

        cout << A << endl;

}
