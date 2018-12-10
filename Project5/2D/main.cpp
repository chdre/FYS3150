#include <armadillo>
#include <cmath>
#include <iostream>
#include "ForwEuler.hpp"
#include "BackEuler.hpp"
//#include "CrankNic.hpp"
//#include "tridiag.hpp"
#include "analytical.hpp"
#include "Jacobi.hpp"

using namespace std;
using namespace arma;

int main() {
        int n = 100;
        double dx = 1.0/120;
        double dt = pow(dx,2)/4.0;

        int tmax = 1/dt;
        cout << tmax << endl;

        double alpha = dt/pow(dx,2);
        cout << "Alpha = " << alpha << endl;

        //FESolver(n, alpha, tmax);
        BESolver(n, alpha, tmax, dx, dt);

        //analytical2D(n, dx, 1);
}
