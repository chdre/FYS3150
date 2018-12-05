#include <armadillo>
#include <cmath>
#include <iostream>
#include "ForwEuler.hpp"
#include "BackEuler.hpp"
#include "CrankNic.hpp"
#include "tridiag.hpp"

using namespace std;
using namespace arma;

int main() {
        int n = 10;
        double dx = 1.0/10;
        double dt = pow(dx,2)/2.0;

        int tmax = 1/dt;
        cout << tmax << endl;

        double alpha = dt/pow(dx,2);
        cout << "Alpha = " << alpha << endl;

        FESolver(n, alpha, tmax);
        //BESolver(n, alpha, tmax);
        //CNSolver(n, alpha, tmax);
}
