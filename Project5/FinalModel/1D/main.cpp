#include <armadillo>
#include <cmath>
#include <iostream>
#include "ForwEuler.hpp"
#include "BackEuler.hpp"
#include "CrankNic.hpp"
#include "tridiag.hpp"
#include "analytical.hpp"

using namespace std;
using namespace arma;

int main() {
        double L = 1;
        int n = 120;
        double dx = L/n;
        double dt = 0.001;//pow(dx,2)/2.0;

        cout << dx << " " << dt << endl;

        int tmax = 1/dt;

        double alpha = dt/pow(dx,2);

        //FESolver(n, alpha, tmax);
        BESolver(n, alpha, tmax);
        //CNSolver(n, alpha, tmax);

        analytical1D(n, L/(n+1), tmax, dt, L);

}
