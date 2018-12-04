#include <armadillo>
#include <cmath>
#include <iostream>
#include "ForwEuler.hpp"
#include "BackEuler.hpp"
#include "CrankNic.hpp"
#include "tridiag.hpp"

using namespace std;
using namespace arma;

void FWSolver(int n, double alpha, int tmax);


int main(int argc, char *argv[]) {
        int n = 1000;
        double dx = atof(argv[0]);
        double dt = 2.0*dx;

        int tmax = (int) dt*n;

        double alpha = pow(dx/dt,2);

        FWSolver(n, alpha, tmax);
        cout << "end of main" << endl;
}
