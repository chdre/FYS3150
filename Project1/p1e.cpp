#include <iostream>
#include <cmath>
#include <armadillo>

using namespace std;
using namespace arma;

inline double f(double x, double h) {
        return pow(h,2)*100.0*exp(-10*x);
}


int main(int argc, char*argv[]) {
        int n; // number of points

        if(argc <= 1) {
                cout << "Missing arguments:" << argv[0] << " value of n" << endl;
                exit(1);
        }
        else{
                n = atoi(argv[1]);
        }
        double h = 1.0/(n+1);
        mat fmat = zeros<mat>(n,1); // matrix of f

        for(int i = 0; i < n+1; i++) {
                double x = i*h;
                fmat[i] = f(x,h);
        }

        mat A = zeros<mat>(n,n);
        A.diag() += 2;
        A.diag(1) -= 1;
        A.diag(-1) -= 1;



        return 0;

}
