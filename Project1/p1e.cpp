#include <iostream>
#include <cmath>
#include <armadillo>

using namespace std;
using namespace arma;

ofstream outfile;

inline double f(double x, double h) {
        return pow(h,2)*100.0*exp(-10*x);
}
inline double closed_form(double x) {
        return 1.0 - (1.0 - exp(-10.0))*x - exp(-10.0*x);
}


int main(int argc, char*argv[]){
        int n; // number of points
        string filename;  // String filename for outputfile

        if(argc <= 1) {
                cout << "Missing arguments:" << argv[0] << " specify output filename and value of n" << endl;
                exit(1);
        }
        else{
                filename = argv[1];
                n = atoi(argv[2]);
        }
        double h = 1.0/(n+1);

        // Matrices and Vectors
        mat fmat = zeros<vec>(n+2);            // matrix of f
        mat A = zeros<mat>(n+2,n+2);
        double*x = new double [n+2];

        for(int i = 1; i < n+2; i++) {
                x[i] = i*h;
                fmat[i] = f(x[i],h);
        }

        // Changing diagonal elements of A
        A.diag() += double(2.0);
        A.diag(1) -= double(1.0);
        A.diag(-1) -= double(1.0);

        vec u = zeros<vec>(n+2);
        u = solve(A,fmat);

        outfile.open(filename);
        //outfile << "  x:        approx:          exact:       relative error:" << endl;
        for(int i=0; i < n+2; i++) {
                outfile << x[i];
                outfile << " " << u[i];
                outfile << " " << closed_form(x[i]) << endl;
        }
        outfile.close();

        return 0;
}
