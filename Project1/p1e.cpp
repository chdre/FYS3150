#include <iostream>
#include <cmath>
#include <armadillo>

using namespace std;
using namespace arma;

ofstream outfile;

inline double f(double x, double h) {
        // The solution for f, including the step length h squared
        return pow(h,2)*100.0*exp(-10*x);
}
inline double closed_form(double x) {
        // Closed form solution
        return 1.0 - (1.0 - exp(-10.0))*x - exp(-10.0*x);
}


int main(int argc, char*argv[]){
        int n;  // Defining number of points as an integer variable
        string filename;  // String filename for outputfile

        if(argc <= 1) {
                cout << "Missing arguments:" << argv[0] << " specify output filename and value of n" << endl;
                exit(1);
        }
        else{
                // reading filename and value of n from command line
                filename = argv[1];
                n = atoi(argv[2]);      // setting n from ascii to integer
        }
        double h = 1.0/(n+1);

        // Matrices and Vectors
        mat fmat = zeros<vec>(n+2);     // matrix of f
        mat A = zeros<mat>(n+2,n+2);    // matrix for A
        double*x = new double [n+2];    // dynamic array for x_i

        for(int i = 1; i < n+2; i++) {
                // Updating x_i's and storing values for h^2*f in fprime array
                x[i] = i*h;
                fmat[i] = f(x[i],h);
        }

        // Changing elements along the diagonal of A
        A.diag() += double(2.0);        // center diagonal
        A.diag(1) -= double(1.0);       // upper diagonal
        A.diag(-1) -= double(1.0);      // lower diagonal

        vec u = zeros<vec>(n+2);        // defining vector for solution to LU decomposition
        u = solve(A,fmat);              // LU decomposing
        
        // writing to file
        outfile.open(filename);
        for(int i=0; i < n+2; i++) {
                outfile << x[i];
                outfile << " " << u[i];
                outfile << " " << closed_form(x[i]) << endl;
        }
        outfile.close();

        return 0;
}
