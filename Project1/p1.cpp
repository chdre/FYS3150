#include <iostream>
#include <cmath>
#include <string>
#include <iomanip>
#include <fstream>
#include <armadillo>
#include <array>

using namespace std;
ofstream outfile;
//Functions used.
inline double f(double x, double h) {
        return pow(h,2)*100.0*exp(-10*x);
}
inline double closed_form(double x) {
        return 1.0 - (1.0 - exp(-10.0))*x - exp(-10.0*x);
}

int main(int argc, char *argv[]){
        int n;  // Defining exponent as an integer variable
        string filename;  // String filename for outputfile
        if(argc <= 1) {
                cout << "Missing arguments:" << argv[0] << " specify output filename and value of n" << endl;
                exit(1);
        }
        else{
                filename = argv[1];
                n = atoi(argv[2]);
        }
        double h = 1.0/(n+1); // Step size

        double*x = new double [n+2]; // x_i's

        // Vectors forming the tridiagonal matrix. The upper and lower diagonal are n-1 long
        double*a = new double [n+1];
        double*b = new double [n+2];
        double*c = new double [n+1];

        // Vectors for calculations
        double*u = new double [n+2];
        double*bt = new double [n+2]; // b tilde, diagonal elements (NOT RHS of equation)
        double*ft = new double [n+2]; // f tilde (RHS of equation)
        double*fprime = new double [n+2];

        clock_t start, finish;
        start = clock();

        for(int i = 0; i < n+2; i++) {
                b[i] = 2.0;
        }
        for(int i = 0; i < n+1; i++) {
                a[i] = -1.0;
                c[i] = -1.0;
        }
        //loop to update
        for(int i = 0; i < n+2; i++) {
                x[i] = double(i)*h;
                fprime[i] = f(x[i],h);
        }

        // Boundary conditions
        u[0] = 0;
        u[n+1] = 0;
        ft[0] = fprime[0];
        ft[1] = fprime[1];
        bt[0] = b[0];
        bt[1] = b[0];

        for(int i = 2; i < n+2; i++) { //n+1 eller n???
                bt[i] = b[i] - a[i-1]*c[i-1]/bt[i-1];
                ft[i] = fprime[i] - ft[i-1]*a[i-1]/bt[i-1];
        }
        finish = clock();
        cout << 1.0*(finish - start)/CLOCKS_PER_SEC << endl;
        outfile.open(filename);
        //outfile << "  x:        approx:          exact:       relative error:" << endl;
        for(int i = n; i > 0; i--) {
                u[i] = (ft[i] - c[i]*u[i+1])/bt[i];
        }
        for(int i=0; i < n+2; i++) {
                outfile << x[i];
                outfile << " " << u[i];
                outfile << " " << closed_form(x[i]) << endl;
        }
        outfile.close();
        return 0;
};
