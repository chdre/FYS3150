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
        // The solution for f, including the step length f squared
        return pow(h,2)*100.0*exp(-10*x);
}
inline double closed_form(double x) {
        // Closed form solution
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
                // reading filename and value of n from command line
                filename = argv[1];
                n = atoi(argv[2]);      // setting n from ascii to integer
        }
        double h = 1.0/(n+1); // Step size

        double*x = new double [n+2]; // x_i's

        // Dynamic arrays forming the tridiagonal matrix. The upper and lower diagonal are n+1 long
        double*a = new double [n+1];
        double*b = new double [n+2];
        double*c = new double [n+1];

        // Dynamical arrays for calculations
        double*u = new double [n+2];
        double*bt = new double [n+2]; // b tilde, diagonal elements (NOT RHS of equation)
        double*ft = new double [n+2]; // f tilde (RHS of equation)
        double*fprime = new double [n+2];       // the function h^2*f = f'

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
                // Updating x_i's and storing values for h^2*f in fprime array
                x[i] = i*h;
                fprime[i] = f(x[i],h);
        }

        // Boundary conditions
        u[0] = 0;
        u[n+1] = 0;
        ft[0] = fprime[0];      // ftilde_0 not included in calculations
        ft[1] = fprime[1];      // ftilde_1 = h^2*f by the general equation
        bt[0] = b[0];   // btilde_0 not included in calculations
        bt[1] = b[0];   // btilde_1 = -2 by the general equation

        for(int i = 2; i < n+2; i++) {
                // updating btilde and ftilde
                bt[i] = b[i] - a[i-1]*c[i-1]/bt[i-1];
                ft[i] = fprime[i] - ft[i-1]*a[i-1]/bt[i-1];
        }
        for(int i = n; i > 0; i--) {
                // updating u
                u[i] = (ft[i] - c[i]*u[i+1])/bt[i];
        }
        
        finish = clock();
        cout << 1.0*(finish - start)/CLOCKS_PER_SEC << endl;
        
        // printing to file
        outfile.open(filename);
        for(int i=0; i < n+2; i++) {
                outfile << x[i];
                outfile << " " << u[i];
                outfile << " " << closed_form(x[i]) << endl;
        }
        outfile.close();
        return 0;
};
