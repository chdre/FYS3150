#include <iostream>
#include <cmath>
#include <string>
#include <iomanip>
#include <fstream>
#include <armadillo>
#include <array>
#include <time.h>
#include <ctime>

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

        // Vectors for calculations
        double*x = new double [n+2];      // x_i's
        double*k = new double [n+2];      // (f tilde)/(b tilde)
        double*u = new double [n+2];
        double*bt = new double [n+2];     // b tilde, diagonal elements (NOT RHS of equation)
        double*ft = new double [n+2];     // f tilde (RHS of equation)
        double*fprime = new double [n+2]; // our function f multiplied by h²
        double*eps = new double [n];      // relative error

        //loop to update
        for(int i = 0; i < n+2; i++) {
                x[i] = double(i)*h;
        }
        for(int i = 0; i < n+2; i++) {
                fprime[i] = f(x[i],h);
        }

        // Boundary conditions
        u[0] = 0;
        u[n+1] = 0;
        ft[0] = fprime[0];
        ft[1] = fprime[1];
        bt[0] = 2.0;  // Diagonal values
        bt[1] = 2.0;  // Diagonal values

        for(int i = 2; i < n+2; i++) {
                double bt_temp = 1.0/bt[i-1];
                bt[i] = 2.0 - bt_temp;
                ft[i] = fprime[i] + ft[i-1]*bt_temp;
        }
        for(int i = n; i > 0; i--) {
                u[i] = (fprime[i] + u[i+1])/bt[i];
        }
        for(int i = 1; i < n+1; i++) {
                eps[i] = log10(abs((u[i] - ft[i])/ft[i]))
        }
        // Writing to file
        outfile.open(filename);
        //outfile << "  x:        approx:          exact:       relative error:" << endl;
        for(int i=0; i < n+2; i++) {
                outfile << x[i];
                outfile << " " << u[i];
                outfile << " " << closed_form(x[i]);
                outfile << " " << eps[i] << endl;
        }
        outfile.close();
        return 0;
};
