#include <iostream>
#include <cmath>
#include <string>
#include <iomanip>
#include <fstream>
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
        string filename; // String filename for outputfile
        int j_max;
        if(argc <= 1) {
                cout << "Missing arguments:" << argv[0] << " specify output filename and max exponent of n" << endl;
                exit(1);
        }
        else{
                filename = argv[1];
                j_max = atoi(argv[2]);
        }
        double*eps_max_arr = new double [j_max];
        double*log10h = new double [j_max];

        // Loop to calculate maximum value of relative error.
        for(int j = 1; j < j_max+1; j++) {
                int n = pow(10,j);
                double h = 1.0/(n+1); // Step size
                log10h[j-1] = log10(h);
                // Vectors for calculations
                double*x = new double [n+2];     // x_i's
                double*k = new double [n+2]; // (f tilde)/(b tilde)
                double*u = new double [n+2];
                double*bt = new double [n+2]; // b tilde, diagonal elements (NOT RHS of equation)
                double*ft = new double [n+2]; // f tilde (RHS of equation)
                double*fprime = new double [n+2]; // our function f multiplied by h²

                //loop to update
                for(int i = 0; i < n+2; i++) {
                        x[i] = i*h;
                        fprime[i] = f(x[i],h);
                }

                // Boundary conditions
                u[0] = 0;
                u[n+1] = 0;
                ft[0] = fprime[0];
                ft[1] = fprime[1];
                bt[0] = 2.0; // Diagonal values
                bt[1] = 2.0; // Diagonal values

                for(int i = 2; i < n+2; i++) {
                        double bt_temp = 1.0/bt[i-1];
                        bt[i] = 2.0 - bt_temp;
                        ft[i] = fprime[i] + ft[i-1]*bt_temp;
                }
                for(int i = n; i > 0; i--) {
                        u[i] = (ft[i] + u[i+1])/bt[i];
                }
                double eps_max = log10(abs((u[n/10] - closed_form(x[n/10]))/closed_form(x[n/10])));
                for(int i = n/10+1; i < (9*n/10)+1; i++) {
                        double eps = log10(abs((u[i] - closed_form(x[i]))/closed_form(x[i])));
                        if(eps > eps_max) {
                                eps_max = eps;
                        }
                }
                eps_max_arr[j-1] = eps_max;
        }
        // Writing to file
        outfile.open(filename);
        //outfile << "  x:        approx:          exact:       relative error:" << endl;
        for(int i=0; i < j_max; i++) {
                outfile << eps_max_arr[i];
                outfile << " " << log10h[i] << endl;
        }
        outfile.close();
        return 0;
}
