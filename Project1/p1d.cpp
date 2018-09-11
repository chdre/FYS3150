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
        // The solution for f, including the step length h squared
        return pow(h,2)*100.0*exp(-10*x);
}
inline double closed_form(double x) {
        // Closed form solution
        return 1.0 - (1.0 - exp(-10.0))*x - exp(-10.0*x);
}

int main(int argc, char *argv[]){
        string filename; // String filename for outputfile
        int j_max;       // max exponent of 10
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
                
                // Dynamical arrays for calculations
                double*x = new double [n+2];     // x_i's
                double*k = new double [n+2];    // (f tilde)/(b tilde)
                double*u = new double [n+2];
                double*bt = new double [n+2];   // b tilde, diagonal elements (NOT RHS of equation)
                double*ft = new double [n+2];   // f tilde (RHS of equation)
                double*fprime = new double [n+2]; // the function h^2*f = f'

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
                bt[0] = 2.0; // btilde_0 not included in calculations
                bt[1] = 2.0; // btilde_1 = 2 by the general equation

                for(int i = 2; i < n+2; i++) {
                        // updating btilde and ftilde
                        double bt_temp = 1.0/bt[i-1];   // temporary value so to not calculate 1/bt_{i-1} twice
                        bt[i] = 2.0 - bt_temp;
                        ft[i] = fprime[i] + ft[i-1]*bt_temp;
                }
                for(int i = n; i > 0; i--) {
                        // updating u
                        u[i] = (ft[i] + u[i+1])/bt[i];
                }
                double eps_max = log10(abs((u[n/10] - closed_form(x[n/10]))/closed_form(x[n/10]))); // setting max value as the first value
                for(int i = n/10+1; i < (9*n/10)+1; i++) {
                        // as to avoid closed_form = 0 we change the limits for the calculation. Example 
                        // n = 10 => i 1:10 instead of 0:11
                        double eps = log10(abs((u[i] - closed_form(x[i]))/closed_form(x[i])));
                        if(eps > eps_max) {
                                // checking if next value > previous, if so define new max value eps_max
                                eps_max = eps;
                        }
                }
                eps_max_arr[j-1] = eps_max;     // setting eps_max to array for a specific h=1/(n+1)
        }
        
        // Writing to file
        outfile.open(filename);
        for(int i=0; i < j_max; i++) {
                outfile << eps_max_arr[i];
                outfile << " " << log10h[i] << endl;
        }
        outfile.close();
        return 0;
}
