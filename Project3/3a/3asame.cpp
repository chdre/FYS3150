#include <armadillo>
#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;
using namespace arma;
ofstream outfile;

main(int argc, char *argv[]){
        int n;      // number of points
        double T;   // simulation time [years]
        string filename;
        int euler;
        if (argc <= 1) {
                cout << "Provide number of steps n and simulation time T" << endl;
                exit(1);
        }
        else{
                n = atoi(argv[1]);
                T = at::of(argv[2]);
                euler = atoi(argv[3]);
                filename = argv[4];

        }
        double AU = 1.0;//1.5e11;           // Astronomical unit [m]
        double M_s = 1.0;//;          // Mass of sun = 1
        double GM_s = 4.0*pow(M_PI,2);
        double M_e = 6.10e24/2.10e30; // Relative mass of the earth [1/M_sun]

        double dt = T/n;


        double vx = 0;
        double vy = 2.0*M_PI;
        double x = 1.0;
        double y = 0;

        if(euler == 1) {
                outfile.open(filename);
                double ax, ay, r;
                for(int i = 0; i <= n; i++) {
                        // calling function and calculating distance earth-sun
                        r = sqrt(pow(x,2) + pow(y,2));
                        ax = -GM_s*x/pow(r,3);
                        ay = -GM_s*y/pow(r,3);

                        vx += ax*dt;
                        vy += ay*dt;

                        x += vx*dt;
                        y += vy*dt;

                        outfile << x;
                        outfile << " " << y << endl;;
                }
        }
        if(euler == 0) {
                outfile.open(filename);
                double ax, ay, r;
                r = sqrt(pow(x,2) + pow(y,2)); //r[i]
                for(int i = 0; i <= n; i++) {
                        ax = -GM_s*x/pow(r,3);    // a[i]
                        ay = -GM_s*y/pow(r,3);

                        x += dt*vx + pow(dt,2)/2.0*ax;  //x[i+1]
                        y += dt*vy + pow(dt,2)/2.0*ay;

                        r = sqrt(pow(x,2) + pow(y,2));  //r[i+1]

                        vx += dt/2.0*(-GM_s*x/pow(r,3) + ax);  // v[i+1]
                        vy += dt/2.0*(-GM_s*y/pow(r,3) + ay);


                        //vx += dt/2.0*(ax + )
                        outfile << x;
                        outfile << " " << y << endl;;
                }
        }
}
