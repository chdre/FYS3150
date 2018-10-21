#include <iostream>
#include <cmath>
#include <armadillo>
#include "planet.h"
#include "solver.h"

using namespace std;

int main(){
        double M_e, M_sun, M_sunval, G, Time, timestep;
        int n;

        Time = 5.0;  // time [years]
        n = 1000;    // steps
        timestep = Time/n;

        // constants
        G = 4.0*pow(M_PI,2);    // gravitational constant [AU³/yr²]

        // masses
        M_sunval = 2e30;
        M_sun = 1.0;
        M_e = 6.0e24/M_sunval;

        planet Sun(M_sun, 0, 0, 0, 0, 0, 0);
        planet Earth(M_e, 1, 0, 0, 0, 2.0*M_PI, 0);

        ofstream outfile;
        outfile.open("data/test.txt");
        for(int i=0; i <= n; i++) {
                solver solve(G, timestep, Earth, Sun);
                solve.VelocityVerlet(G, timestep, Earth, Sun);
                Earth.angularMom = Earth.angularMomentum(Sun);

                outfile << Earth.angularMom << " ";
                outfile << timestep*i << endl;
        }
        outfile.close();
};
