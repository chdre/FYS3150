#include <iostream>
#include <cmath>
#include <armadillo>
#include "planet.h"
#include "solver.h"

using namespace std;

int main(){
        double M_e, G;

        double Time = 2.0;  // years
        int n = 1000;               // steps

        double timestep = Time/n;

        G = 4.0*pow(M_PI,2);                   // gravitational constant
        M_e = 6.0e24/2e30;

        planet Sun(1, 0, 0, 0, 0, 0, 0);
        planet Earth(M_e, 1, 0, 0, 0, 2.0*M_PI, 0);

        ofstream outfile;
        outfile.open("test.txt");
        for(int i=0; i <= n; i++) {
                solver solve(G, timestep, Earth, Sun);
                solve.VelocityVerlet(G, timestep, Earth, Sun);
                outfile << Earth.position[0] << " ";
                outfile << Earth.position[1] << endl;
        }
        outfile.close();
};
