#include <iostream>
#include <cmath>
#include <armadillo>
#include "planet.h"
#include "solver.h"

using namespace std;

int main(){
        double M_e, M_sun, M_j, M_sunval, G, Time, timestep, vfac, v0_e;
        int n;

        //beta = atof(argv[1]);
        //vfac = atof(argv[2]);    // factor to change the velocity

        Time = 24.0;    // time [years]
        n = 10000;    // steps
        timestep = Time/n;

        // constants
        G = 4.0*pow(M_PI,2);    // gravitational constant [AU³/yr²]

        // masses
        M_sunval = 2e30;          // mass of sun
        M_sun = 1.0;              // relative mass of sun
        M_e = 6.0e24/M_sunval;    // relative mass of earth
        M_j = 1.9e27/M_sunval;    // relative mass of jupiter

        planet Sun(M_sun, -0.0049, 0, 0, 0, -0.01, 0);
        planet Earth(M_e, 1-0.0049, 0, 0, 0, 2.0*M_PI, 0);
        planet Jupiter(M_j, 5.2-0.0049, 0, 0, 0, 0.7*M_PI, 0);

        /*ofstream outfile;
           outfile.open("data/3e.txt");
           for(int i=0; i <= n; i++) {
                solver ES(G, timestep, Earth, Sun);
                ES.VelocityVerlet(G, timestep, Earth, Sun);
                // solver SE(G, timestep, Earth, Sun);
                // SE.VelocityVerlet(G, timestep, Sun, Earth)

                outfile << Earth.position[0] << " ";
                outfile << Earth.position[1] << " ";
                outfile << Jupiter.position[0] << " ";
                outfile << Jupiter.position[1] << endl;
           }
           outfile.close();*/
        ofstream outfile;
        outfile.open("data/3f.txt");
        for(int i=0; i <= n; i++) {
                solver earth(G, timestep, Jupiter, Sun);
                earth.VelocityVerletSystem(G, timestep, Earth, Jupiter, Sun);

                solver jupiter(G, timestep, Earth, Sun);
                jupiter.VelocityVerletSystem(G, timestep, Jupiter, Earth, Sun);

                solver sun(G, timestep, Earth, Jupiter);
                sun.VelocityVerletSystem(G, timestep, Sun, Earth, Jupiter);

                outfile << Earth.position[0] << " ";
                outfile << Earth.position[1] << " ";
                outfile << Jupiter.position[0] << " ";
                outfile << Jupiter.position[1] << " ";
                outfile << Sun.position[0] << " ";
                outfile << Sun.position[1] << " ";
                outfile << Sun.angularMomentum(Earth, Jupiter) << " ";
                outfile << Earth.angularMomentum(Sun, Jupiter) << " ";
                outfile << Jupiter.angularMomentum(Earth, Sun) << endl;
        }
        outfile.close();
};
