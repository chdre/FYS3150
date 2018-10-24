#include <iostream>
#include <cmath>
#include <armadillo>
#include "planet.h"
#include "solver.h"

using namespace std;

int main(){
        double M_e, M_sun, M_j, M_sunval, M_m, v_sun0, v_j0, v_e0, initPos_sun;
        double totalMass, CenterOfMass, G, Time, timestep, vfac, v0_e, c;
        double r_old2, r_old, r_new;
        int n;
        vec r_pos;

        vector<double> arcsec;

        Time = 100.0;    // time [years]
        n = 1e9;    // steps
        timestep = Time/n;

        // constants
        G = 4.0*pow(M_PI,2);    // gravitational constant [AU³/yr²]

        // masses
        M_sunval = 2.0e30;          // mass of sun
        M_sun = 1.0;              // relative mass of Sun
        M_e = 6.0e24/M_sunval;    // relative mass of Earth
        M_j = 1.9e27/M_sunval;    // relative mass of Jupiter
        M_m = 3.3e23/M_sunval;    // relative mass of Mercury
        totalMass = M_sun + M_e + M_j;

        v_e0 = 2.0*M_PI;    // initial velocity of Earth
        v_j0 = 0.7*M_PI;    // initial velocity of Jupiter
        v_sun0 = -(v_e0*M_e + v_j0*M_j)/M_sun;  // initial velocuty of Sun

        c = 63240.08034;   // speed of light [AU/yr]


        //initPos_sun = -(5.2*M_j + M_e)/(totalMass*M_sun);   // initial position of Sun found by requiring center of mass = 0

        planet Sun(M_sun, 0, 0, 0, 0, 0, 0);
        planet Mercury(M_m, 0.3075, 0, 0, 0, 12.44, 0);


        /*ofstream outfile;
           outfile.open("data/3g.txt");
           for(int i=0; i <= n; i++) {
                solver mercury(G, timestep, Mercury, Sun);
                ES.VelocityVerlet(G, timestep, Mercury, Sun);
                solver sun(G, timestep, Sun, Mercury);
                SE.VelocityVerlet(G, timestep, Sun, Mercury);

                outfile << Earth.position[0] << " ";
                outfile << Earth.position[1] << " ";
                outfile << Mercury.position[0] << " ";
                outfile << Mercury.position[1] << endl;
           }
           outfile.close();*/
        /*ofstream outfile;
           outfile.open("data/3g.txt");
           double counter = 0;
           for(int i=0; i <= n; i++) {
                if(i=counter) {
                        solver earth(G, timestep, Mercury, Sun);
                        earth.VelocityVerletSystem(G, timestep, Earth, Mercury, Sun);

                        solver mercury(G, timestep, Earth, Sun);
                        mercury.VelocityVerletSystem(G, timestep, Mercury, Earth, Sun);

                        solver sun(G, timestep, Earth, Mercury);
                        sun.VelocityVerletSystem(G, timestep, Sun, Earth, Mercury);

                        //vec CoM = (M_e*Earth.position + M_j*Jupiter.position + M_sun*Sun.position)/totalMass;

                        outfile << Earth.position[0] << " ";
                        outfile << Earth.position[1] << " ";
                        // outfile << Earth.position[2] << " ";

                        outfile << Mercury.position[0] << " ";
                        outfile << Mercury.position[1] << " ";
                        // outfile << Mercury.position[2] << " ";

                        outfile << Sun.position[0] << " ";
                        outfile << Sun.position[1] << endl;
                        // outfile << Sun.position[2] << endl;
                }
                counter += 100;

           }
           outfile.close();*/
        ofstream outfile;
        outfile.open("data/3g.txt");

        r_old2 = Mercury.distance(Sun);
        for(int i=0; i <= n; i++) {
                r_old = Mercury.distance(Sun);
                r_pos = Mercury.position;
                solver mercury(G, timestep, Mercury, Sun);
                mercury.VelocityVerletEinstein(G, c, timestep, Mercury, Sun);

                r_new = Mercury.distance(Sun);

                if (r_old < r_old2 and r_new > r_old) {
                        arcsec.push_back(atan(r_pos[1]/r_pos[0])*206265.806);
                }
                r_old2 = r_old;



                /*outfile << Mercury.position[0] << " ";
                   outfile << Mercury.position[1] << " ";
                   outfile << Sun.position[0] << " ";
                   outfile << Sun.position[1] << endl;*/
        }

        outfile.close();
        for (int i = 0; i < arcsec.size(); i++) {
                cout << arcsec[i] << endl;
        }
};
