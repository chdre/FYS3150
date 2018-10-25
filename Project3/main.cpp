#include <iostream>
#include <cmath>
#include <armadillo>
#include "planet.h"
#include "solver.h"

using namespace std;

void Write3cToFile(int method, double G, double h, int n, double M_e, double M_sun);
void WriteEnergynMomentumToFile(double G, double h, int n, double M_sun, double M_e);
void Write3eToFile(double G, double h, int n, double M_e, double M_j, double M_sun);
void Write3dToFile(double G, double h, double beta, int n, double M_sun, double M_e);
void Write3fToFile(double G, double h, int n, double M_e, double M_m, double M_sun);
void Write3gToFile(double G, double h, int n, double c, double M_sun, double M_m);

int main(){
        double M_e, M_sun, M_j, M_sunval, M_m, v_sun0, v_j0, v_e0, initPos_sun;
        double totalMass, CenterOfMass, G, Time, timestep, vfac, v0_e, c, beta;
        int n;

        Time = 2.0;    // time [years]
        n = 1e3;    // steps
        timestep = Time/n;
        beta = 3;

        // constants
        G = 4.0*pow(M_PI,2);        // gravitational constant [AU³/yr²]
        c = 63240.08034;      // speed of light [AU/yr]

        // masses
        M_sunval = 2.0e30;        // mass of sun
        M_sun = 1.0;              // relative mass of Sun
        M_e = 6.0e24/M_sunval;    // relative mass of Earth
        M_j = 1.9e27/M_sunval;    // relative mass of Jupiter
        M_m = 3.3e23/M_sunval;    // relative mass of Mercury
        totalMass = M_sun + M_e + M_j;

        v_e0 = 2.0*M_PI;    // initial velocity of Earth
        v_j0 = 0.7*M_PI;    // initial velocity of Jupiter
        v_sun0 = -(v_e0*M_e + v_j0*M_j)/M_sun;  // initial velocuty of Sun

        
        // Below are the calls for the functions to write to file. Remove // to call
        //Write3cToFile(1, G, timestep, n, M_e, M_sun);
        WriteEnergynMomentumToFile(G, timestep, n, M_sun, M_e);
        //Write3dToFile(G, timestep, beta, n, M_sun, M_e);
        //Write3eToFile(G, timestep, n, M_e, M_j, M_sun);
        //Write3fToFile(G, timestep, n, M_e, M_m, M_sun);
        //Write3gToFile(G, timestep, n, c, M_sun, M_m);
}

void WriteEnergynMomentumToFile(double G, double h, int n, double M_sun, double M_e){
        /* Writing the enery and angular momentum off the Sun and Earth system. Using Velocity Verlet. */
        planet Sun(M_sun, 0, 0, 0, 0, 0, 0);
        planet Earth(M_e, 1, 0, 0, 0, 2.0*M_PI, 0);

        ofstream outfile;
        outfile.open("data/3c-energy.txt");
        for(int i=0; i <= n; i++) {
                solver solve(G, h, Earth, Sun);
                solve.VelocityVerlet(G, h, Earth, Sun);
                Earth.potential = Earth.potentialEnergy(G,Sun);
                Earth.kinetic = Earth.kineticEnergy();
                Earth.angularMom = Earth.angularMomentum(Sun);


                outfile << Earth.kinetic << " ";
                outfile << Earth.potential << " ";
                outfile << Earth.angularMom << " ";
                outfile << h*i << endl;
        }
        outfile.close();
}

void Write3cToFile(int method, double G, double h, int n, double M_e, double M_sun){
        /* method = 1 => Verlet. method = 2 = > euler. 
        Writing the position of the Earth to file using either Euler or Velocity Verlet.*/
        planet Sun(M_sun, 0, 0, 0, 0, 0, 0);
        planet Earth(M_e, 1, 0, 0, 0, 2.0*M_PI, 0);

        ofstream outfile;
        outfile.open("data/3c-verlet.txt");
        for(int i=0; i <= n; i++) {
                if(method == 1) {
                        // Velocity Verlet
                        solver earth(G, h, Earth, Sun);
                        earth.VelocityVerlet(G, h, Earth, Sun);
                        solver sun(G, h, Sun, Earth);
                        sun.VelocityVerlet(G, h, Sun, Earth);
                }
                if(method == 2) {
                        // Euler
                        solver earth(G, h, Earth, Sun);
                        earth.euler(G, h, Earth, Sun);
                        solver sun(G, h, Sun, Earth);
                        sun.euler(G, h, Sun, Earth);
                }

                outfile << Earth.position[0] << " ";
                outfile << Earth.position[1] << endl;
        }
        outfile.close();
};

void Write3dToFile(double G, double h, double beta, int n, double M_sun, double M_e) {
        /* Writing Earth position to file using the alternative force with 1/r^{beta}. Velocity Verlet.*/
        planet Sun(M_sun, 0, 0, 0, 0, 0, 0);
        planet Earth(M_e, 1, 0, 0, 0, 2.0*M_PI, 0);

        ofstream outfile;
        outfile.open("data/3d-newforce.txt");
        for(int i=0; i <= n; i++) {
                solver earth(G, h, Earth, Sun);
                earth.VelocityVerletAlt(G, h, beta, Earth, Sun);

                outfile << Earth.position[0] << " ";
                outfile << Earth.position[1] << " ";

        }
        outfile.close();
}

void Write3eToFile(double G, double h, int n, double M_e, double M_j, double M_sun){
        /* Writing Earth and  Jupiter's position to file. Three body problem with Earth fixed.
        Using Velocity Verlet. */
        double v0_j = sqrt(4.0*pow(M_PI,2)/5.2);        // initial velocity of Jupiter
        double mass_factor = 1000;      // factor to change Jupiter's mass

        planet Sun(M_sun, 0, 0, 0, 0, 0, 0);
        planet Earth(M_e, 1, 0, 0, 0, 2.0*M_PI, 0);
        planet Jupiter(M_j*mass_factor, 5.2, 0, 0, 0, v0_j, 0);

        ofstream outfile;
        outfile.open("data/3e.txt");
        int counter = 0;        // counter to not write all positions to file, if n is large
        for(int i=0; i <= n; i++) {
                solver earth(G, h, Jupiter, Sun);
                earth.VelocityVerletSystem(G, h, Earth, Jupiter, Sun);

                solver jupiter(G, h, Earth, Sun);
                jupiter.VelocityVerletSystem(G, h, Jupiter, Earth, Sun);

                if (i == counter) {
                        outfile << Earth.position[0] << " ";
                        outfile << Earth.position[1] << " ";

                        outfile << Jupiter.position[0] << " ";
                        outfile << Jupiter.position[1] << endl;
                        counter += 1;   // set this larger if n > 1e4
                }
        }
        outfile.close();
}

void Write3fToFile(double G, double h, int n, double M_e, double M_j, double M_sun){
        /* Writing positon of Earth, Jupiter and Sun (no longer fixex) to file. Calculating using 
        Velocity Verlet. */
        double v0_j = sqrt(4.0*pow(M_PI,2)/5.2);        // initial velocity of Jupiter
        double totalMass = M_j + M_sun + M_e;           // total mass of solar system
        double initPos_sun = -(5.2*M_j + M_e)/(totalMass*M_sun);   // initial position of Sun found by requiring center of mass = 0
        double initVel_sun = -(M_e*2*M_PI + M_j*v0_j)/M_sun;       // initial velocity of Sun found by requiring sum of momentum = 0

        planet Sun(M_sun, initPos_sun, 0, 0, 0, 0, 0);
        planet Earth(M_e, 1, 0, 0, 0, 2.0*M_PI, 0);
        planet Jupiter(M_j, 5.2, 0, 0, 0, v0_j, 0);

        ofstream outfile;
        outfile.open("data/3f.txt");
        int counter = 0;
        for(int i=0; i <= n; i++) {
                solver earth(G, h, Jupiter, Sun);
                earth.VelocityVerletSystem(G, h, Earth, Jupiter, Sun);

                solver jupiter(G, h, Earth, Sun);
                jupiter.VelocityVerletSystem(G, h, Jupiter, Earth, Sun);

                solver sun(G, h, Earth, Jupiter);
                sun.VelocityVerletSystem(G, h, Sun, Earth, Jupiter);

                if (i == counter) {
                        outfile << Earth.position[0] << " ";
                        outfile << Earth.position[1] << " ";

                        outfile << Jupiter.position[0] << " ";
                        outfile << Jupiter.position[1] << " ";

                        outfile << Sun.position[0] << " ";
                        outfile << Sun.position[1] << endl;

                        counter += 100;
                }
        }
        outfile.close();
}

void Write3gToFile(double G, double h, int n, double c, double M_sun, double M_m){
        /* Writing positions of Mercury and Sun to file, OR precession of Mercury in arcseconds.
        Main idea: 
        1) calculate old and old old radius and old position
        2) update and calculate new radius
        3) check if the old radius is the lowest by requiring that the new radius and the old old radius is larger */
        double r_old, r_old2, r_new;
        vec r_pos;
        vector<double> arcsec;

        planet Sun(M_sun, 0, 0, 0, 0, 0, 0);
        planet Mercury(M_m, 0.3075, 0, 0, 0, 12.44, 0);
        //planet Mercury(1.65e-7, 2.268331723540750E-02, -4.527628626217761E-01, -3.976201975514675E-02,
        //               2.244877340585597E-02*365.25, 2.845538482588606E-03*365.25, -1.827605072960171E-03*365.25);

        ofstream outfile;
        outfile.open("data/3g-data.txt");

        r_old2 = Mercury.distance(Sun);
        for(int i=0; i <= n; i++) {
                r_old = Mercury.distance(Sun);
                r_pos = Mercury.position;

                solver mercury(G, h, Mercury, Sun);
                mercury.VelocityVerletEinstein(G, c, h, Mercury, Sun);

                solver sun(G,h,Mercury,Sun);
                sun.VelocityVerletEinstein(G, c, h, Sun, Mercury);

                r_new = Mercury.distance(Sun);
                if (r_old < r_old2 and r_new > r_old) {
                        arcsec.push_back(atan(r_pos[1]/r_pos[0])*206265.806);
                }
                r_old2 = r_old;
                
                // unlock below to write position to file, if interested.
                /*int counter = 0;
                   if (i == counter) {
                        outfile << Mercury.position[0] << " ";
                        outfile << Mercury.position[1] << endl;

                        counter += 1000;
                   }
                   outfile.close();*/

        }
        for (int i = 0; i < arcsec.size(); i++) {
                outfile << arcsec[i] << endl;
                outfile.close();
        }
};
