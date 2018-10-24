#include <iostream>
#include <cmath>
#include <armadillo>
#include "planet.h"
#include "solver.h"

using namespace std;

void Write3cToFile(string method, double G, double h, int n, double M_e, double M_sun);
void Write3eToFile(double G, double h, int n, double M_e, double M_j, double M_sun);
void Write3dToFile(double G, double h, double beta, int n, double M_sun, double M_e);
void Write3fToFile(double G, double h, int n, double M_e, double M_m, double M_sun);
void Write3gToFile(double G, double h, int n, double c, double M_sun, double M_m);

int main(){
        double M_e, M_sun, M_j, M_sunval, M_m, v_sun0, v_j0, v_e0, initPos_sun;
        double totalMass, CenterOfMass, G, Time, timestep, vfac, v0_e, c, beta;
        int n;

        Time = 100.0;    // time [years]
        n = 1e9;    // steps
        timestep = Time/n;
        beta = 2.0;

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

        //initPos_sun = -(5.2*M_j + M_e)/(totalMass*M_sun);   // initial position of Sun found by requiring center of mass = 0
        string euler;
        Write3cToFile(euler, G, timestep, n, M_e, M_sun);
        //Write3eToFile(G, timestep, n, M_e, M_j, M_sun);
        //Write3dToFile(G, timestep, beta, n, M_sun, M_e);
        //Write3fToFile(G, timestep, n, M_e, M_m, M_sun);
        //Write3gToFile(G, timestep, n, c, M_sun, M_m)M
}


void Write3cToFile(string method, double G, double h, int n, double M_e, double M_sun){
        planet Sun(M_sun, 0, 0, 0, 0, 0, 0);
        planet Earth(M_e, 1, 0, 0, 0, 2.0*M_PI, 0);

        ofstream outfile;
        outfile.open("data/testing.txt");
        for(int i=0; i <= n; i++) {
                if(method == "Verlet") {
                        solver earth(G, h, Earth, Sun);
                        earth.VelocityVerlet(G, h, Earth, Sun);
                        solver sun(G, h, Sun, Earth);
                        sun.VelocityVerlet(G, h, Sun, Earth);
                }
                if(method == "Euler") {
                        solver earth(G, h, Earth, Sun);
                        earth.euler(G, h, Earth, Sun);
                        solver sun(G, h, Sun, Earth);
                        sun.euler(G, h, Sun, Earth);
                }

                outfile << Earth.position[0] << " ";
                outfile << Earth.position[1] << " ";
                outfile << Sun.position[0] << " ";
                outfile << Sun.position[1] << endl;
        }
        outfile.close();
};

void Write3dToFile(double G, double h, double beta, int n, double M_sun, double M_e) {
        planet Sun(M_sun, 0, 0, 0, 0, 0, 0);
        planet Earth(M_e, 1, 0, 0, 0, 1.41*2.0*M_PI, 0);

        ofstream outfile;
        outfile.open("data/testing.txt");
        for(int i=0; i <= n; i++) {
                solver earth(G, h, Earth, Sun);
                earth.VelocityVerletAlt(G, h, beta, Earth, Sun);
                solver sun(G, h, Sun, Earth);
                sun.VelocityVerletAlt(G, h, beta, Sun, Earth);

                outfile << Earth.position[0] << " ";
                outfile << Earth.position[1] << " ";
                outfile << Sun.position[0] << " ";
                outfile << Sun.position[1] << endl;
        }
        outfile.close();
}

void Write3eToFile(double G, double h, int n, double M_e, double M_j, double M_sun){
        double v0_j = sqrt(4.0*pow(M_PI,2)/5.2);
        planet Sun(M_sun, 0, 0, 0, 0, 0, 0);
        planet Earth(M_e, 1, 0, 0, 0, 2.0*M_PI, 0);
        planet Jupiter(M_j, 5.2, 0, 0, 0, v0_j, 0);

        ofstream outfile;
        outfile.open("data/testing.txt");
        double counter = 0;
        for(int i=0; i <= n; i++) {
                if(i=counter) {
                        solver earth(G, h, Jupiter, Sun);
                        earth.VelocityVerletSystem(G, h, Earth, Jupiter, Sun);

                        solver jupiter(G, h, Earth, Sun);
                        jupiter.VelocityVerletSystem(G, h, Jupiter, Earth, Sun);

                        solver sun(G, h, Earth, Jupiter);
                        sun.VelocityVerletSystem(G, h, Sun, Earth, Jupiter);

                        //vec CoM = (M_e*Earth.position + M_j*Jupiter.position + M_sun*Sun.position)/totalMass;

                        outfile << Earth.position[0] << " ";
                        outfile << Earth.position[1] << " ";

                        outfile << Jupiter.position[0] << " ";
                        outfile << Jupiter.position[1] << " ";

                        outfile << Sun.position[0] << " ";
                        outfile << Sun.position[1] << endl;
                }
                counter += 100;
        }
        outfile.close();
}

void Write3fToFile(double G, double h, int n, double M_e, double M_m, double M_sun){
        double v0_j = sqrt(4.0*pow(M_PI,2)/5.2);
        planet Sun(M_sun, 0, 0, 0, 0, 0, 0);
        planet Earth(M_e, 1, 0, 0, 0, 2.0*M_PI, 0);
        planet Mercury(M_m, 5.2, 0, 0, 0, v0_j, 0);
        ofstream outfile;
        outfile.open("data/testing.txt");
        for(int i=0; i <= n; i++) {
                solver earth(G, h, Mercury, Sun);
                earth.VelocityVerletSystem(G, h, Earth, Mercury, Sun);

                solver mercury(G, h, Earth, Sun);
                mercury.VelocityVerletSystem(G, h, Mercury, Earth, Sun);

                solver sun(G, h, Earth, Mercury);
                sun.VelocityVerletSystem(G, h, Sun, Earth, Mercury);

                //vec CoM = (M_e*Earth.position + M_j*Jupiter.position + M_sun*Sun.position)/totalMass;

                outfile << Earth.position[0] << " ";
                outfile << Earth.position[1] << " ";

                outfile << Mercury.position[0] << " ";
                outfile << Mercury.position[1] << " ";

                outfile << Sun.position[0] << " ";
                outfile << Sun.position[1] << endl;
        }
        outfile.close();
}


void Write3gToFile(double G, double h, int n, double c, double M_sun, double M_m){
        double r_old, r_old2, r_new;
        vec r_pos;
        vector<double> arcsec;

        planet Sun(M_sun, 0, 0, 0, 0, 0, 0);
        planet Mercury(M_m, 0.3075, 0, 0, 0, 12.44, 0);

        ofstream outfile;
        outfile.open("data/testing.txt");

        r_old2 = Mercury.distance(Sun);
        for(int i=0; i <= n; i++) {
                r_old = Mercury.distance(Sun);
                r_pos = Mercury.position;
                solver mercury(G, h, Mercury, Sun);
                mercury.VelocityVerletEinstein(G, c, h, Mercury, Sun);

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
