#include <cmath>
#include <iostream>
#include <fstream>
#include <random>
#include <armadillo>
#include <random>
#include <map>

using namespace std;
using namespace arma;

ofstream outfile;

inline int PB(int index, int spins, int correction) {
        // returns the index of the element
        return (index + spins + correction)%(spins);
}

void Energy(int n, mat &SMat, double &energy);
map<double, double> energies(double T);
void initialize(double &energy, double &magMoment, int n, mat SMat);
void Metropolis(int n, int mcs, double T, map<double, double> energyDiff);
void WriteToFile(double energy, double magMoment);


int main(int argc, char *argv[]){
        int n, mcc;
        double T;

        n = atoi(argv[1]);    // number of spins
        mcc = atoi(argv[2]);  // number of MC cycles
        T = atof(argv[3]);    // temperature [kT/J]

        Metropolis(n, mcc, T, energies(T));
}

void initialize(double &energy, double &magMoment, int n, mat SMat){
        magMoment = pow(double(n),2);

        for(int x = 0; x < n; x++) {
                for(int y = 0; y < n; y++)
                        energy -= SMat(x,y)*(SMat(x, PB(y,n,-1)) + SMat(PB(x,n,-1),y));
        }
}

void Metropolis(int n, int mcc, double T, map<double, double> energyDiff){
        random_device rd;
        mt19937_64 generator(10);
        uniform_real_distribution<> RNG(0.0,1.0);

        mat SMat = ones<mat>(n,n);    // ground state

        double energy = 0;  // energy
        double magMoment = 0; // magnetic moment (initial value, of ground state = n²)

        // initializing lattice
        initialize(energy, magMoment, n, SMat);
        // calculating energy of lattice
        Energy(n, SMat, energy);

        for(int m = 1; m < mcc; m++) {
                for (int x = 0; x < n; x++) {
                        for (int y = 0; y < n; y++) {
                                int xr = (int) RNG(generator)*(double) n; // indices for random element
                                int yr = (int) RNG(generator)*(double) n;

                                int dE = 2*SMat(xr,yr)*(SMat(xr, PB(yr, n, 1)) + SMat(xr, PB(yr, n, -1))
                                                        + SMat(PB(xr, n, 1),yr) + SMat(PB(xr, n, -1),yr));

                                if (RNG(generator) <= energyDiff.find(dE)->second) {
                                        SMat(xr,yr) *= -1.0;   // flipping spin
                                        magMoment += (double) 2.0*SMat(xr,yr);
                                        energy += (double) dE;
                                }
                        }
                }
                WriteToFile(energy, magMoment);
        } //mc e
}

void Energy(int n, mat &SMat, double &energy){
        // initial energy of lattice
        for(int x = 0; x < n; x++) {
                for(int y = 0; y < n; y++) {
                        energy -= SMat(x, y)*(SMat(PB(x, n, -1), y) + SMat(PB(x, n, 1), y)
                                              + SMat(x, PB(y, n, -1)) + SMat(x, (y, n, 1)));
                }
        }
}

map<double, double> energies(double T){
        /* for comparing the accepted values of the energy with a random number
           in the Metropolis algorithm. Compare key dE to value of energy. */
        map<double, double> acceptE;
        for (int dE = -8; dE <= 8; dE+=4) {
                acceptE.insert(pair<double, double>(dE, exp(-dE/T)));
        }
        return acceptE;
}

void WriteToFile(double energy, double magMoment){
        outfile.open("data/mean_vals.txt", fstream::app);
        outfile << energy << " ";
        outfile << energy*energy << " ";
        outfile << magMoment << " ";
        outfile << magMoment*magMoment << endl;

        outfile.close();
}





















//
