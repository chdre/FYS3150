#include <cmath>
#include <iostream>
#include <fstream>
#include <random>
#include <armadillo>
#include <random>
#include <map>

using namespace std;
using namespace arma;

inline int PB(int index, int spins, int correction) {
        // returns the index of the element
        return (index + spins + correction)%(spins);
}

void Energy(int n, mat &SMat, double &energy);
void Metropolis(int n, int mcs, double T);
map<double, double> energies(double T);


int main(int argc, char *argv[]){
        int n, mcc;
        double T;

        n = atoi(argv[1]);    // number of spins
        mcc = atoi(argv[2]);  // number of MC cycles
        T = atof(argv[3]);    // temperature [kT/J]

}

void Metropolis(int n, int mcc, double T, map<double, double> energyDiff){
        random_device rd;
        mt19937_64 generator(rd());
        uniform_int_distribution<> random_spin(0, 1);
        uniform_real_distribution<> RNG(0.0,1.0);

        mat SMat = ones<mat>(n,n);    // ground state
        mat ExpectVals = zeros(n,5);  // expectation values

        double energy = 0;  // energy
        double magMoment = pow(double(n),2); // magnetic moment (initial value, of ground state = n²)

        // calculating energy of lattice
        Energy(n, SMat, energy);

        for(int m = 1; m < mcc; m++) {
                for (int x = 0; x < n; x++) {
                        for (int y = 0; y < n; y++) {
                                SMat(x,y) = 2*random_spin(generator) - 1;   // randomizing spin at x,y
                                int dE = 2*SMat(x,y)*(SMat(x, PB(x, n, 1)) + SMat(x, PB(x, n, -1))
                                                      + SMat(PB(y, n, 1),y) + SMat(PB(y, n, -1),y));

                                if (RNG(generator) <= energyDiff.find(dE)->second) {
                                        SMat(x,y) *= -1;   // flipping spin
                                        magMoment += (double) 2.0*SMat(x,y);
                                        energy += (double) dE;
                                }
                                else{
                                        if (RNG(generator) <= energyDiff.find(dE)->second) {
                                                SMat(x,y) *= -1; // flipping spin
                                                magMoment += (double) 2.0*SMat(x,y);
                                                energy += (double) dE;
                                        }
                                }
                        }
                }
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
























//
