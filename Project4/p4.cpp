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
void WriteToFile(double energy, double magMoment, int mcc, double T, int n);


int main(int argc, char *argv[]){
        int n, mcc;
        double T;

        n = atoi(argv[1]);    // number of spins
        mcc = atoi(argv[2]);  // number of MC cycles
        T = atof(argv[3]);    // temperature [kT/J]

        Metropolis(n, mcc, T, energies(T));
}

void initialize(double &energy, double &magMoment, int n, mat SMat){
        magMoment += pow(double(n),2);  // initial all spin up

        for(int x = 0; x < n; x++) {
                for(int y = 0; y < n; y++)
                        energy -= (double) SMat(x,y)*(SMat(x, PB(y,n,-1)) + SMat(PB(x,n,-1),y));
        }
}

void Metropolis(int n, int mcc, double T, map<double, double> acceptE){
        random_device rd;
        mt19937_64 generator(10);
        uniform_real_distribution<float> RNG(0.0,1.0);
        uniform_int_distribution<int> RNGSpin(0,n-1);

        mat SMat = ones<mat>(n,n);    // ground state

        double energy = 0.0;  // energy
        double magMoment = 0.0; // magnetic moment (initial value, of ground state = n²)

        // initializing lattice
        initialize(energy, magMoment, n, SMat);

        int counter = 0;
        for(int m = 1; m < mcc; m++) {
                for (int x = 0; x < n; x++) {
                        for (int y = 0; y < n; y++) {
                                int xr = RNGSpin(generator);//int xr = (int) (RNG(generator)*(double) n); // indices for random element
                                int yr = RNGSpin(generator);//int yr = (int) (RNG(generator)*(double) n);//;

                                int dE = 2.0*SMat(xr,yr)*(SMat(xr,PB(yr,n,1)) + SMat(xr,PB(yr,n,-1))
                                                          + SMat(PB(xr,n,1),yr) + SMat(PB(xr,n,-1),yr));

                                cout << acceptE.find(dE)->second << endl;

                                if (RNG(generator) <= acceptE.find(dE)->second) {
                                        counter += 1;
                                        SMat(xr,yr) *= -1.0; // flipping spin
                                        magMoment += (double) 2*SMat(xr,yr);
                                        energy += (double) dE;
                                }
                        }
                }
                WriteToFile(energy, magMoment, mcc, T, n);
        } //mc e
        cout << counter << endl;
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

void WriteToFile(double energy, double magMoment, int mcc, double T, int n){
        outfile.open("data/mean_vals.txt", fstream::app);
        double nfac = 1.0/mcc;

        double E = energy;//*nfac;
        double E2 = E*E;
        double M = magMoment;//*nfac;
        double M2 = M*M;

        double C_V = (E2 - E)/(pow(T,2));//*pow(n,2));
        double chi = (M2 - M)/(T);//*pow(n,2));


        outfile << E << " ";
        outfile << M << " ";
        outfile << C_V << " ";
        outfile << chi << endl;
        outfile.close();

}





















//
