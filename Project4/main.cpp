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

void Energy(int n, mat &S, double &energy);
map<double, double> energies(double T);
void initialize(double &energy, double &magMoment, int n, mat S);
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

void initialize(double &energy, double &magMoment, int n, mat &S, int GS, mt19937_64 generator){
        if (GS==1) {
                S.fill(1.0); // ground state
                magMoment += pow(double(n),2); // initial all spin up
        }
        else {
                uniform_int_distribution<int> RNGSpin(0,1);
                for(int x = 0; x < n; x++) {
                        for (int y = 0; y < n; y++) {
                                S(x,y) = (double) 2*RNGSpin(generator) - 1; // adding random spin states
                        }
                }
        }

        for(int x = 0; x < n; x++) {
                for(int y = 0; y < n; y++) {
                        energy -= (double) S(x,y)*(S(x,PB(y,n,-1)) + S(PB(x,n,-1),y));
                        magMoment += (double) S(x,y);
                }
        }
}

void Metropolis(int n, int mcc, double T, map<double, double> acceptE){
        random_device rd;
        mt19937_64 generator(rd());
        uniform_real_distribution<float> RNG(0.0,1.0);
        uniform_int_distribution<int> RNGPos(0,n-1);

        double energy = 0.0;  // energy
        double magMoment = 0.0; // magnetic moment (initial value, of ground state = n²)

        vec ExpectVals = zeros<vec>(5);
        mat S = zeros<mat>(n,n);

        // initializing lattice
        initialize(energy, magMoment, n, S, 1, generator);

        outfile.open("data/test.txt", fstream::app);

        int counter = 0;
        for(int m = 1; m <= mcc; m++) {
                for (int x = 0; x < n; x++) {
                        for (int y = 0; y < n; y++) {
                                int xr = RNGPos(generator); // indices for random element int) (RNG(generator)*(double) n); //
                                int yr = RNGPos(generator); // (int) (RNG(generator)*(double) n); // RNGSpin(generator);//;

                                int dE = 2.0*S(xr,yr)*(S(xr,PB(yr,n,1)) + S(xr,PB(yr,n,-1))
                                                       + S(PB(xr,n,1),yr) + S(PB(xr,n,-1),yr));
                                if (RNG(generator) <= acceptE.find(dE)->second) {
                                        counter += 1;
                                        S(xr,yr) *= -1.0; // flipping spin
                                        magMoment += (double) 2*S(xr,yr);
                                        energy += (double) dE;
                                }
                        }
                }
                WriteToFile(energy, magMoment, mcc, T, n);
                ExpectVals(0) += energy;
                ExpectVals(1) += energy*energy;
                ExpectVals(2) += magMoment;
                ExpectVals(3) += magMoment*magMoment;
                ExpectVals(4) += fabs(magMoment);

        } //mc e
        outfile.close();
        ExpectVals = ExpectVals/mcc;
        double temp = 1.0;//*pow(n,2);
        cout << "E: " << ExpectVals(0)*temp << " ";
        cout << "M: " << ExpectVals(4)*temp << " ";
        cout << "C_V: " << (ExpectVals(1) - pow(ExpectVals(0),2))/(pow(T,2))*temp << " ";
        cout << "chi: " << (ExpectVals(3) - pow(ExpectVals(2),2))/T*temp << endl;
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
        double nfac = 1.0/mcc;

        double E = energy*nfac;
        double E2 = energy*energy*nfac;
        double M = magMoment*nfac;
        double M2 = magMoment*magMoment*nfac;
        double Mabs = fabs(magMoment)*nfac;

        double C_V = (E2 - E)/(pow(T,2));//*pow(n,2));
        double chi = (M2 - M)/(T);//*pow(n,2));

        outfile << E << " ";
        outfile << M << " ";
        outfile << C_V << " ";
        outfile << chi << " ";
        outfile << Mabs << endl;
}





















//
