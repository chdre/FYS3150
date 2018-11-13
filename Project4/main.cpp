#include "mpi.h"
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
void initialize(double &energy, double &magMoment, int n, mat S, int GS);
void Metropolis(int n, int mcs, double T, map<double, double> acceptE, int numprocs, int my_rank, int myloop_end, int myloop_begin);
map<double, double> energies(double T);
void WriteToFile(mat EnergyMagSave, int mcc, double T);
void AddExpectValsToVec(vec &ExpectVals, double energy, double magMoment);
void PrintExpectVals(vec TotalExpectVals, int mcc, double T);

int main(int argc, char *argv[]){
        int n, mcc, numprocs, my_rank;
        double T;

        T = 1.0;    // temperature [kT/J]

        n = atoi(argv[1]);    // number of spins
        mcc = atoi(argv[2]);  // number of MC cycles

        //  MPI initializations
        MPI_Init (&argc, &argv);
        MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
        MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);

        int no_intervals = mcc/numprocs;
        int myloop_begin = my_rank*no_intervals + 1;
        int myloop_end = (my_rank + 1)*no_intervals;

        if((my_rank == numprocs - 1) && (myloop_end < mcc)) {
                myloop_end = mcc;
        }

        double TimeStart, TimeEnd, TotalTime;
        TimeStart = MPI_Wtime();

        Metropolis(n, mcc, T, energies(T), numprocs, my_rank, myloop_end, myloop_begin);

        TimeEnd = MPI_Wtime();
        TotalTime = TimeEnd - TimeStart;
        if(my_rank == 0) {
                cout << "Time = " << TotalTime << " on number of threads: " << numprocs << endl;
        }
        MPI_Finalize();
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

void Metropolis(int n, int mcc, double T, map<double, double> acceptE, int numprocs, int my_rank, int myloop_end, int myloop_begin){
        random_device rd;
        mt19937_64 generator(rd());
        uniform_real_distribution<float> RNG(0.0,1.0);
        uniform_int_distribution<int> RNGPos(0,n-1);

        double energy = 0.0;  // energy
        double magMoment = 0.0; // magnetic moment (initial value, of ground state = n²)

        vec ExpectVals = zeros<vec>(5);
        vec TotalExpectVals = zeros<vec>(5);
        mat S = zeros<mat>(n,n);
        mat EnergyMagSave = zeros<mat>(mcc,2);     // for storing energy and magnetic moment values
        mat EnergyMagSaveTot = zeros<mat>(mcc,2);     // for storing energy and magnetic moment values


        // initializing lattice
        initialize(energy, magMoment, n, S, 1, generator);

        if(my_rank == 0) {
                outfile.open("data/test.txt", fstream::app);
        }

        int counter = 0;
        for(int m = myloop_begin; m <= myloop_end; m++) {
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
                AddExpectValsToVec(ExpectVals, energy, magMoment);
                EnergyMagSave(m-1,0) = energy;
                EnergyMagSave(m-1,1) = magMoment;

        }  //mc end

        for(int i = 0; i < 5; i++) {
                // Summing the expectvals
                MPI_Reduce(&ExpectVals(i), &TotalExpectVals(i), 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        }

        if(my_rank == 0) {
                // only 1 print to file
                //EnergyMagSave = EnegyMagSave/mcc;
                EnergyMagSave.save("data/test.txt", arma_ascii);
                //  WriteToFile(EnergyMagSave, mcc, T);

        }
        if(my_rank == 0) {
                outfile.close();
                PrintExpectVals(TotalExpectVals, mcc, T);
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

void WriteToFile(mat EnergyMagSave, int mcc, double T){
        double nfac = 1.0/mcc;

        for(int i = 0; i < mcc; i++) {
                double E = EnergyMagSave(i,0)*nfac;
                double E2 = EnergyMagSave(i,0)*EnergyMagSave(i,0)*nfac;
                double M = EnergyMagSave(i,1)*nfac;
                double M2 = EnergyMagSave(i,1)*EnergyMagSave(i,1)*nfac;
                double Mabs = fabs(EnergyMagSave(i,1))*nfac;

                double C_V = (E2 - E)/pow(T,2);//*pow(n,2));
                double chi = (M2 - M)/T;//*pow(n,2));

                outfile << E << " ";
                outfile << M << " ";
                outfile << C_V << " ";
                outfile << chi << " ";
                outfile << Mabs << endl;
        }
}

void AddExpectValsToVec(vec &ExpectVals, double energy, double magMoment){
        ExpectVals(0) += energy;
        ExpectVals(1) += energy*energy;
        ExpectVals(2) += magMoment;
        ExpectVals(3) += magMoment*magMoment;
        ExpectVals(4) += fabs(magMoment);
}

void PrintExpectVals(vec TotalExpectVals, int mcc, double T) {
        TotalExpectVals = TotalExpectVals/mcc;
        cout << "E: " << TotalExpectVals(0) << " ";
        cout << "Mabs: " << TotalExpectVals(4) << " ";
        cout << "C_V: " << (TotalExpectVals(1) - pow(TotalExpectVals(0),2))/(pow(T,2)) << " ";
        cout << "chi: " << (TotalExpectVals(3) - pow(TotalExpectVals(2),2))/T << endl;
}
















//
