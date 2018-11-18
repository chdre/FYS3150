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

void Energy(int L, mat &S, double &energy);
void initialize(double &energy, double &magMoment, int L, mat S, int GS);
void Metropolis(int L, int mcs, double T, map<double, double> acceptE, int numprocs, int my_rank, int myloop_end, int myloop_begin, string filename);
map<double, double> energies(double T);
void WriteToFile(double energy, double magMoment, int mcc, double T, int accepts, int L);
void AddExpectValsToVec(vec &ExpectVals, double energy, double magMoment);
void PrintExpectVals(vec TotalExpectVals, int mcc, double T);
void WriteExpectValsToFile(vec TotalExpectVals, int mcc, double T, int L);


int main(int argc, char *argv[]){
        int L, mcc, numprocs, my_rank;
        //double T;
        string filename;

        //T = atoi(argv[4]);    // temperature [kT/J]
        L = atoi(argv[1]);    // number of spins
        mcc = atoi(argv[2]);  // number of MC cycles
        filename = argv[3];   // name of file

        //  MPI initializations
        MPI_Init (&argc, &argv);
        MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
        MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);

        int no_intervals = mcc/numprocs;
        int myloop_begin = my_rank*no_intervals + 1;
        int myloop_end = (my_rank + 1)*no_intervals;

        if((my_rank == numprocs - 1) && (myloop_end < mcc)) myloop_end = mcc;

        double TimeStart, TimeEnd;
        TimeStart = MPI_Wtime();

        for(double T = 2.0; T <= 2.3; T+=0.03) {
                Metropolis(L, mcc, T, energies(T), numprocs, my_rank, myloop_end, myloop_begin, filename);
        }

        TimeEnd = MPI_Wtime();
        if(my_rank == 0) {
                cout << "Time = " << TimeEnd - TimeStart << " on number of threads: " << numprocs << endl;
        }
        MPI_Finalize();
}

void initialize(double &energy, double &magMoment, int L, mat &S, int GS, mt19937_64 generator){
        if (GS==1) {
                S.fill(1.0); // ground state
                magMoment += pow(double(L),2); // initial all spin up
        }
        else {
                uniform_int_distribution<int> RNGSpin(0,1);
                for(int x = 0; x < L; x++) {
                        for (int y = 0; y < L; y++) {
                                S(x,y) = (double) 2*RNGSpin(generator) - 1; // adding random spin states
                        }
                }
        }

        for(int x = 0; x < L; x++) {
                for(int y = 0; y < L; y++) {
                        energy -= (double) S(x,y)*(S(x,PB(y,L,-1)) + S(PB(x,L,-1),y));
                        magMoment += (double) S(x,y);
                }
        }
}

void Metropolis(int L, int mcc, double T, map<double, double> acceptE, int numprocs, int my_rank, int myloop_end, int myloop_begin, string filename){
        random_device rd;
        mt19937_64 generator(rd());
        uniform_real_distribution<float> RNG(0.0,1.0);
        uniform_int_distribution<int> RNGPos(0,L-1);

        double energy = 0.0;  // energy
        double magMoment = 0.0; // magnetic moment (initial value, of ground state = n²)

        vec ExpectVals = zeros<vec>(5);
        vec TotalExpectVals = zeros<vec>(5);
        mat S = zeros<mat>(L,L);

        // initializing lattice
        initialize(energy, magMoment, L, S, 0, generator);

        if(my_rank == 0) outfile.open(filename, fstream::app);

        int accepts = 0;
        int cutoff = mcc*0.05/numprocs; // 5% cutoff on each thread

        for(int m = myloop_begin; m <= myloop_end; m++) {
                for (int x = 0; x < pow(L,2); x++) {
                        int xr = RNGPos(generator);         // indices for random element
                        int yr = RNGPos(generator);

                        int dE = 2.0*S(xr,yr)*(S(xr,PB(yr,L,1)) + S(xr,PB(yr,L,-1))
                                               + S(PB(xr,L,1),yr) + S(PB(xr,L,-1),yr));
                        if (RNG(generator) <= acceptE.find(dE)->second) {
                                accepts += 1;
                                S(xr,yr) *= -1.0;         // flipping spin
                                magMoment += (double) 2*S(xr,yr);
                                energy += (double) dE;
                        }
                }
                if(m >= myloop_begin+cutoff) AddExpectValsToVec(ExpectVals, energy, magMoment);

                // check if we are running 1 node. If we are, write values to file for plotting
                if(numprocs == 1) WriteToFile(energy, magMoment, mcc, T, accepts, L);

        }  //mc end
        for(int i = 0; i < 5; i++) {
                // Summing the expectvals
                MPI_Reduce(&ExpectVals(i), &TotalExpectVals(i), 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        }

        if(my_rank == 0) {
                PrintExpectVals(TotalExpectVals, mcc - cutoff*numprocs, T);
                WriteExpectValsToFile(TotalExpectVals, mcc - cutoff*numprocs, T, L);
        }
        if(my_rank == 0) outfile.close();
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

void WriteToFile(double energy, double magMoment, int mcc, double T, int accepts, int L){
        double prSpin = 1.0/pow(L,2);

        double E = energy;
        double E2 = energy*energy;
        double M = magMoment;
        double M2 = magMoment*magMoment;
        double Mabs = fabs(magMoment);

        double C_V = (E2 - E)/pow(T,2);
        double chi = (M2 - M)/T;

        outfile << E/pow(L,2) << " ";
        outfile << M/pow(L,2) << " ";
        outfile << C_V/pow(L,2) << " ";
        outfile << chi/pow(L,2) << " ";
        outfile << Mabs/pow(L,2) << " ";
        outfile << accepts << " ";
        outfile << pow(E2 - pow(E,2),2) << endl;
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

void WriteExpectValsToFile(vec TotalExpectVals, int mcc, double T, int L){
        // not using M, only |M|
        TotalExpectVals *= 1.0/mcc;
        outfile << TotalExpectVals(0)/pow(L,2) << " ";
        outfile << TotalExpectVals(4)/pow(L,2) << " ";
        outfile << (TotalExpectVals(1) - pow(TotalExpectVals(0),2))/pow(T*L,2) << " ";
        outfile << (TotalExpectVals(3) - pow(TotalExpectVals(4),2))/(T*pow(L,2)) << " ";
        outfile << T << endl;
}














//
