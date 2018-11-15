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
void Metropolis(int n, int mcs, double T, map<double, double> acceptE, int numprocs, int my_rank, int myloop_end, int myloop_begin, string filename);
map<double, double> energies(double T);
void WriteToFile(double energy, double magMoment, int mcc, double T, int accepts);
void AddExpectValsToVec(vec &ExpectVals, double energy, double magMoment);
void PrintExpectVals(vec TotalExpectVals, int mcc, double T);
void WriteExpectValsToFile(vec TotalExpectVals, int mcc, double T, int n);


int main(int argc, char *argv[]){
        int n, mcc, numprocs, my_rank;
        //double T;
        string filename;

        //T = atoi(argv[4]);    // temperature [kT/J]
        n = atoi(argv[1]);    // number of spins
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
                string temp = to_string(T);
                string temp_filename = filename;
                //temp_filename.append(temp); // unlock to make one file pr T
                temp_filename.append(".txt");
                Metropolis(n, mcc, T, energies(T), numprocs, my_rank, myloop_end, myloop_begin, temp_filename);
        }

        TimeEnd = MPI_Wtime();
        if(my_rank == 0) {
                cout << "Time = " << TimeEnd - TimeStart << " on number of threads: " << numprocs << endl;
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

void Metropolis(int n, int mcc, double T, map<double, double> acceptE, int numprocs, int my_rank, int myloop_end, int myloop_begin, string filename){
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

        // initializing lattice
        initialize(energy, magMoment, n, S, 0, generator);

        if(my_rank == 0) outfile.open(filename, fstream::app);

        int accepts = 0;

        for(int m = myloop_begin; m <= myloop_end; m++) {
                for (int x = 0; x < n; x++) {
                        for (int y = 0; y < n; y++) {
                                int xr = RNGPos(generator); // indices for random element int) (RNG(generator)*(double) n); //
                                int yr = RNGPos(generator); // (int) (RNG(generator)*(double) n); // RNGSpin(generator);//;

                                int dE = 2.0*S(xr,yr)*(S(xr,PB(yr,n,1)) + S(xr,PB(yr,n,-1))
                                                       + S(PB(xr,n,1),yr) + S(PB(xr,n,-1),yr));
                                if (RNG(generator) <= acceptE.find(dE)->second) {
                                        accepts += 1;
                                        S(xr,yr) *= -1.0; // flipping spin
                                        magMoment += (double) 2*S(xr,yr);
                                        energy += (double) dE;
                                }
                        }
                }
                AddExpectValsToVec(ExpectVals, energy, magMoment);

                // check if we are running 1 node. If we are, write values to file for plotting
                //if(numprocs == 1 && m > mcc*0.05) WriteToFile(energy, magMoment, mcc, T, accepts);

        }  //mc end
        for(int i = 0; i < 5; i++) {
                // Summing the expectvals
                MPI_Reduce(&ExpectVals(i), &TotalExpectVals(i), 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        }

        if(my_rank == 0) {
                PrintExpectVals(TotalExpectVals, mcc, T);
                WriteExpectValsToFile(TotalExpectVals, mcc, T, n);
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

void WriteToFile(double energy, double magMoment, int mcc, double T, int accepts){
        double E = energy;
        double E2 = energy*energy;
        double M = magMoment;
        double M2 = magMoment*magMoment;
        double Mabs = fabs(magMoment);

        double C_V = (E2 - E)/pow(T,2);
        double chi = (M2 - M)/T;

        outfile << E << " ";
        outfile << M << " ";
        outfile << C_V << " ";
        outfile << chi << " ";
        outfile << Mabs << " ";
        outfile << accepts << endl;
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

void WriteExpectValsToFile(vec TotalExpectVals, int mcc, double T, int n){
        TotalExpectVals = TotalExpectVals/(mcc*pow(n,2));
        outfile << TotalExpectVals(0) << " ";
        outfile << TotalExpectVals(4) << " ";
        outfile << (TotalExpectVals(1) - pow(TotalExpectVals(0),2))/(pow(T,2)) << " ";
        outfile << (TotalExpectVals(3) - pow(TotalExpectVals(2),2))/T << " ";
        outfile << T << endl;
}














//
