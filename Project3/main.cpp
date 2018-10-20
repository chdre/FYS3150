#include <iostream>
#include <cmath>
#include <armadillo>
#include "planet.h"
#include "euler.h"

using namespace std;


int main(){
        double M_e, G;

        double Time = 1.0;  // years
        int n = 10;            // steps

        double timestep = Time/n;

        G = 6.67e-11;                   // gravitational constant
        M_e = 6.0e24/2e30;

        planet Sun();
        planet Earth(M_e, 1, 0, 0, 0, 2.0*M_PI, 0);

        for(int i=0; i < n; i++) {
                euler(G, Earth, Sun, timestep);
                cout << Earth.position << " ";
                cout << Earth.velocity << " ";
                cout << i << endl;
        }

};
