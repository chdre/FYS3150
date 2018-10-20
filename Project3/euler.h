#ifndef euler_H
#define euler_H

#include <armadillo>
#include <iostream>

using namespace std;
using namespace arma;

class euler {

public:

double h;   // step length
double vel, pos, accel;

//
euler(double accel, double G, planet &otherPlanet);

};

#endif /*euler_H*/
