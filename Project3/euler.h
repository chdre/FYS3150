#ifndef euler_H
#define euler_H

#include "planet.h"
#include <armadillo>
#include <iostream>

using namespace std;
using namespace arma;

class euler {

public:

double h;   // step length
double accel, G;
vec vel, pos;

//
euler(double accel, double G, planet &otherPlanet);

};

#endif /*euler_H*/
