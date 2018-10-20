#ifndef euler_H
#define euler_H

#include "planet.h"
#include <armadillo>
#include <iostream>

using namespace arma;

class euler {

public:

double h;   // step length
double accel, G;
vec *vel, *pos;

// initializers
euler(double G, planet otherPlanet, double h);

};

#endif /*euler_H*/
