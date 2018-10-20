#ifndef euler_H
#define euler_H

#include "planet.h"
#include <armadillo>
#include <iostream>

using namespace arma;

class euler {

public:

double h, G;
vec accel;
vec *vel, *pos;

// initializers
euler(double G, double h, planet &otherPlanet, planet &currentPlanet);

};

#endif /*euler_H*/
