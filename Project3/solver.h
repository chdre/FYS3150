#ifndef solver_H
#define solver_H

#include "planet.h"
#include <armadillo>
#include <iostream>

using namespace arma;

class solver {

public:

double h, G;
vec accel;
vec *vel, *pos;

// initializers
void euler(double G, double h, planet &otherPlanet, planet &currentPlanet);

};

#endif /*solver_H*/
