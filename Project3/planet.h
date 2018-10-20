#ifndef planet_H
#define planet_H

#include <armadillo>
#include <iostream>

using namespace std;
using namespace arma;

class planet {
public:

double mass;
vec position;
vec velocity;

// initializers
planet();
planet(double M, double x, double y, double z, double vx, double vy, double vz);

// functions
double distance(planet otherPlanet);
double acceleration(planet otherPlanet, double G);
};

#endif /*planet_H*/
