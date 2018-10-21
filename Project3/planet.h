#ifndef planet_H
#define planet_H

#include <armadillo>
#include <iostream>

using namespace arma;

class planet {
public:

double mass, potential, kinetic;
vec position;
vec velocity;

// initializers
planet();
planet(double M, double x, double y, double z, double vx, double vy, double vz);

// functions
double distance(planet otherPlanet);
vec acceleration(planet otherPlanet, double G);
double kineticEnergy();
double potentialEnergy(double G, planet otherPlanet);

};

#endif /*planet_H*/
