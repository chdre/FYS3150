#ifndef solver_H
#define solver_H

#include "planet.h"
#include <armadillo>
#include <iostream>

using namespace arma;

class solver {

public:

double h, G;
vec accel, accel_new;
vec *vel, *pos;

//contructor
solver(double G, double h, planet &currentPlanet, planet &otherPlanet);
solver(double G, double h, planet &currentPlanet, planet &otherPlanet1, planet &otherPlanet2);
// function
void euler(double G, double h, planet &otherPlanet, planet &currentPlanet);
void VelocityVerlet(double G, double h, planet &currentPlanet, planet &otherPlanet);
void energy(double G, planet &currentPlanet, planet &otherPlanet);
void VelocityVerletSystem(double G, double h, planet &currentPlanet, planet &otherPlanet1, planet &otherPlanet2);
};

#endif /*solver_H*/
