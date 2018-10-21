#include "planet.h"
#include "solver.h"
#include <cmath>

solver::solver(double G, double h, planet &currentPlanet, planet &otherPlanet) {
        // empty solver object for use with methods
}

void solver::euler(double G, double h, planet &currentPlanet, planet &otherPlanet) {
        /* calculating over the adresse of vel og pos, so that the values in the class
           "planet" are changed. This is again used when calculating the acceleration.
           Extract values in main to plot. */
        vel = &(currentPlanet.velocity);
        pos = &(currentPlanet.position);
        accel = currentPlanet.acceleration(otherPlanet, G);

        // euler-cromer. Change position of vel and pos below to use euler
        (*vel) += accel*h;
        (*pos) += (*vel)*h;
};

void solver::VelocityVerlet(double G, double h, planet &currentPlanet, planet &otherPlanet) {
        pos = &(currentPlanet.position);
        vel = &(currentPlanet.velocity);
        accel = currentPlanet.acceleration(otherPlanet, G);

        (*pos) += (*vel)*h + h*h/2.0*accel;
        accel_new = currentPlanet.acceleration(otherPlanet, G);
        (*vel) += h/2.0*(accel_new + accel);
};

void solver::energy(double G, planet &currentPlanet, planet &otherPlanet) {
        currentPlanet.kinetic = currentPlanet.kineticEnergy();
        currentPlanet.potential = currentPlanet.potentialEnergy(G, otherPlanet);
};

void solver::VelocityVerletSystem(double G, double h, planet &currentPlanet, planet &otherPlanet1, planet &otherPlanet2) {
        pos = &(currentPlanet.position);
        vel = &(currentPlanet.velocity);
        accel = currentPlanet.newton(G, otherPlanet1, otherPlanet2);

        (*pos) += (*vel)*h + h*h/2.0*accel;
        accel_new = currentPlanet.newton(G, otherPlanet1, otherPlanet2);
        (*vel) += h/2.0*(accel_new + accel);
};
