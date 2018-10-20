#include "planet.h"
#include "solver.h"
#include <cmath>

solver::solver(double G, double h, planet &currentPlanet, planet &otherPlanet) {

}

void solver::euler(double G, double h, planet &currentPlanet, planet &otherPlanet) {
        /* calculating over the adresse of vel og pos, so that the values in the class
           "planet" are changed. This is again used when calculating the acceleration.
           Extract values in main to plot. */
        vel = &(currentPlanet.velocity);
        pos = &(currentPlanet.position);
        accel = currentPlanet.acceleration(otherPlanet, G);
        //cout << accel << endl;

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
