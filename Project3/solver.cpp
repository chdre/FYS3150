#include "planet.h"
#include "solver.h"
#include <cmath>

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
