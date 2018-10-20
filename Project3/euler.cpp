#include "planet.h"
#include "euler.h"

#include <cmath>

euler::euler(double G, planet otherPlanet, double h) {
        /* calculating over the adresse of vel og pos, so that the values in the class
           "planet" are changed. This is again used when calculating the acceleration.
           Extract values in main to plot. */
        vel = &(otherPlanet.velocity);
        pos = &(otherPlanet.position);
        accel = otherPlanet.acceleration(otherPlanet, G);
        cout << accel << endl;

        (*vel) += accel*h;
        (*pos) += (*vel)*h;
};
