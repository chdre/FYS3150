#include "euler.h"
#include "planet.h"
#include <cmath>

euler::euler(double accel, double G, planet &otherPlanet) {
        /* calculating over the adresse of vel og pos, so that the values in the class
           "planet" are changed. This is again used when calculating the acceleration.
           Extract values in main to plot. */
        vel = otherPlanet.velocity;
        pos = otherPlanet.position;
        accel = otherPlanet.acceleration(otherPlanet, G);

        vel += accel*h;
        pos += vel*h;
};
