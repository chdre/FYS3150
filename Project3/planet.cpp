#include "planet.h"
#include <cmath>

/*
   1) Call either planet or planet(w/ arguments) to set sun or otherPlanet.
   2) calculate distance.
 */

planet::planet(){
        // sun
        mass = 1.0;
        position = {0.0, 0.0, 0.0};
        velocity = {0.0, 0.0, 0.0};
}

planet::planet(double M, double x, double y, double z, double vx, double vy, double vz)
{
        // for other planets
        mass = M;
        position = vec({x, y, z});
        velocity = vec({vx, vy, vz});
}

double planet::distance(planet otherPlanet) {
        // returning the distance between planets/sun-planet.
        double x, y, z, x1, x2, y1, y2, z1, z2;

        // position of planet (sun)
        x1 = this->position[0];
        y1 = this->position[1];
        z1 = this->position[2];

        // position of other planet
        x2 = otherPlanet.position[0];
        y2 = otherPlanet.position[1];
        z2 = otherPlanet.position[2];

        x = x1 - x2;
        y = y1 - y2;
        z = z1 - z2;
        cout << sqrt(pow(x,2) + pow(y,2) + pow(z,2)) << endl;

        return sqrt(pow(x,2) + pow(y,2) + pow(z,2));
}

double planet::acceleration(planet otherPlanet, double G) {
        // calculating the acceleration a = -GM_sun/r²
        double r = this->distance(otherPlanet);
        if (r != 0) {
                return -G*this->mass/pow(r,2);
        }
        else {
                return 0;
        }
};
