#include "planet.h"
#include "solver.h"
#include <cmath>

solver::solver(double G, double h, planet &currentPlanet, planet &otherPlanet) {
        // empty solver object for use with methods
}

solver::solver(double G, double h, planet &currentPlanet, planet &otherPlanet1, planet &otherPlanet2) {
        // empty solver object for use with methods of three objects
}

void solver::euler(double G, double h, planet &currentPlanet, planet &otherPlanet) {
        /* calculating over the adresse of vel og pos, so that the values in the class
           "planet" are changed. This is again used when calculating the acceleration.
           Extract values in main to plot. */
        vel = &(currentPlanet.velocity);
        pos = &(currentPlanet.position);
        accel = currentPlanet.acceleration(otherPlanet, G);
        
        // euler-cromer => move pos below vel
        (*pos) += (*vel)*h;
        (*vel) += accel*h;

};

void solver::VelocityVerlet(double G, double h, planet &currentPlanet, planet &otherPlanet) {
        /* see disciption above, applied to Velocity Verlet */
        pos = &(currentPlanet.position);
        vel = &(currentPlanet.velocity);
        accel = currentPlanet.acceleration(otherPlanet, G);

        (*pos) += (*vel)*h + h*h/2.0*accel;
        accel_new = currentPlanet.acceleration(otherPlanet, G);
        (*vel) += h/2.0*(accel_new + accel);
};

void solver::VelocityVerletAlt(double G, double h, double beta, planet &currentPlanet, planet &otherPlanet) {
        /* Velocity Verlet with alternative force -> 1/r^{beta} */
        pos = &(currentPlanet.position);
        vel = &(currentPlanet.velocity);
        accel = currentPlanet.accelerationAlt(otherPlanet, G, beta);

        (*pos) += (*vel)*h + h*h/2.0*accel;
        accel_new = currentPlanet.accelerationAlt(otherPlanet, G, beta);
        (*vel) += h/2.0*(accel_new + accel);
};

void solver::energy(double G, planet &currentPlanet, planet &otherPlanet) {
        // for solving the energy of planet relative to other planet
        currentPlanet.kinetic = currentPlanet.kineticEnergy();
        currentPlanet.potential = currentPlanet.potentialEnergy(G, otherPlanet);
};

void solver::VelocityVerletSystem(double G, double h, planet &currentPlanet, planet &otherPlanet1, planet &otherPlanet2) {
        // Velocity Verlet of a solar system
        pos = &(currentPlanet.position);
        vel = &(currentPlanet.velocity);
        accel = currentPlanet.newton(G, otherPlanet1, otherPlanet2);

        (*pos) += (*vel)*h + h*h/2.0*accel;
        accel_new = currentPlanet.newton(G, otherPlanet1, otherPlanet2);
        (*vel) += h/2.0*(accel_new + accel);
};

void solver::VelocityVerletEinstein(double G, double c, double h, planet &currentPlanet, planet &otherPlanet) {
        // Velocity Verlet for a relativistic force term
        pos = &(currentPlanet.position);
        vel = &(currentPlanet.velocity);
        accel = currentPlanet.einstein(G, c, otherPlanet);

        (*pos) += (*vel)*h + h*h/2.0*accel;
        accel_new = currentPlanet.einstein(G, c, otherPlanet);
        (*vel) += h/2.0*(accel_new + accel);
}
