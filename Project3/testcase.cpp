#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "planet.h"
#include "solver.h"


TEST_CASE("Total energy of system"){

        double M_e, M_sun, M_sunval, G, Time, timestep;
        double vfac, v0_e, t_0, k1, p1, k2, p2, energy1, energy2;
        int n;

        Time = 5.0; // time [years]
        n = 1e2; // steps
        timestep = Time/n;

        // constants
        G = 4.0*pow(M_PI,2); // gravitational constant [AU³/yr²]

        // masses
        M_sunval = 2e30;    // mass of sun
        M_sun = 1.0;        // relative mass of sun
        M_e = 6.0e24/M_sunval; // relative mass of earth

        planet Sun(M_sun, 0, 0, 0, 0, 0, 0);
        planet Earth(M_e, 1, 0, 0, 0, 2.0*M_PI, 0);


        k1 = Earth.kineticEnergy();
        p1 = Earth.potentialEnergy(G, Sun);


        energy1 = k1 + p1;
        for(int i=0; i <= n; i++) {
                solver earth(G, timestep, Earth, Sun);
                earth.VelocityVerlet(G, timestep, Earth, Sun);
                k2 = Earth.kineticEnergy();
                p2 = Earth.potentialEnergy(G, Sun);
        }


        energy2 = k2 + p2;
        cout << energy1 << " " << energy2 << endl;

        REQUIRE(energy1 == Approx(energy2));
}
