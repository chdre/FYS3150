#include "planet.h"
#include <cmath>

/*
   1) Call either planet or planet(w/ arguments) to set sun or otherPlanet.
   2) calculate distance.
 */

planet::planet(){
        // sun
        mass = 0.0;
        position = vec({1.0, 0.0, 0.0});
        velocity = vec({0.0, 0.0, 0.0});
        potential = 0;
        kinetic = 0;
        angularMom = 0;
}

planet::planet(double M, double x, double y, double z, double vx, double vy, double vz)
{
        // for other planets
        mass = M;
        position = vec({x, y, z});
        velocity = vec({vx, vy, vz});
        potential = 0;
        kinetic = 0;
        angularMom = 0;
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
        // cout << sqrt(pow(x,2) + pow(y,2) + pow(z,2)) << endl;

        return sqrt(pow(x,2) + pow(y,2) + pow(z,2));
};

vec planet::acceleration(planet otherPlanet, double G) {
        // calculating the acceleration a = -GM_sun/r²
        double r = this->distance(otherPlanet);
        if (r != 0) {
                double otherMass = otherPlanet.mass;
                vec currentPos = this->position;
                vec otherPos = otherPlanet.position;
                return -G*otherMass*(currentPos - otherPos)/pow(r,3);
        }
        else {
                return vec({0,0,0});
        }
};


double planet::kineticEnergy(){
        double vel2 = pow(this->velocity[0],2) + pow(this->velocity[1],2) + pow(this->velocity[2],2);
        return this->mass*vel2/2.0;
};

double planet::potentialEnergy(double G, planet otherPlanet){
        return -G*otherPlanet.mass*this->mass/this->distance(otherPlanet);
};

double planet::angularMomentum(planet otherPlanet){
        double vel = sqrt(pow(this->velocity[0],2) + pow(this->velocity[1],2) + pow(this->velocity[2],2));
        return this->mass*vel*this->distance(otherPlanet);
};

double planet::angularMomentum(planet otherPlanet1,planet otherPlanet2){
        double vel = sqrt(pow(this->velocity[0],2) + pow(this->velocity[1],2) + pow(this->velocity[2],2));
        double p1 = this->mass*vel*this->distance(otherPlanet1);
        double p2 = this->mass*vel*this->distance(otherPlanet2);
        return p1 + p2;
};


vec planet::newton(double G, planet otherPlanet1, planet otherPlanet2) {
        double r1 = this->distance(otherPlanet1);
        double r2 = this->distance(otherPlanet2);
        double currentMass = this->mass;
        vec F1, F2;
        if (r1 != 0) {
                double otherMass = otherPlanet1.mass;
                vec currentPos1 = this->position;
                vec otherPos1 = otherPlanet1.position;
                F1 = -G*otherMass*currentMass*(currentPos1 - otherPos1)/pow(r1,3);
        }
        else{
                F1 = vec({0,0,0});
        }
        if (r2 != 0) {
                double otherMass = otherPlanet2.mass;
                vec currentPos2 = this->position;
                vec otherPos2 = otherPlanet2.position;
                F2 = -G*otherMass*currentMass*(currentPos2 - otherPos2)/pow(r2,3);
        }
        else{
                F2 = vec({0,0,0});
        }
        return (F1 + F2)/currentMass;;
}

vec planet::einstein(double G, double c, planet otherPlanet) {
        double r = this->distance(otherPlanet);
        double currentMass = this->mass;
        double otherMass = otherPlanet.mass;
        vec currentPos = this->position;
        vec currentVel = this->velocity;
        vec otherPos = otherPlanet.position;
        vec l = cross((currentPos - otherPos),currentVel);

        vec F = (-G*otherMass*(currentPos-otherPos)/pow(r,3))%(1 + (3*pow(l,2)/(pow(r,2)*pow(c,2))));
        return F/currentMass;
}


//
