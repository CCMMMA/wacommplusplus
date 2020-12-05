//
// Created by Raffaele Montella on 12/5/20.
//

#include "Particle.hpp"

Particle::Particle(double x, double y, double z, double health, double tpart) {
    _x=x;
    _y=y;
    _z=z;
    _health=health;
    _tpart=tpart;
}

Particle::~Particle() {
}


