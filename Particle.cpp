//
// Created by Raffaele Montella on 12/5/20.
//

#include "Particle.hpp"

Particle::Particle(double j, double k, double i, double health, double tpart) {
    _k=k;
    _j=j;
    _i=i;
    _health=health;
    _tpart=tpart;
}

Particle::~Particle() {
}

bool Particle::isAlive() {
    return _health==1;
}

double Particle::k() { return _k; }

void Particle::move(int i, Array::Array2<double> mask, Array::Array3<double> u, Array::Array3<double> v,
                    Array::Array3<double> w, Array::Array3<double> akt) {

}

double Particle::j() { return _j; }
double Particle::i() { return _i; }


