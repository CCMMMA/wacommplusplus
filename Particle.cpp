//
// Created by Raffaele Montella on 12/5/20.
//

#include "Particle.hpp"

Particle::Particle(double j, double k, double i,
                   double health, double tpart):
                   j(j), k(k), i(i),
                   health(health),tpart(tpart)
{
}

Particle::~Particle() {
}

bool Particle::isAlive() {
    return health==1;
}



void Particle::move(int i, Array::Array2<double> mask, Array::Array3<double> u, Array::Array3<double> v,
                    Array::Array3<double> w, Array::Array3<double> akt) {

}

double Particle::K() { return k; }
double Particle::J() { return j; }
double Particle::I() { return i; }


