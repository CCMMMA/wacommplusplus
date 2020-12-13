//
// Created by Raffaele Montella on 07/12/20.
//

#include "Source.hpp"

Source::Source(int id, double k, double j, double i, int start, int end, int particlesPerHour, int mode):
    id(id),k(k),j(j),i(i),start(start),end(end),particlesPerHour(particlesPerHour),mode(mode) {

    logger = log4cplus::Logger::getInstance(LOG4CPLUS_TEXT("WaComM"));
}

Source::~Source() = default;

void Source::emit(std::shared_ptr<Particles> particles, double tpart) {
    // Check if the source is active
    if (mode > 0) {
        for (int i=0; i<particlesPerHour;i++) {
            particles->push_back(Particle(
                    k+ gen()*0.5 - 0.25,
                    j+ gen()*0.5 - 0.25,
                    i+ gen()*0.5 - 0.25, tpart));
        }
    }
}

double Source::gen()
{
    return Random::get<double>(0.0, 1.0);
}