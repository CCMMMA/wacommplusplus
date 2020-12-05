//
// Created by Raffaele Montella on 12/5/20.
//

#ifndef WACOMMPLUSPLUS_PARTICLES_HPP
#define WACOMMPLUSPLUS_PARTICLES_HPP

// log4cplus - https://github.com/log4cplus/log4cplus
#include "log4cplus/configurator.h"
#include "log4cplus/logger.h"
#include "log4cplus/loggingmacros.h"

#include "Particle.hpp"

using namespace std;

class Particles : public vector<Particle> {

public:
    Particles();
    Particles(string fileName);
    ~Particles();

    void save(string fileName);

private:
    log4cplus::Logger logger;

};


#endif //WACOMMPLUSPLUS_PARTICLES_HPP
