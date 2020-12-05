//
// Created by Raffaele Montella on 12/5/20.
//

#ifndef WACOMMPLUSPLUS_PARTICLE_HPP
#define WACOMMPLUSPLUS_PARTICLE_HPP

#include <fstream>
#include <iostream>
#include <netcdf>
#include <math.h>
#include <stdlib.h> /* getenv */
#include <sstream>
#include <string>

// log4cplus - https://github.com/log4cplus/log4cplus
#include "log4cplus/configurator.h"
#include "log4cplus/logger.h"
#include "log4cplus/loggingmacros.h"

using namespace std;

class Particle {
    public:
        Particle(double x, double y, double z, double health, double tpart);
        ~Particle();

    private:
        log4cplus::Logger logger;

        double _x;
        double _y;
        double _z;
        double _health;
        double _tpart;
        double _pstatus;
};


#endif //WACOMMPLUSPLUS_PARTICLE_HPP
