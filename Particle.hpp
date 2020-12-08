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
#include "Array.h"
#include "Config.hpp"

using namespace std;



class Particle {
    public:
        Particle(double k, double j, double i, double health, double tpart);
        ~Particle();

    bool isAlive();

    void move(Config &config,
              int ocean_time_idx,
              Array::Array1<double> &depth, Array::Array2<double> &zeta,
              Array::Array2<double> &lon, Array::Array2<double> &lat,
              Array::Array2<double> &mask,
              Array::Array4<double> &u, Array::Array4<double> &v, Array::Array4<double> &w,
              Array::Array4<double> &akt);

    double K();
    double J();
    double I();

    double gen();

private:
        log4cplus::Logger logger;

        double k;
        double j;
        double i;
        double health;
        double tpart;
        double pstatus;


};


#endif //WACOMMPLUSPLUS_PARTICLE_HPP
