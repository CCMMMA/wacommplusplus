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

using namespace std;

class Particle {
    public:
        Particle(double k, double j, double i, double health, double tpart);
        ~Particle();

    bool isAlive();

    void move(int i, Array::Array2<double> array2, Array::Array3<double> u, Array::Array3<double> v,
              Array::Array3<double> w, Array::Array3<double> akt);

    double K();
    double J();
    double I();

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
