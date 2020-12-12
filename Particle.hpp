//
// Created by Raffaele Montella on 12/5/20.
//

#ifndef WACOMMPLUSPLUS_PARTICLE_HPP
#define WACOMMPLUSPLUS_PARTICLE_HPP

#include <fstream>
#include <iostream>
#include <random>
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
#include "OceanModelAdapter.hpp"

using namespace std;



class Particle {
    public:
        Particle(double k, double j, double i, double health, double tpart);
        Particle(double k, double j, double i);

        ~Particle();

        bool isAlive() const;

        void move(const std::shared_ptr<Config>& config, int ocean_time_idx,
                  const std::shared_ptr<OceanModelAdapter>& oceanModelAdapter);

        double K() const;
        double J() const;
        double I() const;

        double Health() const;
        double TPart() const;

        int KasInt() const;
        int JasInt() const;
        int IasInt() const;

        static double gen();

        std::string to_string() const;

private:
        log4cplus::Logger logger;



        double k;
        double j;
        double i;
        double health;
        double tpart;
        double pstatus{};

        static double sgn(double a);
        static double mod(double a, double p);
        static double sign(double a, double b);

        int time{};
        double health0{};


};


#endif //WACOMMPLUSPLUS_PARTICLE_HPP
