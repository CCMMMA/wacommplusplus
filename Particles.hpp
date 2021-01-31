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

class Particles : private vector<Particle> {

public:
    Particles();
    Particles(string fileName);
    ~Particles();

    using vector::push_back;
    using vector::operator[];
    using vector::size;
    using vector::at;
    using vector::empty;
    using vector::begin;
    using vector::end;
    using vector::erase;

    void loadFromTxt(const string& fileName);

    void saveAsTxt(const string& fileName);
    void saveAsJson(const string &fileName, double particleTime, std::shared_ptr<OceanModelAdapter> oceanModelAdapter);
    void saveAsNetCDF(const string &fileName, double particleTime, std::shared_ptr<OceanModelAdapter> oceanModelAdapter);

private:
    log4cplus::Logger logger;
};


#endif //WACOMMPLUSPLUS_PARTICLES_HPP
