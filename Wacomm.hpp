//
// Created by Raffaele Montella on 12/5/20.
//

#ifndef WACOMMPLUSPLUS_WACOMM_HPP
#define WACOMMPLUSPLUS_WACOMM_HPP

// log4cplus - https://github.com/log4cplus/log4cplus
#include <memory>
#include <cstring>

#include "log4cplus/configurator.h"
#include "log4cplus/logger.h"
#include "log4cplus/loggingmacros.h"

#include "Sources.hpp"
#include "Particles.hpp"
#include "Array.h"
#include "OceanModelAdapter.hpp"

using namespace std;
using namespace Array;

class ROMSAdapter;
class Wacomm {
    public:

        Wacomm(std::shared_ptr<Config> config,
               std::shared_ptr<OceanModelAdapter> oceanModelAdapter,
               std::shared_ptr<Sources> sources,
               std::shared_ptr<Particles> particles);
        ~Wacomm();

        void run();

    private:
        friend class Wacomm;

        log4cplus::Logger logger;

        std::shared_ptr<Config> config;
        std::shared_ptr<OceanModelAdapter> oceanModelAdapter;

        std::shared_ptr<Sources> sources;
        std::shared_ptr<Particles> particles;

        void save(Array4<float> &conc);
};




#endif //WACOMMPLUSPLUS_WACOMM_HPP
