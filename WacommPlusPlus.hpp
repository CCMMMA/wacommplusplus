//
// Created by Raffaele Montella on 14/12/20.
//

#ifndef WACOMMPLUSPLUS_WACOMMPLUSPLUS_HPP
#define WACOMMPLUSPLUS_WACOMMPLUSPLUS_HPP

#include "OceanModelAdapters/ROMSAdapter.hpp"
#include "Wacomm.hpp"

class WacommPlusPlus {
public:
    WacommPlusPlus(std::shared_ptr<Config> config);

    ~WacommPlusPlus();

    void run();

private:
    log4cplus::Logger logger;

    std::shared_ptr<Config> config;
    std::shared_ptr<Sources> sources;
    std::shared_ptr<Particles> particles;
};


#endif //WACOMMPLUSPLUS_WACOMMPLUSPLUS_HPP
