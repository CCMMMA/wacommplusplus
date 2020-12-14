//
// Created by Raffaele Montella on 14/12/20.
//

#include "WacommPlusPlus.hpp"


WacommPlusPlus::~WacommPlusPlus() = default;

WacommPlusPlus::WacommPlusPlus(std::shared_ptr<Config> config): config(config) {
    logger = log4cplus::Logger::getInstance(LOG4CPLUS_TEXT("WaComM"));

    if (config->UseRestart()) {
        particles = std::make_shared<Particles>(config->RestartFile());
    } else {
        particles = std::make_shared<Particles>();
    }
    if (config->UseSources()) {
        sources = std::make_shared<Sources>(config->ConfigFile());
    } else {
        sources =  std::make_shared<Sources>();
    }

}



void WacommPlusPlus::run() {
    LOG4CPLUS_INFO(logger,"External loop...");
    int idx=0;
    for (auto & ncInput : config->NcInputs()) {
        if (idx>=config->StartTimeIndex() && idx<config->StartTimeIndex()+config->NumberOfInputs()) {
            LOG4CPLUS_INFO(logger, "Input from Ocean Model: " << ncInput);
            ROMSAdapter romsAdapter(ncInput);
            romsAdapter.process();
            shared_ptr<OceanModelAdapter> oceanModelAdapter;
            oceanModelAdapter = make_shared<OceanModelAdapter>(romsAdapter);
            Wacomm wacomm(config, oceanModelAdapter, sources, particles);
            wacomm.run();
        }
    }
}
