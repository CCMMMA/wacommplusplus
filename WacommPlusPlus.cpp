//
// Created by Raffaele Montella on 14/12/20.
//

#include "WacommPlusPlus.hpp"

#ifdef USE_MPI
#include <mpi.h>
#endif

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

    int world_size=1, world_rank=0;

#ifdef USE_MPI
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
#endif

    int idx=0;
    for (auto & ncInput : config->NcInputs()) {
        if (idx>=config->StartTimeIndex() && idx<config->StartTimeIndex()+config->NumberOfInputs()) {
            LOG4CPLUS_INFO(logger, "Input from Ocean Model: " << ncInput);
            ROMSAdapter romsAdapter(ncInput);
            romsAdapter.process();
            shared_ptr<OceanModelAdapter> oceanModelAdapter;
            oceanModelAdapter = make_shared<OceanModelAdapter>(romsAdapter);

            if (world_rank==0) {
                string inputFilename = "input.nc";
                oceanModelAdapter->saveAsNetCDF(inputFilename);
            }
            Wacomm wacomm(config, oceanModelAdapter, sources, particles);
            wacomm.run();
        }
        idx++;
    }
}
