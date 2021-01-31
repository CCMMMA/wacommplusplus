//
// Created by Raffaele Montella on 14/12/20.
//

#include "WacommPlusPlus.hpp"
#include "JulianDate.hpp"
#include "OceanModelAdapters/WacommAdapter.hpp"
#include "OceanModelAdapters/ROMSAdapter.hpp"
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
    sources =  std::make_shared<Sources>();
}



void WacommPlusPlus::run() {
    LOG4CPLUS_DEBUG(logger, "External loop...");

    int world_size = 1, world_rank = 0;

#ifdef USE_MPI
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
#endif

    int idx = 0;
    for (auto &ncInput : config->NcInputs()) {

        if (world_rank == 0) {
            LOG4CPLUS_INFO(logger, world_rank << ": Input from Ocean Model: " << ncInput);
        }

        OceanModelAdapter *pOceanModelAdapter;
        if (config->OceanModel() == "ROMS") {
            ROMSAdapter romsAdapter(ncInput);
            pOceanModelAdapter=&romsAdapter;
        } else {
            WacommAdapter wacommAdapter(ncInput);
            pOceanModelAdapter=&wacommAdapter;
        }
        pOceanModelAdapter->process();
        auto oceanModelAdapter = make_shared<OceanModelAdapter>(*pOceanModelAdapter);

        cout << "WacommPlusPlus::run --oceanModelAdapter->OceanTime()(0): " << oceanModelAdapter->OceanTime()(0) << endl;
        cout << "WacommPlusPlus::run --oceanModelAdapter->H()(650,550):" << oceanModelAdapter->H()(650,550) << endl;

        Calendar cal;

        // Time in "seconds since 1968-05-23 00:00:00"
        double modJulian=oceanModelAdapter->OceanTime()(0);

        // Convert time in days based
        modJulian=modJulian/86400;

        JulianDate::fromModJulian(modJulian, cal);

        // Check if it is needed to load the sources
        if (config->UseSources() && sources->empty()) {
            string fileName = config->SourcesFile();
            if (fileName.empty()) {
                sources->loadFromNamelist(config->ConfigFile());
            } else {
                if (fileName.substr(fileName.find_last_of('.') + 1) == "json") {
                    // The configuration is a json
                    sources->loadFromJson(fileName, oceanModelAdapter);
                } else {
                    // the configuration is a fortran style namelist
                    sources->loadFromNamelist(fileName);
                }
            }
        }

        cout << "WacommPlusPlus::run ++oceanModelAdapter->OceanTime()(0): " << oceanModelAdapter->OceanTime()(0) << endl;
        cout << "WacommPlusPlus::run ++oceanModelAdapter->H()(650,550):" << oceanModelAdapter->H()(650,550) << endl;



        if (world_rank == 0) {
            if (config->SaveInput()) {
                string inputFilename = config->NcInputRoot() + cal.asNCEPdate() + ".nc";
                oceanModelAdapter->saveAsNetCDF(inputFilename);
            }
        }


        Wacomm wacomm(config, oceanModelAdapter, sources, particles);
        wacomm.run();

        idx++;
    }
}
