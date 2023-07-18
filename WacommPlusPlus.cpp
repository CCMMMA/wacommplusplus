//
// Created by Raffaele Montella on 14/12/20.
//

#include "WacommPlusPlus.hpp"
#include "JulianDate.hpp"
#include "OceanModelAdapters/WacommAdapter.hpp"
#include "OceanModelAdapters/ROMSAdapter.hpp"

#if defined(USE_MPI) || defined(USE_EMPI)
#define OMPI_SKIP_MPICXX
#include <mpi.h>
#endif

#ifdef USE_EMPI
extern "C" {
#include <empi.h>
}
#endif

WacommPlusPlus::~WacommPlusPlus() = default;

WacommPlusPlus::WacommPlusPlus(std::shared_ptr<Config> config): config(config) {
    logger = log4cplus::Logger::getInstance(LOG4CPLUS_TEXT("WaComM"));

    // Create the particles
    particles = std::make_shared<Particles>();

    // Create the sources
    sources =  std::make_shared<Sources>();
}



void WacommPlusPlus::run() {
    LOG4CPLUS_DEBUG(logger, "External loop...");

    int world_size = 1, world_rank = 0;
    char mpi_name[128];
    int len, nParticles = 0;

    char bin[1024];
    sprintf(bin,"%s","wacomm");

    int idx = 0;
    double time_average = 0.0;
    double part_average = 0.0;
    double cuda_average = 0.0;

#ifdef USE_MPI
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
#endif

#ifdef USE_EMPI
    MPI_Comm_size(ADM_COMM_WORLD, &world_size);
    MPI_Comm_rank(ADM_COMM_WORLD, &world_rank);

    MPI_Get_processor_name (mpi_name, &len);
    
    if (world_rank == 0) {
        nParticles = -1;
    }

    // get process type
    int proctype;
    ADM_GetSysAttributesInt ("ADM_GLOBAL_PROCESS_TYPE", &proctype);

    // if process is native
    if (proctype == ADM_NATIVE) {
        printf ("Rank(%d/%d): Process native\n", world_rank, world_size);
    // if process is spawned
    } else {
        printf ("Rank(%d/%d): Process spawned\n", world_rank, world_size);
    }

    ADM_GetSysAttributesInt ("ADM_GLOBAL_ITERATION", &idx);
    ADM_GetSysAttributesInt ("ADM_GLOBAL_PARTICLES", &nParticles);

    /* starting monitoring service */
    ADM_MonitoringService (ADM_SERVICE_START);
#endif

    for (auto &ncInput : config->NcInputs()) {

        if (world_rank == 0) {
            LOG4CPLUS_INFO(logger, world_rank << ": Input from Ocean Model: " << ncInput);
        }

        shared_ptr<OceanModelAdapter> oceanModelAdapter;
        if (config->OceanModel() == "ROMS") {
            oceanModelAdapter = make_shared<ROMSAdapter>(ncInput);

        } else {
            oceanModelAdapter = make_shared<WacommAdapter>(ncInput);
        }
        oceanModelAdapter->process();

        Calendar cal;

        // Time in "seconds since 1968-05-23 00:00:00"
        double modJulian=oceanModelAdapter->OceanTime()(0);

        // Convert time in days based
        modJulian=modJulian/86400;

        JulianDate::fromModJulian(modJulian, cal);

        // Check if the rank is 0
        if (world_rank == 0) {

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

            // Check if the processed input must be saved
            if (config->SaveInput()) {

                // Create the filename
                string inputFilename = config->NcInputRoot() + cal.asNCEPdate() + ".nc";

                // Show a information message
                LOG4CPLUS_INFO(logger,  "Saving processed output: " << inputFilename);

                // Save the processed input
                oceanModelAdapter->saveAsNetCDF(inputFilename);
            }

            // Check if using the restart file (only if it is the first iteration)
            if (idx==0 && config->UseRestart() && !config->RestartFile().empty()) {

                // Get the file name
                string fileName = config->RestartFile();

                // Check if the restart is a NetCDF
                if (fileName.substr(fileName.find_last_of('.') + 1) == "nc") {

                    // The restart is a NetCDF
                    particles->loadFromNetCDF(fileName);

                    // Check if the restart is a geojson
                } else if (fileName.substr(fileName.find_last_of('.') + 1) == "json") {

                    // The restart is a json
                    particles->loadFromJson(fileName);
                } else {

                    // the restart is a fortran style text file
                    particles->loadFromTxt(fileName);
                }

            }
        }

        // Check if it is a dry run
        if (!config->Dry()) {
            // Create a new Wacomm object
            Wacomm wacomm(config, oceanModelAdapter, sources, particles);

            // Run the model
            double t=0, p=0, cuda=0;

            int status = wacomm.run(t,p,cuda,nParticles,idx);
            
            // check if process ended after malleable region
            if (status == ADM_ACTIVE) {
                // updata world_rank and size
                MPI_Comm_rank(ADM_COMM_WORLD, &world_rank);
                MPI_Comm_size(ADM_COMM_WORLD, &world_size);

                printf("world_size: %d\n", world_size);
            } else {
                // end the process
                break;
            }

            time_average += t;
            part_average += p;
            cuda_average += cuda;
        }
        // Go to the next input file
        idx++;
    }

#ifdef USE_EMPI
    /* ending monitoring service */
    ADM_MonitoringService (ADM_SERVICE_STOP);
#endif

    // LOG4CPLUS_INFO(logger,  "Outer Cycle Time (Average): " << time_average / (idx-1));
    // LOG4CPLUS_INFO(logger,  "Outer Cycle Particles/sec (Average): " << part_average / (idx-1));
    // LOG4CPLUS_INFO(logger,  "Inner Cycle (Average) sec: " << cuda_average / (idx-1));
}
