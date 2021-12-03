#include <iostream>
#include <stdlib.h> /* getenv */
#include <string>

#ifdef USE_CUDA
#include <cuda.h>
#include <cuda_runtime_api.h>
#endif

#ifdef USE_MPI
#define OMPI_SKIP_MPICXX
#include <mpi.h>
#endif

#ifdef USE_OMP
#include <omp.h>
#endif

#include "Utils.hpp"

using namespace std;


// log4cplus - https://github.com/log4cplus/log4cplus
#include "log4cplus/configurator.h"
#include "log4cplus/logger.h"
#include "log4cplus/loggingmacros.h"
#include "log4cplus/initializer.h"
#include "log4cplus/consoleappender.h"
#include "log4cplus/layout.h"

#include "WacommPlusPlus.hpp"
#include "Particles.hpp"
#include "Sources.hpp"
#include "OceanModelAdapters/ROMSAdapter.hpp"
#include "JulianDate.hpp"

#include <iostream>


log4cplus::Logger logger;

int main(int argc, char **argv) {
    int num_gpus = 0;
    int world_size=1, world_rank=0;
    int ompMaxThreads=1;

#ifdef USE_OMP
    // The the number of threads
    ompMaxThreads=omp_get_max_threads();
#endif

#ifdef USE_MPI
    // Initialize MPI
    MPI_Init(&argc, &argv);

    // Get the number of involved processes
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the number of the current process (world_rank=0 is for the main process)
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
#endif

#ifdef USE_CUDA
    cudaGetDeviceCount(&num_gpus);
#endif

    // This variables can be set via the command line.
    std::string configFile = "namelist.wacomm";
    configFile = "wacomm.json";

    // Inizitalizer
    log4cplus::Initializer initializer;

    // Basic configuration
    log4cplus::BasicConfigurator basicConfigurator;
    basicConfigurator.configure();

    //Create an appender pointing to the console
    log4cplus::SharedAppenderPtr appender(new log4cplus::ConsoleAppender());

    //Set the layout of the appender
    appender->setName(LOG4CPLUS_TEXT("console"));

    log4cplus::tstring pattern = LOG4CPLUS_TEXT("%D{%y-%m-%d %H:%M:%S,%Q} %-5p %c");
    appender->setLayout(std::unique_ptr<log4cplus::Layout>(new log4cplus::PatternLayout(pattern)));

    logger = log4cplus::Logger::getInstance(LOG4CPLUS_TEXT("WaComM"));

    logger.addAppender(appender);

    // Set the logging level
#ifdef DEBUG
    log4cplus::LogLevel logLevel = log4cplus::DEBUG_LOG_LEVEL;
#else
    log4cplus::LogLevel logLevel = log4cplus::INFO_LOG_LEVEL;
#endif
    std::string logLevelString = Utils::getEnvVar("WACOMM_LOGLEVEL");
    if (!logLevelString.empty()) {
        logLevel = std::stoi(logLevelString);
    }
    logger.setLogLevel(logLevel);

    if (argc==2) {
        configFile=string(argv[1]);
    }

    // Check if the process is the main one
    if (world_rank == 0) {
        LOG4CPLUS_INFO(logger, "WaComM - C++ Version");

#ifdef USE_MPI
        LOG4CPLUS_INFO(logger, "Parallel: Distributed Memory");
#endif
#ifdef USE_OMP
        LOG4CPLUS_INFO(logger, "Parallel: Shared Memory");
#endif
#ifndef USE_OMP
#ifndef USE_MPI
        LOG4CPLUS_INFO(logger, "Parallel: None");
#endif
#endif
#ifdef USE_OPENACC
        LOG4CPLUS_INFO(logger, "Acceleration: OpenAcc");
#elif USE_CUDA
        LOG4CPLUS_INFO(logger, "Acceleration: CUDA " << num_gpus << " device(s)");
        for (int i=0; i<num_gpus; i++){
		cudaDeviceProp prop;
    		cudaGetDeviceProperties(&prop, i);
    		LOG4CPLUS_INFO(logger, "Device Number: " << i);
    		LOG4CPLUS_INFO(logger, "Device name: " <<  prop.name);
    		LOG4CPLUS_INFO(logger, "Memory Clock Rate (KHz): " << prop.memoryClockRate);
    		LOG4CPLUS_INFO(logger, "Memory Bus Width (bits): " << prop.memoryBusWidth);
    		LOG4CPLUS_INFO(logger, "Peak Memory Bandwidth (GB/s): " << 2.0*prop.memoryClockRate*(prop.memoryBusWidth/8)/1.0e6);
	}
#else
        LOG4CPLUS_INFO(logger, "Acceleration: None");
#endif

        LOG4CPLUS_INFO(logger, world_rank << ": Using 1/" << world_size << " processes, each on " << ompMaxThreads
                                          << " threads.");
        LOG4CPLUS_INFO(logger, "Configuration: " << configFile);
    }

    // Load the configuration file
    auto config = std::make_shared<Config>(configFile);

    // Check if it is the main process
    if (world_rank == 0) {

        // Show a message if it is a dry run
        if (config->Dry()) {
            LOG4CPLUS_INFO(logger, "*** Dry run ***");
        }
    }

    WacommPlusPlus wacommPlusPlus(config);

    // Record start time
    auto start_whole = std::chrono::high_resolution_clock::now();

    wacommPlusPlus.run();
    auto end_whole = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsedLocal = end_whole - start_whole;
    LOG4CPLUS_INFO(logger,  "Global Execution Time (sec): " << elapsedLocal.count() / 24);

#ifdef USE_MPI
    // Finalize MPI
    MPI_Finalize();
#endif
    return 0;

}
