#include <iostream>
#include <filesystem>
#include <stdlib.h> /* getenv */
#include <string>

#ifdef USE_MPI
#include <mpi.h>
#endif

#include "CommandLine.hpp"
#include "Utils.hpp"

using namespace std;
namespace fs = std::filesystem;

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
#include "ROMSAdapter.hpp"

log4cplus::Logger logger;

int main(int argc, char **argv) {

#ifdef USE_MPI
    // Initialize MPI
    MPI_Init(&argc, &argv);
#endif

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
    log4cplus::LogLevel logLevel = log4cplus::INFO_LOG_LEVEL;
    std::string logLevelString = Utils::getEnvVar("WACOMM_LOGLEVEL");
    if (!logLevelString.empty()) {
        logLevel = std::stoi(logLevelString);
    }
    logger.setLogLevel(logLevel);

    LOG4CPLUS_INFO(logger, "WaComM - C++ Version");

#ifdef USE_MPI
    LOG4CPLUS_INFO(logger, "Parallel: Distributed Memory");
#endif
#ifdef USE_OMP
    LOG4CPLUS_INFO(logger, "Parallel: Shared Memory");
#endif

#ifdef USE_OPENACC
    LOG4CPLUS_INFO(logger, "Acceleration: OpenAcc");
#elif USE_OPENCL
    LOG4CPLUS_INFO(logger, "Acceleration: OpenCL");
#elif USE_CUDA
    LOG4CPLUS_INFO(logger, "Acceleration: CUDA");
#else
    LOG4CPLUS_INFO(logger, "Acceleration: None");
#endif

    LOG4CPLUS_INFO(logger,"Current directory: " << fs::current_path());

    // This variables can be set via the command line.
    std::string configFile = "namelist.wacomm";
    bool        oPrintHelp = false;

    // First configure all possible command line options.
    CommandLine args("WaComM++");
    args.addArgument({"-c", "--config"},   &configFile,   "Fortran style namelist or JSON configuration file");
    args.addArgument({"-h", "--help"},     &oPrintHelp,
      "Print this help. This help message is actually so long "
      "that it requires a line break!");

    // Then do the actual parsing.
    try {
        args.parse(argc, argv);
    } catch (std::runtime_error const& e) {
        std::cout << e.what() << std::endl;
        return -1;
    }

    // When oPrintHelp was set to true, we print a help message and exit.
    if (oPrintHelp) {
        args.printHelp();
        return 0;
    }

    LOG4CPLUS_INFO(logger,"Configuration: " << configFile);

    auto config = std::make_shared<Config>(configFile);
    config->NumberOfInputs(1);
    config->UseSources(false);
    config->Dry(false);

    WacommPlusPlus wacommPlusPlus(config);
    wacommPlusPlus.run();

#ifdef USE_MPI
    // Finalize MPI
    MPI_Finalize();
#endif
    return 0;

}
