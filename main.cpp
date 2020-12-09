#include <iostream>
#include <filesystem>
#include <stdlib.h> /* getenv */
#include <string>

#include "CommandLine.hpp"

using namespace std;
namespace fs = std::filesystem;

// log4cplus - https://github.com/log4cplus/log4cplus
#include "log4cplus/configurator.h"
#include "log4cplus/logger.h"
#include "log4cplus/loggingmacros.h"

#include "Wacomm.hpp"
#include "Particles.hpp"
#include "Sources.hpp"
#include "ROMS2Wacomm.hpp"

log4cplus::Logger logger;

std::string getEnvVar(std::string const &key) {
    char *val = getenv(key.c_str());
    return val == NULL ? std::string("") : std::string(val);
}

int main(int argc, char **argv) {
    // Logger configuration
    log4cplus::BasicConfigurator basicConfig;
    basicConfig.configure();
    logger = log4cplus::Logger::getInstance(LOG4CPLUS_TEXT("WaComM"));

    // Set the logging level
    log4cplus::LogLevel logLevel = log4cplus::INFO_LOG_LEVEL;
    std::string logLevelString = getEnvVar("WACOMM_LOGLEVEL");
    if (logLevelString != "") {
        logLevel = std::stoi(logLevelString);
    }
    logger.setLogLevel(logLevel);

    LOG4CPLUS_INFO(logger, "🛈  - WaComM - C++ Version");
    LOG4CPLUS_INFO(logger,"Current directory: " << fs::current_path());



    // This variables can be set via the command line.
    std::string configFile = "namelist.wacomm";
    std::string restartFileName = "WACOMM_rst_20201130Z00.txt";
    std::string sourcesFileName = "input/sources.txt";
    std::string netcdfFileName = "input/rms3_d03_20201130Z0000.nc";
    bool        oPrintHelp = false;

    // First configure all possible command line options.
    CommandLine args("WaComM++");
    args.addArgument({"-c", "--config"},   &configFile,   "Fortran style namelist or JSON configuration file");
    args.addArgument({"-rst", "--restart"},  &restartFileName,  "Restart file");
    args.addArgument({"-src", "--sources"},  &sourcesFileName,  "Sources file");
    args.addArgument({"-nc", "--netcdf"},  &netcdfFileName,  "NetCDF file");
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

    Config config(configFile);

    Particles particles(restartFileName);
    Sources sources;

    auto roms2Wacomm = std::make_shared<ROMS2Wacomm>(netcdfFileName);
    roms2Wacomm->process();

    Wacomm wacomm(config, roms2Wacomm, sources, particles);
    wacomm.run();

    // ROMS2Wacomm roms2Wacomm(fileName, ocean_time, ucomp, vcomp, wcomp, aktcomp);
    // Array4<double> conc=wacomm.run(ocean_time, ucomp, vcomp, wcomp, aktcomp);
    // Results results(conc);
    // results.saveAsNetCDF(concFileName);
    // Restart restart(particles);
    // restart.saveAsTextfile(restartFileName);
    return 0;

}