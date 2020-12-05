#include <iostream>
#include <fstream>
#include <iostream>
#include <netcdf>
#include <math.h>
#include <stdlib.h> /* getenv */
#include <sstream>
#include <string>

using namespace std;


// log4cplus - https://github.com/log4cplus/log4cplus
#include "log4cplus/configurator.h"
#include "log4cplus/logger.h"
#include "log4cplus/loggingmacros.h"

#include "Wacomm.hpp"
#include "Particles.hpp"

log4cplus::Logger logger;

std::string getEnvVar(std::string const &key) {
    char *val = getenv(key.c_str());
    return val == NULL ? std::string("") : std::string(val);
}

int main() {
    // Logger configuration
    log4cplus::BasicConfigurator config;
    config.configure();
    logger = log4cplus::Logger::getInstance(LOG4CPLUS_TEXT("WaComM"));

    // Set the logging level
    log4cplus::LogLevel logLevel = log4cplus::INFO_LOG_LEVEL;
    std::string logLevelString = getEnvVar("WACOMM_LOGLEVEL");
    if (logLevelString != "") {
        logLevel = std::stoi(logLevelString);
    }
    logger.setLogLevel(logLevel);

    LOG4CPLUS_INFO(logger, "🛈  - WaComM - C++ Version");

    string restartFileName="/projects/wacomm/data/input/WACOMM_rst_20201130Z00.txt";
    string fileName="/projects/wacomm/data/gnu-8.3.1/input/rms3_d03_20201130Z0000.nc";


    //Particles particles(restartFileName);
    Particles particles;
    Wacomm wacomm(fileName,particles);

    // ROMS2Wacomm roms2Wacomm(fileName, ocean_time, ucomp, vcomp, wcomp, aktcomp);
    // Array4<double> conc=wacomm.run(ocean_time, ucomp, vcomp, wcomp, aktcomp);
    // Results results(conc);
    // results.saveAsNetCDF(concFileName);
    // Restart restart(particles);
    // restart.saveAsTextfile(restartFileName);
    return 0;

}