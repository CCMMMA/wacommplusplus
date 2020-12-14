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

#include "WacommPlusPlus.hpp"
#include "Particles.hpp"
#include "Sources.hpp"
#include "ROMSAdapter.hpp"

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

    return 0;

}