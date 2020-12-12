//
// Created by Raffaele Montella on 08/12/20.
//

#include "Config.hpp"

Config::Config() = default;

Config::Config(string &fileName) {
    setDefault();

    log4cplus::BasicConfigurator basicConfig;
    basicConfig.configure();
    logger = log4cplus::Logger::getInstance(LOG4CPLUS_TEXT("WaComM"));

    LOG4CPLUS_INFO(logger, "Reading config file:" + fileName);

    std::ifstream infile(fileName);
}

Config::~Config() {
    setDefault();
}


void Config::setDefault() {
    // Is the model dry mode?
    dry = false;

    // Input time step in s
    deltat=3600;

    // Integration time in s
    dti=30;

    // To be clarified
    tau0=86400.0;

    // Probability of particle surviving
    survprob=1.0e-4;

    // Sedimentation velocity
    sv=0;

    // Reduction Coefficient
    crid=1;
}

void Config::fromNamelist(ifstream &ifs) {

}

void Config::formJson(ifstream &ifs) {

}

double Config::Dti() const { return dti; }

double Config::Deltat() const { return deltat; }

double Config::Survprob() const {
    return survprob;
}

double Config::Tau0() const {
    return tau0;
}

double Config::SedimentationVelocity() const {
    return sv;
}

bool Config::isDry() const {
    return dry;
}

void Config::setDry(bool value) {
    dry=value;
}

double Config::ReductionCoefficient() const {
    return crid;
}
