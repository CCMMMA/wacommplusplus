//
// Created by Raffaele Montella on 08/12/20.
//

#include "Config.hpp"

Config::Config() {

}

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
    deltat=3600;
    dti=30;
    tau0=86400.0;
    survprob=1.0e-4;
}

void Config::fromNamelist(ifstream &ifs) {

}

void Config::formJson(ifstream &ifs) {

}

double Config::Dti() { return dti; }

double Config::Deltat() { return deltat; }

double Config::Survprob() {
    return survprob;
}

double Config::Tau0() {
    return tau0;
}
