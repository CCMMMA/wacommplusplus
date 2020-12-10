//
// Created by Raffaele Montella on 08/12/20.
//

#ifndef WACOMMPLUSPLUS_CONFIG_HPP
#define WACOMMPLUSPLUS_CONFIG_HPP

#include <string>
#include <fstream>

using namespace std;

// log4cplus - https://github.com/log4cplus/log4cplus
#include "log4cplus/configurator.h"
#include "log4cplus/logger.h"
#include "log4cplus/loggingmacros.h"

class Config {
public:
    Config();
    Config(string &fileName);
    ~Config();

    double Survprob();
    double Tau0();
    double Dti();
    double Deltat();

private:
    log4cplus::Logger logger;

    double deltat;
    double dti;
    double survprob;
    double tau0;

    void setDefault();
    void fromNamelist(ifstream &ifs);
    void formJson(ifstream &ifs);
};


#endif //WACOMMPLUSPLUS_CONFIG_HPP
