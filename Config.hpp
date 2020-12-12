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
    explicit Config(string &fileName);
    ~Config();

    void setDry(bool value);
    [[nodiscard]] bool isDry() const;

    [[nodiscard]] double ReductionCoefficient() const;
    [[nodiscard]] double Survprob() const;
    [[nodiscard]] double Tau0() const;
    [[nodiscard]] double Dti() const;
    [[nodiscard]] double Deltat() const;
    [[nodiscard]]double SedimentationVelocity() const;

private:
    log4cplus::Logger logger;

    bool dry{};

    double deltat{};
    double dti{};
    double survprob{};
    double tau0{};
    double sv{};
    double crid{};

    void setDefault();
    void fromNamelist(ifstream &ifs);
    void formJson(ifstream &ifs);
};


#endif //WACOMMPLUSPLUS_CONFIG_HPP
