//
// Created by Raffaele Montella on 08/12/20.
//

#ifndef WACOMMPLUSPLUS_CONFIG_HPP
#define WACOMMPLUSPLUS_CONFIG_HPP

#include <string>
#include <fstream>

#include "Utils.hpp"

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

    string &ConfigFile();
    vector<string> &NcInputs();

    void Dry(bool value);
    bool Dry() const;

    void UseSources(bool value);
    bool UseSources() const;

    void Random(bool value);
    bool Random() const;

    double ReductionCoefficient() const;
    double Survprob() const;
    double Tau0() const;
    double Dti() const;
    double Deltat() const;
    double SedimentationVelocity() const;

    bool UseRestart() const;
    string RestartFile() const;

    void StartTimeIndex(int value);
    int StartTimeIndex();

    void NumberOfInputs(int value);
    int NumberOfInputs();

private:
    log4cplus::Logger logger;

    string configFile;

    bool dry{};
    bool useSources{};
    bool random{};

    double deltat{};
    double dti{};
    double survprob{};
    double tau0{};
    double sv{};
    double crid{};

    vector<string> ncInputs;
    string ncOutputRoot;
    double timeStep;
    int nHour;
    int startTime;
    bool useRestart;
    string restartFile;
    double restartInterval;
    bool saveHistory;
    string historyFile;
    string initHistoryTime;
    double historyInterval;

    void setDefault();
    void fromNamelist(ifstream &ifs);
    void fromJson(ifstream &ifs);

    void namelistParseIo(ifstream &ifstream);

    void namelistParseChm(ifstream &infile);

    void namelistParseRst(ifstream &infile);

    void namelistParseHst(ifstream &infile);
};


#endif //WACOMMPLUSPLUS_CONFIG_HPP
