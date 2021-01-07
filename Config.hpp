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

struct config_data {
    bool random{};
    double deltat{};
    double dti{};
    double survprob{};
    double tau0{};
    double sv{};
    double crid{};
};

class Config {
public:
    Config();
    explicit Config(string &fileName);
    ~Config();

    config_data *dataptr();

    string &ConfigFile();
    vector<string> &NcInputs();

    void Dry(bool value);
    bool Dry() const;



    void Random(bool value);
    bool Random() const;

    double ReductionCoefficient() const;
    double Survprob() const;
    double Tau0() const;
    double Dti() const;
    double Deltat() const;
    double SedimentationVelocity() const;

    void SaveHistory(bool value);
    bool SaveHistory() const;
    string HistoryFile() const;
    void HistoryFile(string value);

    string NcOutputRoot() const;
    void NcOutputRoot(string value);

    void SaveInput(bool value);
    bool SaveInput() const;
    string NcInputRoot() const;
    void NcInputRoot(string value);

    void UseRestart(bool value);
    bool UseRestart() const;
    string RestartFile() const;
    void RestartFile(string value);

    bool UseSources() const;
    void UseSources(bool value);
    string SourcesFile() const;
    void SourcesFile(string value);

    void StartTimeIndex(int value);
    int StartTimeIndex();

    void NumberOfInputs(int value);
    int NumberOfInputs();

    void saveAsJson(const string &fileName);
    void loadFromJson(const string &fileName);
    void loadFromNamelist(const string &fileName);


private:
    log4cplus::Logger logger;

    string configFile;

    bool dry{};

    config_data _data;

    string name;
    string institution;
    string url;

    double julianStart;
    double julianEnd;
    double julianRef;

    string oceanModel;
    string ncBasePath;
    vector<string> ncInputs;
    string ncOutputRoot;

    bool saveInput;
    string ncInputRoot;

    double timeStep;
    int nHour;
    int startTime;

    bool useRestart;
    string restartFile;
    double restartInterval;

    bool useSources;
    string sourcesFile;

    bool saveHistory;
    string historyFile;
    string initHistoryTime;
    double historyInterval;

    void setDefault();

    void namelistParseIo(ifstream &ifstream);
    void namelistParseChm(ifstream &infile);
    void namelistParseRst(ifstream &infile);
    void namelistParseHst(ifstream &infile);

};


#endif //WACOMMPLUSPLUS_CONFIG_HPP
