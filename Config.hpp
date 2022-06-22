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
    bool random;
    bool randomSources;
    double deltat;
    double dti;
    double survprob;
    double tau0;
    double sv;
    double crid;
    double sigma;
    double shoreLimit;
    int upperClosure;
    int lowerClosure;
    int horizontalClosure;
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

    void JulianStart(double value);
    double JulianStart() const;
    void JulianEnd(double value);
    double JulianEnd() const;

    int UpperClosure() const;
    int LowerClosure() const;
    int HorizontalClosure() const;

    // Random variable
    void Random(bool value);
    bool Random() const;

    // Sources random variable
    void RandomSources(bool value);
    bool RandomSources() const;

    double Sigma() const;

    double ReductionCoefficient() const;
    double Survprob() const;
    double Tau0() const;
    double Dti() const;
    double Deltat() const;
    double ShoreLimit() const;
    double SedimentationVelocity() const;

    void SaveHistory(string value);
    string SaveHistory() const;
    string HistoryRoot() const;
    void HistoryRoot(string value);

    string NcOutputRoot() const;
    void NcOutputRoot(string value);

    bool EmbeddedHistory() const;
    void EmbeddedHistroy(bool value);

    void SaveInput(bool value);
    bool SaveInput() const;
    string NcInputRoot() const;
    void NcInputRoot(string value);

    void UseRestart(bool value);
    bool UseRestart() const;
    string RestartFile() const;
    void RestartFile(string value);
    int RestartInterval() const;
    void RestartInterval(int value);

    bool MaskOutput() const;
    void MaskOutput(bool value);

    bool UseSources() const;
    void UseSources(bool value);
    string SourcesFile() const;
    void SourcesFile(string value);

    void StartTimeIndex(int value);
    int StartTimeIndex();

    void NumberOfInputs(int value);
    int NumberOfInputs();

    string OceanModel() const;
    void OceanModel(string value);

    void saveAsJson(const string &fileName);
    void loadFromJson(const string &fileName);
    void loadFromNamelist(const string &fileName);

    static const int CLOSURE_MODE_CONSTRAINT=1;
    static const int CLOSURE_MODE_KILL=2;
    static const int CLOSURE_MODE_REFLECTION=3;

private:
    log4cplus::Logger logger;

    std::map<std::string, int> dictionary;

    string configFile;

    bool dry;
    bool maskOutput;
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

    string saveHistory;
    string historyRoot;

    bool embeddedHistory;

    void setDefault();

    void namelistParseIo(ifstream &ifstream);
    void namelistParseChm(ifstream &infile);
    void namelistParseRst(ifstream &infile);
    void namelistParseHst(ifstream &infile);


};


#endif //WACOMMPLUSPLUS_CONFIG_HPP
