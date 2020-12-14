//
// Created by Raffaele Montella on 08/12/20.
//

#include "Config.hpp"

Config::Config() = default;

Config::Config(string &fileName): configFile(fileName) {
    setDefault();

    log4cplus::BasicConfigurator basicConfig;
    basicConfig.configure();
    logger = log4cplus::Logger::getInstance(LOG4CPLUS_TEXT("WaComM"));

    LOG4CPLUS_INFO(logger, "Reading config file:" + fileName);
    std::ifstream infile(fileName);
    if(fileName.substr(fileName.find_last_of(".") + 1) == "json") {
        // The configuration is a json
        fromJson(infile);
    } else {
        // the configuration is a fortran style namelist
        fromNamelist(infile);
    }
    LOG4CPLUS_INFO(logger, "...done!");

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

    // File system root where the outputs will be saved ( default current directory )
    ncOutputRoot = "";

    timeStep = 3600;

    // Number of hours to calculate ( default 1)
    nHour = 1;

    // Starting time (default 0)
    startTime = 0;

    // Use restart (default false)
    useRestart = false;

    // Name of the file used for restarts
    restartFile = "WACOMM_rst_.txt";

    // Restart interval
    restartInterval = 3600;

    // Save history (to be used as restart) default true
    saveHistory = true;

    // Default file name for the history
    historyFile = "WACOMM_rst_";

    initHistoryTime = "";

    // Save history each seconds (default 3600, 1h)
    historyInterval = 3600;

    // Use sources (defult true)
    useSources = true;
}

void Config::fromNamelist(ifstream &infile) {
    std::string line;
    while (std::getline(infile, line))
    {
        // Trim the line
        line = Utils::trim(line," \t\r");

        // Select the section parser
        if (line == "&io") {
            // s starts with prefix &
            namelistParseIo(infile);
        } else if (line == "&chm") {
            // s starts with prefix &
            namelistParseChm(infile);
        } if (line == "&rst") {
            // s starts with prefix &
            namelistParseRst(infile);
        } if (line == "&hst") {
            // s starts with prefix &
            namelistParseHst(infile);
        }
    }
}

void Config::namelistParseIo(ifstream &infile) {
    std::string line;

    LOG4CPLUS_INFO(logger, "Parsing io section");

    // Start the input/output section
    while (std::getline(infile, line))
    {
        // Trim the line
        line = Utils::trim(line, " \t\r");

        // Check for the end of the section
        if (line=="/") { break; }

        vector<string> keyValues;
        Utils::tokenize(line,'=',keyValues);

        if (keyValues.size()==2) {
            string key=Utils::trim(keyValues.at(0)," \t\r");
            if (key == "nc_inputs") {
                // The list of the input files
                vector<string> ncInputs;
                Utils::tokenize(keyValues.at(1),',',ncInputs);

                for(string ncInput: ncInputs) {
                    ncInput=Utils::trim(ncInput," '");
                    this->ncInputs.push_back(ncInput);
                }
            } else if (key == "nc_output_root") {
                this->ncOutputRoot = keyValues.at(1);
            } else if (key == "starttime") {
                this->startTime = stoi(keyValues.at(1));
            } else if (key == "nHour") {
                this->nHour = stoi(keyValues.at(1));
            } else if (key == "timestep") {
                this->timeStep = stod(keyValues.at(1));
            }
        }
    }
}

void Config::namelistParseChm(ifstream &infile) {
    std::string line;

    LOG4CPLUS_INFO(logger, "Parsing chm section");

    // Start the Chem section
    while (std::getline(infile, line))
    {
        // Trim the line
        line = Utils::trim(line, " \t\r");

        // Check for the end of the section
        if (line=="/") { break; }

        vector<string> keyValues;
        Utils::tokenize(line,'=',keyValues);

        if (keyValues.size()==2) {
            string key=Utils::trim(keyValues.at(0)," \t\r");
            if (key == "tau0") {
                this->tau0 = stod(keyValues.at(1));
            } else if (key == "survprob") {
                this->survprob = stod(keyValues.at(1));
            }
        }
    }
}

void Config::namelistParseRst(ifstream &infile) {
    std::string line;

    LOG4CPLUS_INFO(logger, "Parsing rst section");

    // Start the Restart section
    while (std::getline(infile, line))
    {
        // Trim the line
        line = Utils::trim(line, " \t\r");

        // Check for the end of the section
        if (line=="/") { break; }

        vector<string> keyValues;
        Utils::tokenize(line,'=',keyValues);

        if (keyValues.size()==2) {
            string key=Utils::trim(keyValues.at(0)," \t\r");
            if (key == "restart") {
                keyValues.at(1) = Utils::trim(keyValues.at(1), " \t");
                if (keyValues.at(1) == ".true.") {
                    this->useRestart=true;
                } else {
                    this->useRestart=false;
                }
            } else if (key == "restartfile") {
                this->restartFile = Utils::trim(keyValues.at(1)," '\t\r");;
            } else if (key == "interval") {
                this->restartInterval = stod(keyValues.at(1));
            }
        }
    }
}

void Config::namelistParseHst(ifstream &infile) {
    std::string line;

    LOG4CPLUS_INFO(logger, "Parsing hst section");

    // Start the History section
    while (std::getline(infile, line))
    {
        // Trim the line
        line = Utils::trim(line, " \t\r");

        // Check for the end of the section
        if (line=="/") { break; }

        vector<string> keyValues;
        Utils::tokenize(line,'=',keyValues);

        if (keyValues.size()==2) {
            string key=Utils::trim(keyValues.at(0)," \t\r");
            if (key == "history") {
                keyValues.at(1) = Utils::trim(keyValues.at(1), " \t");
                if (keyValues.at(1) == ".true.") {
                    this->saveHistory=true;
                } else {
                    this->saveHistory=false;
                }
            } else if (key == "historyfile") {
                this->historyFile = keyValues.at(1);
            } else if (key == "initthsttime") {
                this->initHistoryTime = keyValues.at(1);
            } else if (key == "outfreq") {
                this->historyInterval = stod(keyValues.at(1));
            }
        }
    }
}

void Config::fromJson(ifstream &infile) {

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

bool Config::Dry() const {
    return dry;
}

void Config::Dry(bool value) {
    dry=value;
}

double Config::ReductionCoefficient() const {
    return crid;
}

string Config::RestartFile() const {
    return restartFile;
}

bool Config::UseRestart() const {
    return useRestart;
}

string &Config::ConfigFile() {
    return configFile;
}

vector<string> &Config::NcInputs() {
    return ncInputs;
}

void Config::StartTimeIndex(int value) {
    startTime=value;
}

int Config::StartTimeIndex() {
    return startTime;
}

void Config::NumberOfInputs(int value) {
    nHour=value;
}

int Config::NumberOfInputs() {
    return nHour;
}

bool Config::UseSources() const {
    return useSources;
}

void Config::UseSources(bool value) {
    useSources=value;
}


