//
// Created by Raffaele Montella on 07/12/20.
//

#include "Sources.hpp"

Sources::Sources() {
    // Logger configuration
    logger = log4cplus::Logger::getInstance(LOG4CPLUS_TEXT("WaComM"));

    LOG4CPLUS_INFO(logger, "Empty Sources");
}

Sources::Sources(string &fileName) {
   logger = log4cplus::Logger::getInstance(LOG4CPLUS_TEXT("WaComM"));

    LOG4CPLUS_INFO(logger, "Reading sources file:" + fileName);

    std::ifstream infile(fileName);
    if(fileName.substr(fileName.find_last_of(".") + 1) == "json") {
        // The configuration is a json
    } else {
        // the configuration is a fortran style namelist
        fromNamelist(infile);
    }
}

void Sources::fromNamelist(ifstream &infile) {
    std::string line;
    while (std::getline(infile, line))
    {
        // Trim the line
        line = Utils::trim(line," \t\r");

        // Select the section parser
        if (line == "&ems") {
            // s starts with prefix &
            namelistParseEms(infile);
        }
    }
}

void Sources::namelistParseEms(ifstream &infile) {
    std::string line;

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
            if (key == "nsources") {
                int nSources = stoi(keyValues.at(1));
                for (int idx=0;idx<nSources;idx++) {
                    this->push_back(Source());
                }
            } else if (key == "id_source") {
                // The list of the input files
                vector<string> idSources;
                Utils::tokenize(keyValues.at(1),',',idSources);

                int idx=0;
                for(string idSource: idSources) {
                    idSource=Utils::trim(idSource," '");
                    this->at(idx).Id(idSource);
                    idx++;
                }
            } else if (key == "i_source") {
                // The list of the input files
                vector<string> iSources;
                Utils::tokenize(keyValues.at(1), ',', iSources);

                int idx = 0;
                for (string iSource: iSources) {
                    iSource = Utils::trim(iSource, " '");
                    this->at(idx).I(stod(iSource));
                    idx++;
                }
            } else if (key == "j_source") {
                // The list of the input files
                vector<string> jSources;
                Utils::tokenize(keyValues.at(1), ',', jSources);

                int idx = 0;
                for (string jSource: jSources) {
                    jSource = Utils::trim(jSource, " '");
                    this->at(idx).J(stod(jSource));
                    idx++;
                }
            } else if (key == "k_source") {
                // The list of the input files
                vector<string> kSources;
                Utils::tokenize(keyValues.at(1), ',', kSources);

                int idx = 0;
                for (string kSource: kSources) {
                    kSource = Utils::trim(kSource, " '");
                    this->at(idx).K(stod(kSource));
                    idx++;
                }
            } else if (key == "nPartsPerHour") {
                // The list of the input files
                vector<string> nPartsPerHours;
                Utils::tokenize(keyValues.at(1), ',', nPartsPerHours);

                int idx = 0;
                for (string nPartsPerHour: nPartsPerHours) {
                    nPartsPerHour = Utils::trim(nPartsPerHour, " '");
                    this->at(idx).ParticlesPerHour(stoi(nPartsPerHour));
                    idx++;
                }
            } else if (key == "mode") {
                // The list of the input files
                vector<string> modes;
                Utils::tokenize(keyValues.at(1), ',', modes);

                int idx = 0;
                for (string mode: modes) {
                    mode = Utils::trim(mode, " '");
                    this->at(idx).Mode(stoi(mode));
                    idx++;
                }
            } else if (key == "source_start") {
                // The list of the input files
                vector<string> sourceStarts;
                Utils::tokenize(keyValues.at(1), ',', sourceStarts);

                int idx = 0;
                for (string sourceStart: sourceStarts) {
                    sourceStart = Utils::trim(sourceStart, " '");
                    this->at(idx).Start(stoi(sourceStart));
                    idx++;
                }
            } else if (key == "source_end") {
                // The list of the input files
                vector<string> sourceEnds;
                Utils::tokenize(keyValues.at(1), ',', sourceEnds);

                int idx = 0;
                for (string sourceEnd: sourceEnds) {
                    sourceEnd = Utils::trim(sourceEnd, " '");
                    this->at(idx).End(stoi(sourceEnd));
                    idx++;
                }
            }
        }
    }
}


void Sources::save(string &fileName)
{

}

Sources::~Sources() {
}