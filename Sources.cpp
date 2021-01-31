//
// Created by Raffaele Montella on 07/12/20.
//

#include "Sources.hpp"
#include <iostream>
#include <fstream>
#include <nlohmann/json.hpp>
#include <iomanip>

// for convenience
using json = nlohmann::json;

Sources::Sources() {
    // Logger configuration
    logger = log4cplus::Logger::getInstance(LOG4CPLUS_TEXT("WaComM"));

    LOG4CPLUS_DEBUG(logger, "Empty Sources");
}

void Sources::loadFromNamelist(string &fileName) {
    LOG4CPLUS_DEBUG(logger, "Reading from namelist:" + fileName);

    std::ifstream infile(fileName);

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
                    this->at(idx).Start(stof(sourceStart));
                    idx++;
                }
            } else if (key == "source_end") {
                // The list of the input files
                vector<string> sourceEnds;
                Utils::tokenize(keyValues.at(1), ',', sourceEnds);

                int idx = 0;
                for (string sourceEnd: sourceEnds) {
                    sourceEnd = Utils::trim(sourceEnd, " '");
                    this->at(idx).End(stof(sourceEnd));
                    idx++;
                }
            }
        }
    }
}


void Sources::saveAsJson(string &fileName, shared_ptr<OceanModelAdapter> oceanModelAdapter)
{
    json features = json::array();

    for (auto source:*this) {

        if (source.J() < 0 || source.I() < 0 || source.K() > 0) {
            continue;
        }

        json coordinates = json::array();
        double dep, lat, lon;
        oceanModelAdapter->kji2deplatlon(source.K(), source.J(), source.I(), dep, lat, lon);

        if (dep == 1e37 || lat == 1e37 || lon == 1e37) {
            continue;
        }
        coordinates.push_back(lon);
        coordinates.push_back(lat);

        json properties = {
                { "id", source.Id()},
                { "k", source.K()},
                { "j", source.J()},
                { "i", source.I()},
                { "start", source.Start()},
                { "end", source.Start()},
                { "particlesPerHour", source.ParticlesPerHour()},
                { "mode", source.Mode()},
                { "depth", dep}
        };

        json geometry = {
                { "type", "Point"},
                { "coordinates", coordinates }
        };

        json feature = {
                {"type", "Feature"},
                {"properties", properties},
                {"geometry", geometry}
        };

        features.push_back(feature);
    }

    json featureCollection = {
            {"type", "FeatureCollection"},
            {"features", features}
    };

    // write prettified JSON to another file
    std::ofstream o(fileName);
    o << std::setw(4) << featureCollection << std::endl;
}

Sources::~Sources() {
}

void Sources::loadFromJson(string &fileName, shared_ptr<OceanModelAdapter> oceanModelAdapter) {

    LOG4CPLUS_INFO(logger, "Reading from json:" << fileName);

    std::ifstream infile(fileName);

    try {
        json featureCollection;
        infile >> featureCollection;

        if (featureCollection.contains("features") && featureCollection["features"].is_array()) {

            int count = 1;
            for (auto feature:featureCollection["features"]) {
                double k = 0, j = 1e37, i = 1e37;
                double startOceanTime = -1, endOceanTime = -1;
                double dep = 0;
                int mode = 1, particlesPerHour = 100;
                string id = "source_" + to_string(count);
                if (feature.contains("properties")) {
                    auto properties = feature["properties"];
                    if (properties.contains("id")) { id = properties["id"]; }
                    if (properties.contains("start")) { startOceanTime = properties["start"]; }
                    if (properties.contains("end")) { endOceanTime = properties["end"]; }
                    if (properties.contains("k")) { k = atof(to_string(properties["k"]).c_str()); }
                    if (properties.contains("j")) { j = atof(to_string(properties["j"]).c_str()); }
                    if (properties.contains("i")) { i = atof(to_string(properties["i"]).c_str()); }
                    if (properties.contains("dep")) { dep = atof(to_string(properties["dep"]).c_str()); }
                    if (properties.contains("particlesPerHour")) {
                        particlesPerHour = stoi(to_string(properties["particlesPerHour"]));
                    }
                }

                // convert lat/lon and depth in k,j,i
                if (j == 1e37 || i == 1e37) {
                    double lat = 1e37, lon = 1e37;
                    if (feature.contains("geometry")) {
                        auto geometry = feature["geometry"];
                        if (geometry.contains("type") && geometry["type"] == "Point") {
                            if (geometry.contains("coordinates") && geometry["coordinates"].is_array()) {
                                auto coordinates = geometry["coordinates"];
                                if (coordinates.size() >= 2) {
                                    lon = coordinates[0];
                                    lat = coordinates[1];
                                    if (coordinates.size() == 3 && dep == 0) {
                                        dep = coordinates[2];
                                    }
                                    //oceanModelAdapter->deplatlon2kji(dep, lat, lon, k, j, i);
                                }
                            }
                        }
                    }
                }

                this->push_back(Source(id, k, j, i, startOceanTime, endOceanTime, particlesPerHour, mode));
                count++;
            }
        }
    } catch (const nlohmann::json::parse_error& e) {
        LOG4CPLUS_ERROR(logger,e.what());
    }

}
