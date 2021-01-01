//
// Created by Raffaele Montella on 07/12/20.
//

#ifndef WACOMMPLUSPLUS_SOURCES_HPP
#define WACOMMPLUSPLUS_SOURCES_HPP

#include <fstream>

// log4cplus - https://github.com/log4cplus/log4cplus
#include "log4cplus/configurator.h"
#include "log4cplus/logger.h"
#include "log4cplus/loggingmacros.h"

#include "Source.hpp"

using namespace std;

class Sources: public vector<Source> {
public:
    Sources();
    ~Sources();

    void loadFromNamelist(string &fileName);
    void loadFromJson(string &fileName, std::shared_ptr<OceanModelAdapter> oceanModelAdapter);
    void saveAsJson(string &fileName, std::shared_ptr<OceanModelAdapter> oceanModelAdapter);

private:
    log4cplus::Logger logger;

    void namelistParseEms(ifstream &ifstream);
};


#endif //WACOMMPLUSPLUS_SOURCES_HPP
