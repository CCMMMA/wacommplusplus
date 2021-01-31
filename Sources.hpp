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

class Sources: private vector<Source> {
public:
    Sources();
    ~Sources();

    using vector::push_back;
    using vector::operator[];
    using vector::size;
    using vector::at;
    using vector::empty;

    void loadFromNamelist(string &fileName);
    void loadFromJson(string &fileName, shared_ptr<OceanModelAdapter> oceanModelAdapter);
    void saveAsJson(string &fileName, shared_ptr<OceanModelAdapter> oceanModelAdapter);


private:
    log4cplus::Logger logger;

    void namelistParseEms(ifstream &ifstream);
};


#endif //WACOMMPLUSPLUS_SOURCES_HPP
