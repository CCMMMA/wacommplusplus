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

    LOG4CPLUS_INFO(logger, "Reading restart file:" + fileName);

    std::ifstream infile(fileName);
}

void Sources::save(string &fileName)
{

}

Sources::~Sources() {
}