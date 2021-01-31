//
// Created by Raffaele Montella on 29/01/21.
//

#ifndef WACOMMPLUSPLUS_WACOMMADAPTER_HPP
#define WACOMMPLUSPLUS_WACOMMADAPTER_HPP



#include <string>
#include "netcdf"

// log4cplus - https://github.com/log4cplus/log4cplus
#include "log4cplus/configurator.h"
#include "log4cplus/logger.h"
#include "log4cplus/loggingmacros.h"

#include "../Array.h"
#include "../OceanModelAdapter.hpp"

using namespace std;
using namespace Array;
using namespace netCDF;

class Wacomm;
class WacommAdapter: public OceanModelAdapter {
public:
    WacommAdapter(string &fileName);

    ~WacommAdapter();

    void process() override;

private:
    log4cplus::Logger logger;

    string &fileName;



};




#endif //WACOMMPLUSPLUS_WACOMMADAPTER_HPP
