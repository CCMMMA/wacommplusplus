//
// Created by Raffaele Montella on 07/12/20.
//

#ifndef WACOMMPLUSPLUS_ROMSADAPTER_HPP
#define WACOMMPLUSPLUS_ROMSADAPTER_HPP

#include <string>
#include "netcdf"

// log4cplus - https://github.com/log4cplus/log4cplus
#include "log4cplus/configurator.h"
#include "log4cplus/logger.h"
#include "log4cplus/loggingmacros.h"

#include "Array.h"
#include "OceanModelAdapter.hpp"

using namespace std;
using namespace Array;
using namespace netCDF;

class Wacomm;
class ROMSAdapter: public OceanModelAdapter {
public:
    ROMSAdapter(string &fileName);

    ~ROMSAdapter();

    void process();

private:

    log4cplus::Logger logger;





    string &fileName;

    void uv2rho(Array2<double>& mask_rho, Array2<double>& mask_u, Array2<double>& mask_v,
                Array4<float>& u, Array4<float>& v);

    void wakt2rho(Array2<double>& mask_rho, Array2<double>& mask_u, Array2<double>& mask_v,
                  Array4<float>& w, Array4<float>& akt);
};


#endif //WACOMMPLUSPLUS_ROMSADAPTER_HPP
