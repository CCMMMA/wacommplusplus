//
// Created by Raffaele Montella on 07/12/20.
//

#ifndef WACOMMPLUSPLUS_ROMS2WACOMM_HPP
#define WACOMMPLUSPLUS_ROMS2WACOMM_HPP

#include <string>
#include "netcdf"

// log4cplus - https://github.com/log4cplus/log4cplus
#include "log4cplus/configurator.h"
#include "log4cplus/logger.h"
#include "log4cplus/loggingmacros.h"

#include "Array.h"

using namespace std;
using namespace Array;
using namespace netCDF;

class Wacomm;
class ROMS2Wacomm {
public:
    ROMS2Wacomm(string &fileName);

    ~ROMS2Wacomm();

    void process();

private:
    friend class Wacomm;
    log4cplus::Logger logger;

    Array1<double> depth;
    Array2<double> zeta;
    Array2<double> mask_rho;
    Array2<double> lon;
    Array2<double> lat;

    Array4<double> ucomp, vcomp, wcomp;
    Array4<double> aktcomp;

    string &fileName;

    void uv2rho(Array2<double>& mask_u, Array2<double>& mask_v,
                Array4<double>& u, Array4<double>& v);

    void wakt2rho(Array2<double>& mask_u, Array2<double>& mask_v,
                  Array4<double>& w, Array4<double>& akt);
};


#endif //WACOMMPLUSPLUS_ROMS2WACOMM_HPP
