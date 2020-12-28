//
// Created by Raffaele Montella on 10/12/20.
//

#ifndef WACOMMPLUSPLUS_OCEANMODELADAPTER_HPP
#define WACOMMPLUSPLUS_OCEANMODELADAPTER_HPP

// log4cplus - https://github.com/log4cplus/log4cplus
#include "log4cplus/configurator.h"
#include "log4cplus/logger.h"
#include "log4cplus/loggingmacros.h"

#include "Array.h"
#include "netcdf"

using namespace Array;
using namespace std;
using namespace netCDF;

struct oceanmodel_data {
    Array1<double> oceanTime;
    Array1<double> sRho;
    Array1<double> sW;
    Array1<double> depth;
    Array2<double> h;
    Array2<double> mask;
    Array2<double> lon;
    Array2<double> lat;
    Array2<double> lonRad;
    Array2<double> latRad;
    Array3<float> zeta;
    Array4<float> u, v, w;
    Array4<float> akt;
};

class OceanModelAdapter {

public:
    OceanModelAdapter();

    void saveAsNetCDF(std::string &fileName);

    virtual void process();

    oceanmodel_data *dataptr();

    Array1<double> &OceanTime();
    Array1<double> &SRho();
    Array1<double> &SW();


    Array2<double> &Mask();
    Array2<double> &Lon();
    Array2<double> &Lat();

    Array2<double> &LonRad();
    Array2<double> &LatRad();
    Array1<double> &Depth();
    Array2<double> &H();
    Array3<float> &Zeta();
    Array4<float> &U();
    Array4<float> &V();
    Array4<float> &W();
    Array4<float> &AKT();

private:
    log4cplus::Logger logger;

    friend class WacommPlusPlus;
    oceanmodel_data _data{};

};


#endif //WACOMMPLUSPLUS_OCEANMODELADAPTER_HPP
