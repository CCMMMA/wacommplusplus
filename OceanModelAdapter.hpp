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
    Array1<double> depthIntervals;
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

    virtual void process()=0;

    void kji2deplatlon(double k, double j, double i, double &dep, double &lat, double &lon);
    void deplatlon2kji(double dep, double lat, double lon, double &k, double &j, double &i);
    void allocateMemory(size_t ocean_time, size_t s_rho, size_t s_w, size_t eta_rho, size_t xi_rho);

    oceanmodel_data *dataptr();

    Array1<double> &OceanTime();
    Array1<double> &SRho();
    Array1<double> &SW();

    Array2<double> &Mask();
    Array2<double> &Lon();
    Array2<double> &Lat();

    Array2<double> &LonRad();
    Array2<double> &LatRad();
    Array1<double> &DepthIntervals();
    Array2<double> &H();
    Array3<float> &Zeta();
    Array4<float> &U();
    Array4<float> &V();
    Array4<float> &W();
    Array4<float> &AKT();

    Array1<double> &Latitude();
    Array1<double> &Longitude();
    Array1<double> &Depth();

private:
    log4cplus::Logger logger;
    oceanmodel_data _data;

    Array1<double> latitude;
    Array1<double> longitude;
    Array1<double> depth;

    double sgn(double a);
};


#endif //WACOMMPLUSPLUS_OCEANMODELADAPTER_HPP
