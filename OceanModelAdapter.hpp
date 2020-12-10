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

using namespace Array;
using namespace std;

class OceanModelAdapter {

public:
    OceanModelAdapter();

    Array1<double> &OceanTime();
    Array1<double> &Depth();
    Array1<double> &SRho();
    Array2<double> &H();

    Array2<double> &Mask();
    Array2<double> &Lon();
    Array2<double> &Lat();
    Array2<double> &LonRad();
    Array2<double> &LatRad();

    Array3<double> &Zeta();
    Array4<double> &U();
    Array4<double> &V();
    Array4<double> &W();
    Array4<double> &AKT();

private:
    friend class Wacomm;

    Array1<double> oceanTime;
    Array1<double> depth;
    Array1<double> sRho;
    Array2<double> h;
    Array2<double> mask;
    Array2<double> lon;
    Array2<double> lat;
    Array2<double> lonRad;
    Array2<double> latRad;
    Array3<double> zeta;
    Array4<double> u, v, w;
    Array4<double> akt;
};


#endif //WACOMMPLUSPLUS_OCEANMODELADAPTER_HPP
