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

    std::shared_ptr<Array1<double>> Depth();
    std::shared_ptr<Array2<double>> Zeta();
    std::shared_ptr<Array2<double>> Mask();
    std::shared_ptr<Array2<double>> Lon();
    std::shared_ptr<Array2<double>> Lat();

    std::shared_ptr<Array4<double>> U();
    std::shared_ptr<Array4<double>> V();
    std::shared_ptr<Array4<double>> W();
    std::shared_ptr<Array4<double>> AKT();
private:
    friend class Wacomm;

    Array1<double> depth;
    Array2<double> zeta;
    Array2<double> mask;
    Array2<double> lon;
    Array2<double> lat;

    Array4<double> u, v, w;
    Array4<double> akt;
};


#endif //WACOMMPLUSPLUS_OCEANMODELADAPTER_HPP
