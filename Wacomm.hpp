//
// Created by Raffaele Montella on 12/5/20.
//

#ifndef WACOMMPLUSPLUS_WACOMM_HPP
#define WACOMMPLUSPLUS_WACOMM_HPP

// log4cplus - https://github.com/log4cplus/log4cplus
#include "log4cplus/configurator.h"
#include "log4cplus/logger.h"
#include "log4cplus/loggingmacros.h"

#include "Sources.hpp"
#include "Particles.hpp"
#include "Array.h"

using namespace std;
using namespace Array;

class Wacomm {
    public:

        Wacomm(double dti, Array2<double> &mask, Array4<double> &u, Array4<double> &v, Array4<double> &w, Array4<double> &akt, Sources &sources, Particles &particles);
        ~Wacomm();

        void run();

    private:
        log4cplus::Logger logger;

        double dti;

        Array2<double> &mask;

        Array4<double> &u, &v, &w;
        Array4<double> &akt;
        Array4<double> conc;

        Sources &sources;
        Particles &particles;


};


#endif //WACOMMPLUSPLUS_WACOMM_HPP
