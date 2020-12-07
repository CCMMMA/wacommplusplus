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
        Wacomm();
        Wacomm(string fileName, Sources& sources, Particles& particles);
        ~Wacomm();

    private:
        log4cplus::Logger logger;

        Array2<double> mask_rho;

        Array4<double> ucomp, vcomp, wcomp;
        Array4<double> aktcomp;

        Sources sources;
        Particles particles;

        void load(string fileName,
                  Array2<double>& mask_u, Array2<double>& mask_v,
                  Array4<double>& u, Array4<double>& v, Array4<double>& w,
                  Array4<double>& akt);

        void uv2rho(Array2<double>& mask_u, Array2<double>& mask_v,
                    Array4<double>& u, Array4<double>& v);

        void wakt2rho(Array2<double>& mask_u, Array2<double>& mask_v,
                      Array4<double>& w, Array4<double>& akt);
};


#endif //WACOMMPLUSPLUS_WACOMM_HPP
