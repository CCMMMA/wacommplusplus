//
// Created by Raffaele Montella on 10/12/20.
//

#include "OceanModelAdapter.hpp"

OceanModelAdapter::OceanModelAdapter() {
    logger = log4cplus::Logger::getInstance(LOG4CPLUS_TEXT("WaComM"));
}

Array1<double> &OceanModelAdapter::OceanTime() { return oceanTime; }
Array1<double> &OceanModelAdapter::Depth() { return depth; }
Array1<double> &OceanModelAdapter::SRho() { return s_rho; }
Array1<double> &OceanModelAdapter::SW() { return s_w; }
Array2<double> &OceanModelAdapter::H() { return h; }
Array2<double> &OceanModelAdapter::Mask() { return mask; }
Array2<double> &OceanModelAdapter::Lon() { return lon; }
Array2<double> &OceanModelAdapter::Lat() { return lat; }
Array2<double> &OceanModelAdapter::LonRad() { return lonRad; }
Array2<double> &OceanModelAdapter::LatRad() { return latRad; }
Array3<float> &OceanModelAdapter::Zeta() { return zeta; }
Array4<float> &OceanModelAdapter::U() { return u; }
Array4<float> &OceanModelAdapter::V() { return v; }
Array4<float> &OceanModelAdapter::W() { return w; }
Array4<float> &OceanModelAdapter::AKT() { return akt; }

double OceanModelAdapter::HCorrectedByZeta(int ocean_time, int j, int i) {
    return h(j,i)+zeta(ocean_time,j,i);
}

void OceanModelAdapter::process() {
    LOG4CPLUS_ERROR(logger,"OceanModelAdapter::process must be implemented in a concrete adapter!");
    exit(-1);
}


