//
// Created by Raffaele Montella on 10/12/20.
//

#include "OceanModelAdapter.hpp"

OceanModelAdapter::OceanModelAdapter() {
}

Array1<double> &OceanModelAdapter::OceanTime() { return oceanTime; }
Array1<double> &OceanModelAdapter::Depth() { return depth; }
Array1<double> &OceanModelAdapter::SRho() { return sRho; }
Array2<double> &OceanModelAdapter::H() { return h; }
Array2<double> &OceanModelAdapter::Mask() { return mask; }
Array2<double> &OceanModelAdapter::Lon() { return lon; }
Array2<double> &OceanModelAdapter::Lat() { return lat; }
Array2<double> &OceanModelAdapter::LonRad() { return lonRad; }
Array2<double> &OceanModelAdapter::LatRad() { return latRad; }
Array3<double> &OceanModelAdapter::Zeta() { return zeta; }
Array4<double> &OceanModelAdapter::U() { return u; }
Array4<double> &OceanModelAdapter::V() { return v; }
Array4<double> &OceanModelAdapter::W() { return w; }
Array4<double> &OceanModelAdapter::AKT() { return akt; }
