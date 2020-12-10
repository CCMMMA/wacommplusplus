//
// Created by Raffaele Montella on 10/12/20.
//

#include "OceanModelAdapter.hpp"

OceanModelAdapter::OceanModelAdapter() {
}

std::shared_ptr<Array1<double>> OceanModelAdapter::Depth() {
    return std::shared_ptr<Array1<double>>(&depth);
}

std::shared_ptr<Array2<double>> OceanModelAdapter::Zeta() {
    return std::shared_ptr<Array2<double>>(&zeta);
}

std::shared_ptr<Array2<double>> OceanModelAdapter::Mask() {
    return std::shared_ptr<Array2<double>>(&mask);
}

std::shared_ptr<Array2<double>> OceanModelAdapter::Lon() {
    return std::shared_ptr<Array2<double>>(&lon);
}

std::shared_ptr<Array2<double>> OceanModelAdapter::Lat() {
    return std::shared_ptr<Array2<double>>(&lat);
}

std::shared_ptr<Array4<double>> OceanModelAdapter::U() {
    return std::shared_ptr<Array4<double>>(&u);
}

std::shared_ptr<Array4<double>> OceanModelAdapter::V() {
    return std::shared_ptr<Array4<double>>(&v);
}

std::shared_ptr<Array4<double>> OceanModelAdapter::W() {
    return std::shared_ptr<Array4<double>>(&w);
}

std::shared_ptr<Array4<double>> OceanModelAdapter::AKT() {
    return std::shared_ptr<Array4<double>>(&akt);
}


