//
// Created by Raffaele Montella on 07/12/20.
//

#ifndef WACOMMPLUSPLUS_SOURCE_HPP
#define WACOMMPLUSPLUS_SOURCE_HPP


#include "Particles.hpp"

class Source {
public:
    Source(int id, double k, double j, double i, int start, int end, int particlesPerHour, int mode);
    ~Source();

    void emit(std::shared_ptr<Particles> particles);
private:
    log4cplus::Logger logger;

    int id;
    double k;
    double j;
    double i;
    int start;
    int end;
    int particlesPerHour;
    int mode;


};


#endif //WACOMMPLUSPLUS_SOURCE_HPP
