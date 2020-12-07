//
// Created by Raffaele Montella on 07/12/20.
//

#ifndef WACOMMPLUSPLUS_SOURCE_HPP
#define WACOMMPLUSPLUS_SOURCE_HPP


#include "Particles.hpp"

class Source {
private:
    int _id;
    int _j;
    int _i;
    int _start;
    int _end;
    int _particlesPerHour;
    int _mode;

public:
    Source(int id, int j, int i, int start, int end, int particlesPerHour, int mode);
    ~Source();

    void emit(Particles &particles);
};


#endif //WACOMMPLUSPLUS_SOURCE_HPP
