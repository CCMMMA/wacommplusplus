//
// Created by Raffaele Montella on 07/12/20.
//

#include "Source.hpp"

Source::Source(int id, int j, int i, int start, int end, int particlesPerHour, int mode) {
    _id=id;
    _j=0;
    _i=i;
    _start=start;
    _end=end;
    _particlesPerHour=particlesPerHour;
    _mode=mode;
}

Source::~Source()
{

}

void Source::emit(Particles &particles) {

}
