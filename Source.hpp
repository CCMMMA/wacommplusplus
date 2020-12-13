//
// Created by Raffaele Montella on 07/12/20.
//

#ifndef WACOMMPLUSPLUS_SOURCE_HPP
#define WACOMMPLUSPLUS_SOURCE_HPP

// random https://github.com/effolkronium/random
#include "effolkronium/random.hpp"
// get base random alias which is auto seeded and has static API and internal state
using Random = effolkronium::random_static;

#include "Particles.hpp"

class Source {
public:
    Source(int id, double k, double j, double i, int start, int end, int particlesPerHour, int mode);
    ~Source();

    void emit(std::shared_ptr<Particles> particles, double tpart);
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

    double gen();
};


#endif //WACOMMPLUSPLUS_SOURCE_HPP
