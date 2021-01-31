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
    Source();
    Source(string id, double k, double j, double i, double start, double end, int particlesPerHour, int mode);
    ~Source();

    void emit(const std::shared_ptr<Config>& config, std::shared_ptr<Particles> particles, double currentOceanTime);

    string Id();
    double K();
    double J();
    double I();
    double Start();
    double End();
    int ParticlesPerHour();
    int Mode();

    void Id(string value);
    void K(double value);
    void J(double value);
    void I(double value);
    void Start(double value);
    void End(double value);
    void ParticlesPerHour(int value);
    void Mode(int value);



private:
    log4cplus::Logger logger;

    string id;
    double k;
    double j;
    double i;
    double start;
    double end;
    int particlesPerHour;
    int mode;

    static double gen();
};


#endif //WACOMMPLUSPLUS_SOURCE_HPP
