//
// Created by Raffaele Montella on 07/12/20.
//

#include "Source.hpp"
#include <random>

Source::Source(string id, double k, double j, double i, double start, double end, int particlesPerHour, int mode):
    id(id),k(k),j(j),i(i),start(start),end(end),particlesPerHour(particlesPerHour),mode(mode) {

    logger = log4cplus::Logger::getInstance(LOG4CPLUS_TEXT("WaComM"));
}

Source::~Source() = default;

void Source::emit(const std::shared_ptr<Config>& config, std::shared_ptr<Particles> particles, double currentOceanTime) {

    // Create a random number generator
    std::default_random_engine generator;

    // Create a distribution probability with mean=0 and stddev=0.25
    std::normal_distribution<double> distribution(0.0,0.25);

    // Check if the source is active
    if (mode>0) {
        unsigned long id=0;
        if (particles->size()>0) {
            id = particles->at(particles->size() - 1).Id() + 1;
        }
        // Check if the particle have to be released
        if ((mode == 1) && (start<0 || start <= currentOceanTime) && (end<0 || end >= currentOceanTime) ) {

            // Release the particles
            for (int idx = 0; idx < particlesPerHour; idx++) {
                double kk = k;
                double jj = j;
                double ii = i;

                if (config->RandomSources()) {

                    kk = k + distribution(generator);
                    jj = j + distribution(generator);
                    ii = i + distribution(generator);
                }

                particles->push_back(Particle(id, kk, jj, ii, currentOceanTime));
                id++;
            }
        }
    }
}

Source::Source() {
    logger = log4cplus::Logger::getInstance(LOG4CPLUS_TEXT("WaComM"));
}

void Source::Id(string value) { id = value; }
void Source::K(double value) { k = value; }
void Source::J(double value) { j = value; }
void Source::I(double value) { i = value; }
void Source::Start(double value) { start = value; }
void Source::End(double value) { end = value; }
void Source::ParticlesPerHour(int value) { particlesPerHour = value; }
void Source::Mode(int value) { mode = value; }

string Source::Id() { return id; }

double Source::K() { return k; }
double Source::J() { return j; }
double Source::I() { return i; }

double Source::Start() { return start; }
double Source::End() { return end; }

int Source::ParticlesPerHour() { return particlesPerHour; }

int Source::Mode() { return mode; }