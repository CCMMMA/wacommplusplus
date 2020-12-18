//
// Created by Raffaele Montella on 07/12/20.
//

#include "Source.hpp"

Source::Source(string id, double k, double j, double i, int start, int end, int particlesPerHour, int mode):
    id(id),k(k),j(j),i(i),start(start),end(end),particlesPerHour(particlesPerHour),mode(mode) {

    logger = log4cplus::Logger::getInstance(LOG4CPLUS_TEXT("WaComM"));
}

Source::~Source() = default;

void Source::emit(const std::shared_ptr<Config>& config, std::shared_ptr<Particles> particles, double tpart) {
    // Check if the source is active
    if (mode>0) {
        // Check if the particle have to be released
        if (mode == 1) {
            for (int i = 0; i < particlesPerHour; i++) {
                double kk = k;
                double jj = j;
                double ii = i;
                if (config->Random()) {
                    kk = k + gen() * 0.5 - 0.25;
                    jj = j + gen() * 0.5 - 0.25;
                    ii = i + gen() * 0.5 - 0.25;
                }
                particles->push_back(Particle(kk, jj, ii, tpart));
            }
        }
    }
}

double Source::gen()
{
    return Random::get<double>(0.0, 1.0);
}

Source::Source() {
    logger = log4cplus::Logger::getInstance(LOG4CPLUS_TEXT("WaComM"));
}

void Source::Id(string value) { id = value; }
void Source::K(double value) { k = value; }
void Source::J(double value) { j = value; }
void Source::I(double value) { i = value; }
void Source::Start(int value) { start = value; }
void Source::End(int value) { end = value; }
void Source::ParticlesPerHour(int value) { particlesPerHour = value; }
void Source::Mode(int value) { mode = value; }
