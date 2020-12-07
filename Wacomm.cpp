//
// Created by Raffaele Montella on 12/5/20.
//

#include "Wacomm.hpp"



Wacomm::Wacomm() {
    log4cplus::BasicConfigurator config;
    config.configure();
    logger = log4cplus::Logger::getInstance(LOG4CPLUS_TEXT("WaComM"));
}

Wacomm::Wacomm(double dti, Array4<double>& u, Array4<double>& v, Array4<double>& w, Array4<double> &akt, Sources& sources, Particles& particles) {
    log4cplus::BasicConfigurator config;
    config.configure();
    logger = log4cplus::Logger::getInstance(LOG4CPLUS_TEXT("WaComM"));

    this->dti = dti;

    this->sources=sources;
    this->particles=particles;

    this->u = u;
    this->v = v;
    this->w = w;
    this->akt = akt;

    size_t ocean_time=u.Nx();
    size_t s_rho=u.Ny();
    size_t s_w=w.Ny();
    size_t eta_rho=mask.Nx();
    size_t xi_rho=mask.Ny();


    LOG4CPLUS_INFO(logger,"ocean_time:" + std::to_string(ocean_time));
    LOG4CPLUS_INFO(logger,"s_rho:" + std::to_string(s_rho));
    LOG4CPLUS_INFO(logger,"s_w:" + std::to_string(s_w));
    LOG4CPLUS_INFO(logger,"eta_rho:" + std::to_string(eta_rho) + " xi_rho:" + std::to_string(xi_rho));


    conc.Allocate(ocean_time,s_rho,eta_rho,xi_rho);
}

void Wacomm::run()
{
    for(Source source: sources) {
        source.emit(particles);
    }

    size_t ocean_time=u.Nx();
    for (int ocean_time_idx=0; ocean_time_idx<ocean_time; ocean_time++) {

        for (Particle particle: particles) {
            particle.move(ocean_time_idx, mask, u[ocean_time_idx], v[ocean_time_idx], w[ocean_time_idx], akt[ocean_time_idx]);
        }

        particles.save("");

        for (Particle particle: particles) {
            if (particle.isAlive()) {
                conc(ocean_time_idx, particle.k(), particle.j(), particle.i())++;
            }
        }
    }
};
Wacomm::~Wacomm()
{
}

