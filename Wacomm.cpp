//
// Created by Raffaele Montella on 12/5/20.
//

#include "Wacomm.hpp"





Wacomm::Wacomm(double dti,
               Array2<double> &mask,
               Array4<double> &u, Array4<double> &v, Array4<double> &w, Array4<double> &akt,
               Sources &sources, Particles &particles) :
               dti(dti),
               mask(mask),
               u(u),v(v),w(w),akt(akt),
               sources(sources), particles(particles) {
    log4cplus::BasicConfigurator config;
    config.configure();
    logger = log4cplus::Logger::getInstance(LOG4CPLUS_TEXT("WaComM"));

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
    size_t ocean_time=u.Nx();

    for (int ocean_time_idx=0; ocean_time_idx<ocean_time; ocean_time_idx++) {
        LOG4CPLUS_INFO(logger,"Running on:" << ocean_time_idx);
        LOG4CPLUS_INFO(logger,"Sources:" << sources.size());
        for(Source source: sources) {
            source.emit(particles);
        }

        LOG4CPLUS_INFO(logger,"Particles:" << particles.size());
        for (Particle particle: particles) {
            particle.move(ocean_time_idx, mask, u[ocean_time_idx], v[ocean_time_idx], w[ocean_time_idx], akt[ocean_time_idx]);
        }

        LOG4CPLUS_INFO(logger,"Saving restart:" << "");
        particles.save("");

        LOG4CPLUS_INFO(logger,"Evaluate concentration");
        for (Particle particle: particles) {
            if (particle.isAlive()) {
                conc(ocean_time_idx, particle.K(), particle.J(), particle.I())++;
            }
        }
    }
};
Wacomm::~Wacomm()
{
}

