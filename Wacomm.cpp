//
// Created by Raffaele Montella on 12/5/20.
//

#include "Wacomm.hpp"

#include <utility>
#include "ROMS2Wacomm.hpp"

Wacomm::Wacomm(const Config &config,
               std::shared_ptr<ROMS2Wacomm> roms2Wacomm,
               Sources &sources, Particles &particles) :
               config(config),
               roms2Wacomm_(std::move(roms2Wacomm)),
               depth(roms2Wacomm_->depth),
               zeta(roms2Wacomm_->zeta),
               lon(roms2Wacomm_->lon),
               lat(roms2Wacomm_->lat),
               mask(roms2Wacomm_->mask_rho),
               u(roms2Wacomm_->ucomp),
               v(roms2Wacomm_->vcomp),
               w(roms2Wacomm_->wcomp),
               akt(roms2Wacomm_->aktcomp),
               sources(sources),
               particles(particles) {
    log4cplus::BasicConfigurator basicConfig;
    basicConfig.configure();
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
            particle.move(config,
                          ocean_time_idx,
                          depth, zeta,
                          lon, lat,
                          mask,
                          u, v, w, akt);
        }

        // Remove dead particles
        particles.erase(std::remove_if(particles.begin(), particles.end(),
                               [](Particle particle) { return !particle.isAlive(); }), particles.end());

        LOG4CPLUS_INFO(logger,"Saving restart:" << "");
        particles.save("");

        LOG4CPLUS_INFO(logger,"Evaluate concentration");
        for (Particle particle: particles) {
            if (particle.isAlive()) {
                conc(ocean_time_idx, particle.iK(), particle.iJ(), particle.iI())++;
            }
        }
    }
};
Wacomm::~Wacomm()
{
}

