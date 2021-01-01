//
// Created by Raffaele Montella on 12/5/20.
//

#include "Particles.hpp"

Particles::Particles() {

    // Logger configuration
    log4cplus::BasicConfigurator basicConfig;
    basicConfig.configure();
    logger = log4cplus::Logger::getInstance(LOG4CPLUS_TEXT("WaComM"));

    LOG4CPLUS_DEBUG(logger, "Empty Particles");
}

Particles::Particles(string fileName) {

    logger = log4cplus::Logger::getInstance(LOG4CPLUS_TEXT("WaComM"));

    LOG4CPLUS_DEBUG(logger, "Reading restart file:"+fileName);

    std::ifstream infile(fileName);

    int count=0;
    std::string line;
    while (std::getline(infile, line))
    {
        std::istringstream iss(line);
        if (count > 0)
        {
            double k, j, i, health, tpart;
            if (!(iss >> i >> j >> k >> health >> tpart )) { break; } // error
            Particle particle(k, j, i, health, tpart, -1);
            push_back(particle);
        }
        count++;
    }

}

void Particles::save(const string& fileName)
{
    LOG4CPLUS_DEBUG(logger, "Saving restart file:"+fileName);
    std::ofstream outfile(fileName);
    int count = size();
    outfile << "\t" << count << endl;
    for(int idx=0; idx<count; idx++) {
        Particle particle = at(idx);
        outfile << particle.to_string() << endl;
    }
}

Particles::~Particles() = default;
