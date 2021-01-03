//
// Created by Raffaele Montella on 12/5/20.
//

#include "Particles.hpp"
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <nlohmann/json.hpp>

// for convenience
using json = nlohmann::json;

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

    loadFromTxt(fileName);

}

void Particles::saveAsTxt(const string& fileName)
{
    LOG4CPLUS_DEBUG(logger, "Saving restart file:"+fileName);
    std::ofstream outfile(fileName);
    int count = size();
    outfile << "\t" << count << endl;
    for(int idx=0; idx<count; idx++) {
        Particle particle = at(idx);
        outfile << particle.Id() << " " << particle.I() << " "  << particle.J() << " " <<  particle.K() << " " << particle.Health() << " " <<  particle.Age() << " " << particle.Time() << endl;
    }
}

void Particles::saveAsJson(const string &fileName, std::shared_ptr<OceanModelAdapter> oceanModelAdapter)
{
    json features = json::array();

    for (auto particle:*this) {

        if (particle.J() < 0 || particle.I() < 0 || particle.K() > 0) {
            continue;
        }

        json coordinates = json::array();
        double dep, lat, lon;
        oceanModelAdapter->kji2deplatlon(particle.K(), particle.J(), particle.I(), dep, lat, lon);

        if (dep == 1e37 || lat == 1e37 || lon == 1e37) {
            continue;
        }
        coordinates.push_back(lon);
        coordinates.push_back(lat);

        json properties = {
                { "id", particle.Id()},
                { "k", particle.K()},
                { "j", particle.J()},
                { "i", particle.I()},
                { "health", particle.Health()},
                { "age", particle.Age()},
                { "time", particle.Time()},
                { "depth", dep}
        };

        json geometry = {
                { "type", "Point"},
                { "coordinates", coordinates }
        };

        json feature = {
                {"type", "Feature"},
                {"properties", properties},
                {"geometry", geometry}
        };

        features.push_back(feature);
    }

    json featureCollection = {
            {"type", "FeatureCollection"},
            {"features", features}
    };

    // write prettified JSON to another file
    std::ofstream o(fileName);
    o << std::setw(4) << featureCollection << std::endl;
}

void Particles::saveAsNetCDF(const string &fileName, std::shared_ptr<OceanModelAdapter> oceanModelAdapter)
{
    // https://github.com/NOAA-ORR-ERD/nc_particles

    size_t ocean_time=oceanModelAdapter->OceanTime().Nx();

    Array1<unsigned long> id(this->size());
    Array2<double> lat(this->size(),ocean_time);
    Array2<double> lon(this->size(),ocean_time);
    Array2<double> dep(this->size(),ocean_time);
    Array2<double> j(this->size(),ocean_time);
    Array2<double> i(this->size(),ocean_time);
    Array2<double> k(this->size(),ocean_time);
    Array2<double> health(this->size(),ocean_time);
    Array2<double> age(this->size(),ocean_time);
    Array2<double> time(this->size(),ocean_time);

    LOG4CPLUS_INFO(logger,"Preparing NetCDF...");

    int count=0;
    for (const auto& particle:*this) {

        if (particle.J() < 0 || particle.I() < 0 || particle.K() > 0) {
            continue;
        }

        json coordinates = json::array();
        double depth, latitude, longitude;
        oceanModelAdapter->kji2deplatlon(particle.K(), particle.J(), particle.I(), depth, latitude, longitude);

        if (depth == 1e37 || latitude == 1e37 || longitude == 1e37) {
            continue;
        }

        id(count)=particle.Id();
        lat(count,0)=latitude;
        lon(count,0)=longitude;
        dep(count,0)=depth;
        j(count,0)=particle.J();
        i(count,0)=particle.I();
        k(count,0)=particle.K();
        health(count,0)=particle.Health();
        age(count,0)=particle.Age();
        time(count,0)=particle.Time();

        count++;
    }


    LOG4CPLUS_INFO(logger,"Saving NetCDF in: " << fileName);

    // Open the file for read access
    netCDF::NcFile dataFile(fileName, NcFile::replace);

    NcDim particlesDim = dataFile.addDim("particles", this->size());
    NcDim oceanTimeDim = dataFile.addDim("ocean_time", ocean_time);

    NcVar oceanTimeVar = dataFile.addVar("ocean_time", ncDouble, oceanTimeDim);
    oceanTimeVar.putAtt("long_name","time since initialization");
    oceanTimeVar.putAtt("units","seconds since 1968-05-23 00:00:00 GMT");
    oceanTimeVar.putAtt("calendar","gregorian");
    oceanTimeVar.putAtt("field","time, scalar, series");
    oceanTimeVar.putAtt("_CoordinateAxisType","Time");
    oceanTimeVar.putVar(oceanModelAdapter->OceanTime()());


    NcVar idVar = dataFile.addVar("id", ncUint64, particlesDim);
    idVar.putAtt("long_name","id");
    idVar.putAtt("cf_role", "trajectory_id");
    idVar.putVar(id());

    vector<NcDim> particlesOceanTimeDims;
    particlesOceanTimeDims.push_back(particlesDim);
    particlesOceanTimeDims.push_back(oceanTimeDim);

    NcVar latVar = dataFile.addVar("lat", ncDouble, particlesOceanTimeDims);
    latVar.putAtt("long_name","latitude of particle");
    latVar.putAtt("unit","degree_north");
    latVar.putAtt("standard_name","latitude");
    //latVar.putAtt("field","lat, scalar");
    //latVar.putAtt("_coordinateaxistype","lat");
    latVar.putVar(lat());

    NcVar lonVar = dataFile.addVar("lon", ncDouble, particlesOceanTimeDims);
    lonVar.putAtt("long_name","longitude of particle");
    lonVar.putAtt("unit","degree_east");
    lonVar.putAtt("standard_name","longitude");
    //lonVar.putAtt("field","lon, scalar");
    //lonVar.putAtt("_coordinateaxistype","lon");
    lonVar.putVar(lon());

    NcVar depVar = dataFile.addVar("depth", ncDouble, particlesOceanTimeDims);
    depVar.putAtt("long_name","particle depth below sea surface");
    depVar.putAtt("unit","meters");
    depVar.putAtt("standard_name","depth");
    depVar.putAtt("positive","up");
    depVar.putAtt("axis","Z");
    //depVar.putAtt("_coordinateaxistype","depth");
    depVar.putVar(dep());

    NcVar KVar = dataFile.addVar("k", ncDouble, particlesOceanTimeDims);
    KVar.putAtt("long_name","depth fractional index");
    KVar.putVar(k());

    NcVar JVar = dataFile.addVar("j", ncDouble, particlesOceanTimeDims);
    JVar.putAtt("long_name","latitude fractional index");
    JVar.putVar(j());

    NcVar IVar = dataFile.addVar("i", ncDouble, particlesOceanTimeDims);
    IVar.putAtt("long_name","longitude fractional index");
    IVar.putVar(i());

    NcVar HealthVar = dataFile.addVar("health", ncDouble, particlesOceanTimeDims);
    HealthVar.putAtt("long_name","health of the particle");
    HealthVar.putAtt("units","1e-9") ;
    HealthVar.putAtt("coordinates","ocean_time lon lat depth") ;
    HealthVar.putVar(health());

    NcVar AgeVar = dataFile.addVar("age", ncDouble, particlesOceanTimeDims);
    AgeVar.putAtt("long_name","age of the particle");
    AgeVar.putAtt("units","seconds since emission");
    AgeVar.putAtt("coordinates","ocean_time lon lat depth") ;
    AgeVar.putVar(age());

    NcVar TimeVar = dataFile.addVar("time", ncDouble, particlesOceanTimeDims);
    TimeVar.putAtt("long_name","emission time of the particle");
    TimeVar.putAtt("units","seconds since 1968-05-23 00:00:00 GMT");
    TimeVar.putAtt("calendar","gregorian");
    TimeVar.putAtt("coordinates","ocean_time lon lat depth") ;
    //TimeVar.putAtt("field","time, scalar, series");
    //TimeVar.putAtt("_CoordinateAxisType","Time");
    TimeVar.putVar(time());

}

void Particles::loadFromTxt(const string &fileName) {
    std::ifstream infile(fileName);

    unsigned long count=0;
    std::string line;
    while (std::getline(infile, line))
    {
        std::istringstream iss(line);
        if (count > 0)
        {
            unsigned long id=count;
            double k, j, i, health, tpart, emitOceanTime=-1;
            //if (!(iss >> i >> j >> k >> health >> tpart )) { break; } // error

            vector<string> tokens{istream_iterator<string>{iss},
                                  istream_iterator<string>{}};

            // Check if the text file has been produced by fortran WaComM
            if (tokens.size()==5) {
                // Fortran style restart text file
                i=atof(tokens[0].c_str());
                j=atof(tokens[1].c_str());
                k=atof(tokens[2].c_str());
                health=atof(tokens[3].c_str());
                tpart=atof(tokens[4].c_str());
            } else {
                // C++ style restart text file
                id=atoi(tokens[0].c_str());
                i=atof(tokens[1].c_str());
                j=atof(tokens[2].c_str());
                k=atof(tokens[3].c_str());
                health=atof(tokens[4].c_str());
                tpart=atof(tokens[5].c_str());
                emitOceanTime=atof(tokens[6].c_str());
            }
            Particle particle(id, k, j, i, health, tpart, emitOceanTime);
            push_back(particle);
        }
        count++;
    }
}

Particles::~Particles() = default;
