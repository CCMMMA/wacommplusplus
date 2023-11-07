//
// Created by Raffaele Montella on 12/5/20.
//

#include "Particles.hpp"
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <nlohmann/json.hpp>
#include <iomanip>

// for convenience
using json = nlohmann::json;

Particles::Particles() {

    // Logger configuration
    log4cplus::BasicConfigurator basicConfig;
    basicConfig.configure();
    logger = log4cplus::Logger::getInstance(LOG4CPLUS_TEXT("WaComM"));

    LOG4CPLUS_DEBUG(logger, "Empty Particles");
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

void Particles::saveAsJson(const string &fileName, double particleTime, std::shared_ptr<OceanModelAdapter> oceanModelAdapter)
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

void Particles::saveAsNetCDF(const string &fileName, double particleTime, std::shared_ptr<OceanModelAdapter> oceanModelAdapter)
{
    // https://github.com/NOAA-ORR-ERD/nc_particles
    // https://cfconventions.org/Data/cf-conventions/cf-conventions-1.7/build/aphs04.html
    
    size_t particle_time=1;
    size_t ocean_time=1;

    Array1<unsigned long> id(this->size());
    Array2<double> lat(this->size(),particle_time);
    Array2<double> lon(this->size(),particle_time);
    Array2<double> dep(this->size(),particle_time);
    Array2<double> j(this->size(),particle_time);
    Array2<double> i(this->size(),particle_time);
    Array2<double> k(this->size(),particle_time);
    Array2<double> health(this->size(),particle_time);
    Array2<double> age(this->size(),particle_time);
    Array2<double> time(this->size(),particle_time);

    /*
    Array4<float> conc(particle_time,
                        oceanModelAdapter->Depth().Nx(),
                        oceanModelAdapter->Latitude().Nx(),
                        oceanModelAdapter->Longitude().Nx());

    for (int k=0; k<oceanModelAdapter->Depth().Nx(); k++) {
        for (int j=0; k<oceanModelAdapter->Latitude().Nx(); j++) {
            for (int i=0; k<oceanModelAdapter->Longitude().Nx(); i++) {
                conc(0,k,j,i)=0;
            }
        }
    }
    */

    LOG4CPLUS_INFO(logger,"Preparing NetCDF...");

    //double maxDepth=oceanModelAdapter->Depth()[0];
    //double minDepth=oceanModelAdapter->Depth()[oceanModelAdapter->Depth().Nx()-1];
    //double minLat=oceanModelAdapter->Latitude()[0];
    //double maxLat=oceanModelAdapter->Latitude()[oceanModelAdapter->Latitude().Nx()-1];
    //double minLon=oceanModelAdapter->Longitude()[0];
    //double maxLon=oceanModelAdapter->Longitude()[oceanModelAdapter->Longitude().Nx()-1];


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

        //double concK=oceanModelAdapter->Depth().Nx()*depth/(maxDepth-minDepth);
        //double concJ=oceanModelAdapter->Latitude().Nx()*latitude/(maxLat-minLat);
        //double concI=oceanModelAdapter->Longitude().Nx()*longitude/(maxLon-minLon);

        //conc(ocean_time,concK,concJ,concI )++;

        count++;
    }


    LOG4CPLUS_INFO(logger,"Saving NetCDF in: " << fileName);

    // Open the file for read access
    netCDF::NcFile dataFile(fileName, NcFile::replace,NcFile::nc4);

    NcDim particlesDim = dataFile.addDim("particles", this->size());
    NcDim particleTimeDim = dataFile.addDim("particle_time", particle_time);

    NcVar particleTimeVar = dataFile.addVar("particle_time", ncDouble, particleTimeDim);
    particleTimeVar.putAtt("long_name","time since initialization");
    particleTimeVar.putAtt("units","seconds since 1968-05-23 00:00:00 GMT");
    particleTimeVar.putAtt("calendar","gregorian");
    particleTimeVar.putAtt("field","time, scalar, series");
    particleTimeVar.putAtt("_CoordinateAxisType","Time");
    particleTimeVar.putVar(&particleTime);


    NcVar idVar = dataFile.addVar("id", ncUint64, particlesDim);
    idVar.putAtt("long_name","id");
    idVar.putAtt("cf_role", "trajectory_id");
    idVar.putVar(id());

    vector<NcDim> particlesParticleTimeDims;
    particlesParticleTimeDims.push_back(particlesDim);
    particlesParticleTimeDims.push_back(particleTimeDim);

    NcVar latVar = dataFile.addVar("lat", ncDouble, particlesParticleTimeDims);
    latVar.putAtt("long_name","latitude of particle");
    latVar.putAtt("unit","degree_north");
    latVar.putAtt("standard_name","latitude");
    //latVar.putAtt("field","lat, scalar");
    //latVar.putAtt("_coordinateaxistype","lat");
    latVar.putVar(lat());

    NcVar lonVar = dataFile.addVar("lon", ncDouble, particlesParticleTimeDims);
    lonVar.putAtt("long_name","longitude of particle");
    lonVar.putAtt("unit","degree_east");
    lonVar.putAtt("standard_name","longitude");
    //lonVar.putAtt("field","lon, scalar");
    //lonVar.putAtt("_coordinateaxistype","lon");
    lonVar.putVar(lon());

    NcVar depVar = dataFile.addVar("depth", ncDouble, particlesParticleTimeDims);
    depVar.putAtt("long_name","particle depth below sea surface");
    depVar.putAtt("unit","meters");
    depVar.putAtt("standard_name","depth");
    depVar.putAtt("positive","up");
    depVar.putAtt("axis","Z");
    //depVar.putAtt("_coordinateaxistype","depth");
    depVar.putVar(dep());

    NcVar KVar = dataFile.addVar("k", ncDouble, particlesParticleTimeDims);
    KVar.putAtt("long_name","depth fractional index");
    KVar.putVar(k());

    NcVar JVar = dataFile.addVar("j", ncDouble, particlesParticleTimeDims);
    JVar.putAtt("long_name","latitude fractional index");
    JVar.putVar(j());

    NcVar IVar = dataFile.addVar("i", ncDouble, particlesParticleTimeDims);
    IVar.putAtt("long_name","longitude fractional index");
    IVar.putVar(i());

    NcVar HealthVar = dataFile.addVar("health", ncDouble, particlesParticleTimeDims);
    HealthVar.putAtt("long_name","health of the particle");
    HealthVar.putAtt("units","1e-9") ;
    HealthVar.putAtt("coordinates","ocean_time lon lat depth") ;
    HealthVar.putVar(health());

    NcVar AgeVar = dataFile.addVar("age", ncDouble, particlesParticleTimeDims);
    AgeVar.putAtt("long_name","age of the particle");
    AgeVar.putAtt("units","seconds since emission");
    AgeVar.putAtt("coordinates","ocean_time lon lat depth") ;
    AgeVar.putVar(age());

    NcVar TimeVar = dataFile.addVar("time", ncDouble, particlesParticleTimeDims);
    TimeVar.putAtt("long_name","emission time of the particle");
    TimeVar.putAtt("units","seconds since 1968-05-23 00:00:00 GMT");
    TimeVar.putAtt("calendar","gregorian");
    TimeVar.putAtt("coordinates","ocean_time lon lat depth") ;
    //TimeVar.putAtt("field","time, scalar, series");
    //TimeVar.putAtt("_CoordinateAxisType","Time");
    TimeVar.putVar(time());

/*
    NcDim oceanTimeDim = dataFile.addDim("ocean_time", ocean_time);
    NcVar oceanTimeVar = dataFile.addVar("ocean_time", ncDouble, oceanTimeDim);
    oceanTimeVar.putAtt("long_name","time since initialization");
    oceanTimeVar.putAtt("units","seconds since 1968-05-23 00:00:00 GMT");
    oceanTimeVar.putAtt("calendar","gregorian");
    oceanTimeVar.putAtt("field","time, scalar, series");
    oceanTimeVar.putAtt("_CoordinateAxisType","Time");
    oceanTimeVar.putVar(oceanModelAdapter->OceanTime()());

    NcDim latitudeDim = dataFile.addDim("latitude", oceanModelAdapter->Latitude().Nx());
    NcVar laitudeVar = dataFile.addVar("latitude", ncDouble, latitudeDim);
    laitudeVar.putAtt("long_name","latitude o");
    laitudeVar.putAtt("unit","degree_north");
    laitudeVar.putAtt("standard_name","latitude");
    laitudeVar.putAtt("field","latitude, scalar");
    laitudeVar.putAtt("_coordinateaxistype","lat");
    laitudeVar.putVar(oceanModelAdapter->Latitude()());

    NcDim longitudeDim = dataFile.addDim("longitude", oceanModelAdapter->Longitude().Nx());
    NcVar longitudeVar = dataFile.addVar("longitude", ncDouble, longitudeDim);
    longitudeVar.putAtt("long_name","longitude o");
    longitudeVar.putAtt("unit","degree_north");
    longitudeVar.putAtt("standard_name","longitude");
    longitudeVar.putAtt("field","latitude, scalar");
    longitudeVar.putAtt("_coordinateaxistype","lon");
    longitudeVar.putVar(oceanModelAdapter->Longitude()());

    vector<NcDim> oceanTimeDepthLatitudeLongitudeDims;
    oceanTimeDepthLatitudeLongitudeDims.push_back(oceanTimeDim);
    oceanTimeDepthLatitudeLongitudeDims.push_back(depthDim);
    oceanTimeDepthLatitudeLongitudeDims.push_back(latitudeDim);
    oceanTimeDepthLatitudeLongitudeDims.push_back(longitudeDim);

    NcVar concVar = dataFile.addVar("conc", ncFloat, oceanTimeDepthLatitudeLongitudeDims);
    concVar.putAtt("long_name","concentration_of_suspended_matter_in_sea_water");
    concVar.putAtt("units","1");
    concVar.putAtt("coordinates","longitude latitude depth ocean_time");
    concVar.putAtt("field","");
    concVar.putAtt("time","ocean_time");
    concVar.putAtt("_FillValue",ncFloat, 9.99999993e+36);
    concVar.putVar(conc());
    */
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

void Particles::loadFromJson(const string &fileName) {
    std::ifstream infile(fileName);
}

void Particles::loadFromNetCDF(const string &fileName) {
    LOG4CPLUS_INFO(logger,"Reading restart file: " << fileName);
    std::ifstream infile(fileName);
    
    // Open the file for read access
    netCDF::NcFile dataFile(fileName, NcFile::read);

    // Retrieve the variable named "id"
    NcVar varId = dataFile.getVar("id");
    size_t size = varId.getDim(0).getSize();
    Array1<double> Id(size);
    varId.getVar(Id());

    // Retrieve the variable named "i"
    NcVar varI = dataFile.getVar("i");
    Array1<double> I(size);
    varI.getVar(I());

    // Retrieve the variable named "j"
    NcVar varJ = dataFile.getVar("j");
    Array1<double> J(size);
    varJ.getVar(J());

    // Retrieve the variable named "k"
    NcVar varK = dataFile.getVar("k");
    Array1<double> K(size);
    varK.getVar(K());

    // Retrieve the variable named "health"
    NcVar varHealth = dataFile.getVar("health");
    Array1<double> health(size);
    varHealth.getVar(health());

    // Retrieve the variable named "age"
    NcVar varAge = dataFile.getVar("age");
    Array1<double> age(size);
    varAge.getVar(age());

    // Retrieve the variable named "time"
    NcVar varTime = dataFile.getVar("time");
    Array1<double> time(size);
    varTime.getVar(time());

    for (int i=0; i<size; i++) {
        Particle particle(Id[i], K[i], J[i], I[i], health[i], age[i], time[i]);
        push_back(particle);
    }
}

Particles::~Particles() = default;
