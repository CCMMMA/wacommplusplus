//
// Created by Raffaele Montella on 12/5/20.
//

#include "Wacomm.hpp"

#include <utility>
#include "ROMSAdapter.hpp"

Wacomm::Wacomm(std::shared_ptr<Config> config,
               std::shared_ptr<OceanModelAdapter> oceanModelAdapter,
               std::shared_ptr<Sources> sources,
               std::shared_ptr<Particles> particles) :
               config(std::move(config)),
               //oceanModelAdapter(std::move(oceanModelAdapter)),
               oceanModelAdapter(oceanModelAdapter),
               sources(std::move(sources)),
               particles(std::move(particles)) {

    logger = log4cplus::Logger::getInstance(LOG4CPLUS_TEXT("WaComM"));







}

void Wacomm::run()
{
    size_t ocean_time=oceanModelAdapter->U().Nx();
    size_t s_rho=oceanModelAdapter->U().Ny();
    size_t s_w=oceanModelAdapter->W().Ny();
    size_t eta_rho=oceanModelAdapter->Mask().Nx();
    size_t xi_rho=oceanModelAdapter->Mask().Ny();


    LOG4CPLUS_INFO(logger,"ocean_time:" + std::to_string(ocean_time));
    LOG4CPLUS_INFO(logger,"s_rho:" + std::to_string(s_rho));
    LOG4CPLUS_INFO(logger,"s_w:" + std::to_string(s_w));
    LOG4CPLUS_INFO(logger,"eta_rho:" + std::to_string(eta_rho) + " xi_rho:" + std::to_string(xi_rho));

    conc.Allocate(ocean_time,s_rho,eta_rho,xi_rho);


    for (int ocean_time_idx=0; ocean_time_idx<ocean_time; ocean_time_idx++) {
        LOG4CPLUS_INFO(logger,"Running on:" << ocean_time_idx);
        LOG4CPLUS_INFO(logger,"Sources:" << sources->size());
        for(Source source: *sources) {
            source.emit(particles);
        }

        LOG4CPLUS_INFO(logger,"Particles:" << particles->size());

        int nParticles = particles->size();

        #pragma omp parallel for
        for (int idx=0; idx<nParticles;idx++) {
            particles->at(idx).move(config,
                          ocean_time_idx,
                          oceanModelAdapter);
        }

        // Remove dead particles
        particles->erase(std::remove_if(particles->begin(), particles->end(),
                               [](Particle particle) { return !particle.isAlive(); }), particles->end());

        LOG4CPLUS_INFO(logger,"Particles:" << particles->size());

        LOG4CPLUS_INFO(logger,"Saving restart:" << "");
        particles->save("");

        LOG4CPLUS_INFO(logger,"Evaluate concentration");
        int count=0;
        for (Particle particle: *particles) {
            conc(ocean_time_idx, particle.iK(), particle.iJ(), particle.iI())++;
            count++;
        }
    }

    save(conc);
};
Wacomm::~Wacomm()
{
}

void Wacomm::save(Array4<double> &conc) {

    size_t ocean_time=conc.Nx();
    size_t s_rho=conc.Ny();
    size_t eta_rho=conc.Nz();
    size_t xi_rho=conc.N4();



    std:string fileName="out.nc";

    LOG4CPLUS_INFO(logger,"Saving in: " << fileName);

    // Open the file for read access
    netCDF::NcFile dataFile(fileName, NcFile::replace);

    NcDim oceanTimeDim = dataFile.addDim("ocean_time", ocean_time);
    NcDim sRhoDim = dataFile.addDim("s_rho", s_rho);
    NcDim etaRhoDim = dataFile.addDim("eta_rho", eta_rho);
    NcDim xiRhoDim = dataFile.addDim("eta_xi", xi_rho);

    NcVar oceanTimeVar = dataFile.addVar("ocean_time", ncDouble, oceanTimeDim);
    oceanTimeVar.putVar(oceanModelAdapter->OceanTime()());
    oceanTimeVar.putAtt("long_name","tim e since initializatio");
    oceanTimeVar.putAtt("units","seconds since 1968-05-23 00:00:00 GMT");
    oceanTimeVar.putAtt("calendar","gregorian");
    oceanTimeVar.putAtt("field","time, scalar, series");
    oceanTimeVar.putAtt("_CoordinateAxisType","Time");

    NcVar sRhoVar = dataFile.addVar("s_rho", ncDouble, sRhoDim);
    sRhoVar.putVar(oceanModelAdapter->SRho()());
    sRhoVar.putAtt("long_name","S-coordinate at RHO-points");
    sRhoVar.putAtt("positive","up");
    sRhoVar.putAtt("standard_name","ocean_sigma_coordinate");
    sRhoVar.putAtt("field","s_rho, scalar");
    sRhoVar.putAtt("_CoordinateTransformType","Vertical");
    sRhoVar.putAtt("_CoordinateAxes","s_rho");
    sRhoVar.putAtt("_CoordinateAxisType","GeoZ");
    sRhoVar.putAtt("_CoordinateZisPositive","up");

    vector<NcDim> etaRhoXiRhoDims;
    etaRhoXiRhoDims.push_back(etaRhoDim);
    etaRhoXiRhoDims.push_back(xiRhoDim);
    NcVar latRhoVar = dataFile.addVar("lat_rho", ncDouble, etaRhoXiRhoDims);
    latRhoVar.putAtt("long_name","latitude of rho-points");
    latRhoVar.putAtt("unit","degree_north");
    latRhoVar.putAtt("standard_name","latitude");
    latRhoVar.putAtt("field","lat_rho, scalar");
    latRhoVar.putAtt("_coordinateaxistype","lat");

    NcVar lonRhoVar = dataFile.addVar("lon_rho", ncDouble, etaRhoXiRhoDims);
    lonRhoVar.putAtt("long_name","longitude of rho-points");
    lonRhoVar.putAtt("unit","degree_east");
    lonRhoVar.putAtt("standard_name","longitude");
    lonRhoVar.putAtt("field","lon_rho, scalar");
    lonRhoVar.putAtt("_coordinateaxistype","lon");

    NcVar maskRhoVar = dataFile.addVar("mask_rho", ncDouble, etaRhoXiRhoDims);
    maskRhoVar.putVar(oceanModelAdapter->Mask()());
    maskRhoVar.putAtt("long_name","mask on RHO-points");
    maskRhoVar.putAtt("coordinates","lon_rho lat_rho");
    maskRhoVar.putAtt("units","1");

    NcVar hVar = dataFile.addVar("h", ncDouble, etaRhoXiRhoDims);
    hVar.putVar(oceanModelAdapter->H()());
    hVar.putAtt("long_name","bathymetry a t RHO-point");
    hVar.putAtt("units","meter");
    hVar.putAtt("coordinates","lon_rho lat_rho");
    hVar.putAtt("field","bath, scalar");

    vector<NcDim> oceanTimeEtaRhoXiRhoDims;
    oceanTimeEtaRhoXiRhoDims.push_back(oceanTimeDim);
    oceanTimeEtaRhoXiRhoDims.push_back(etaRhoDim);
    oceanTimeEtaRhoXiRhoDims.push_back(xiRhoDim);

    NcVar hZeta = dataFile.addVar("zeta", ncFloat, oceanTimeEtaRhoXiRhoDims);
    hZeta.putVar(oceanModelAdapter->Zeta()());
    hZeta.putAtt("long_name","free-surface");
    hZeta.putAtt("units","meter");
    hZeta.putAtt("time","ocean_time");
    hZeta.putAtt("coordinates","lon_rho lat_rho ocean_time");
    hZeta.putAtt("field","free-surface, scalar, series");
    //hZeta.putAtt("_FillValue",9.99999993e+36);

    vector<NcDim> oceanTimeSRhoEtaRhoXiRhoDims;
    oceanTimeSRhoEtaRhoXiRhoDims.push_back(oceanTimeDim);
    oceanTimeSRhoEtaRhoXiRhoDims.push_back(sRhoDim);
    oceanTimeSRhoEtaRhoXiRhoDims.push_back(etaRhoDim);
    oceanTimeSRhoEtaRhoXiRhoDims.push_back(xiRhoDim);
    NcVar concVar = dataFile.addVar("conc", ncFloat, oceanTimeSRhoEtaRhoXiRhoDims);
    concVar.putVar(conc());
    concVar.putAtt("long_name","concentration_of_suspended_matter_in_sea_water");
    concVar.putAtt("units","1");
    concVar.putAtt("coordinates","lon_rho lat_rho s_rho");
    concVar.putAtt("field","");
    concVar.putAtt("time","ocean_time");

}

