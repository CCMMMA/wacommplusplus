//
// Created by Raffaele Montella on 10/12/20.
//

#include "OceanModelAdapter.hpp"

OceanModelAdapter::OceanModelAdapter() {
    logger = log4cplus::Logger::getInstance(LOG4CPLUS_TEXT("WaComM"));
}

Array1<double> &OceanModelAdapter::OceanTime() { return oceanTime; }
Array1<double> &OceanModelAdapter::SRho() { return sRho; }
Array1<double> &OceanModelAdapter::SW() { return sW; }
Array2<double> &OceanModelAdapter::Mask() { return mask; }
Array2<double> &OceanModelAdapter::Lon() { return lon; }
Array2<double> &OceanModelAdapter::Lat() { return lat; }
Array2<double> *OceanModelAdapter::LonRad() { return &lonRad; }
Array2<double> *OceanModelAdapter::LatRad() { return &latRad; }
Array1<double> *OceanModelAdapter::Depth() { return &depth; }
Array2<double> *OceanModelAdapter::H() { return &h; }
Array3<float> *OceanModelAdapter::Zeta() { return &zeta; }
Array4<float> *OceanModelAdapter::U() { return &u; }
Array4<float> *OceanModelAdapter::V() { return &v; }
Array4<float> *OceanModelAdapter::W() { return &w; }
Array4<float> *OceanModelAdapter::AKT() { return &akt; }

double OceanModelAdapter::HCorrectedByZeta(int ocean_time, int j, int i) {
    return h(j,i)+zeta(ocean_time,j,i);
}

void OceanModelAdapter::saveAsNetCDF(std::string &fileName) {

    size_t ocean_time=oceanTime.Nx();
    size_t s_rho=sRho.Nx();
    size_t s_w=sW.Nx();
    size_t eta_rho=mask.Nx();
    size_t xi_rho=mask.Ny();

    LOG4CPLUS_INFO(logger,"Saving in: " << fileName);

    // Open the file for read access
    netCDF::NcFile dataFile(fileName, NcFile::replace);

    NcDim oceanTimeDim = dataFile.addDim("ocean_time", ocean_time);
    NcDim sRhoDim = dataFile.addDim("s_rho", s_rho);
    NcDim sWDim = dataFile.addDim("s_w", s_w);
    NcDim etaRhoDim = dataFile.addDim("eta_rho", eta_rho);
    NcDim xiRhoDim = dataFile.addDim("eta_xi", xi_rho);

    NcVar oceanTimeVar = dataFile.addVar("ocean_time", ncDouble, oceanTimeDim);
    oceanTimeVar.putAtt("long_name","time since initialization");
    oceanTimeVar.putAtt("units","seconds since 1968-05-23 00:00:00 GMT");
    oceanTimeVar.putAtt("calendar","gregorian");
    oceanTimeVar.putAtt("field","time, scalar, series");
    oceanTimeVar.putAtt("_CoordinateAxisType","Time");
    oceanTimeVar.putVar(oceanTime());

    NcVar sRhoVar = dataFile.addVar("s_rho", ncDouble, sRhoDim);
    sRhoVar.putAtt("long_name","S-coordinate at RHO-points");
    sRhoVar.putAtt("valid_min",ncDouble, -1.0);
    sRhoVar.putAtt("valid_max",ncDouble, 0.0);
    sRhoVar.putAtt("positive","up");
    sRhoVar.putAtt("standard_name","ocean_sigma_coordinates");
    sRhoVar.putAtt("field","s_rho, scalar");
    sRhoVar.putAtt("_CoordinateTransformType","Vertical");
    sRhoVar.putAtt("_CoordinateAxes","s_rho");
    sRhoVar.putAtt("_CoordinateAxisType","GeoZ");
    sRhoVar.putAtt("_CoordinateZisPositive","up");
    sRhoVar.putVar(sRho());

    NcVar sWVar = dataFile.addVar("s_w", ncDouble, sWDim);
    sWVar.putAtt("long_name","S-coordinate at W-points");
    sWVar.putAtt("valid_min",ncDouble, -1.0);
    sWVar.putAtt("valid_max",ncDouble, 0.0);
    sWVar.putAtt("positive","up");
    sWVar.putAtt("standard_name","ocean_sigma_coordinates");
    sWVar.putAtt("field","s_w, scalar");
    sWVar.putAtt("_CoordinateTransformType","Vertical");
    sWVar.putAtt("_CoordinateAxes","s_w");
    sWVar.putAtt("_CoordinateAxisType","GeoZ");
    sWVar.putAtt("_CoordinateZisPositive","up");
    sWVar.putVar(sW());

    vector<NcDim> etaRhoXiRhoDims;
    etaRhoXiRhoDims.push_back(etaRhoDim);
    etaRhoXiRhoDims.push_back(xiRhoDim);
    NcVar latRhoVar = dataFile.addVar("lat_rho", ncDouble, etaRhoXiRhoDims);
    latRhoVar.putAtt("long_name","latitude of rho-points");
    latRhoVar.putAtt("unit","degree_north");
    latRhoVar.putAtt("standard_name","latitude");
    latRhoVar.putAtt("field","lat_rho, scalar");
    latRhoVar.putAtt("_coordinateaxistype","lat");
    latRhoVar.putVar(lat());

    NcVar lonRhoVar = dataFile.addVar("lon_rho", ncDouble, etaRhoXiRhoDims);
    lonRhoVar.putAtt("long_name","longitude of rho-points");
    lonRhoVar.putAtt("unit","degree_east");
    lonRhoVar.putAtt("standard_name","longitude");
    lonRhoVar.putAtt("field","lon_rho, scalar");
    lonRhoVar.putAtt("_coordinateaxistype","lon");
    lonRhoVar.putVar(lon());

    NcVar maskRhoVar = dataFile.addVar("mask_rho", ncDouble, etaRhoXiRhoDims);
    maskRhoVar.putAtt("long_name","mask on RHO-points");
    maskRhoVar.putAtt("coordinates","lon_rho lat_rho");
    maskRhoVar.putAtt("units","1");
    maskRhoVar.putVar(mask());

    NcVar hVar = dataFile.addVar("h", ncDouble, etaRhoXiRhoDims);
    hVar.putAtt("long_name","bathymetry at RHO-point");
    hVar.putAtt("units","meter");
    hVar.putAtt("coordinates","lon_rho lat_rho");
    hVar.putAtt("field","bath, scalar");
    hVar.putVar(h());

    vector<NcDim> oceanTimeEtaRhoXiRhoDims;
    oceanTimeEtaRhoXiRhoDims.push_back(oceanTimeDim);
    oceanTimeEtaRhoXiRhoDims.push_back(etaRhoDim);
    oceanTimeEtaRhoXiRhoDims.push_back(xiRhoDim);

    NcVar varZeta = dataFile.addVar("zeta", ncFloat, oceanTimeEtaRhoXiRhoDims);
    varZeta.putAtt("long_name","free-surface");
    varZeta.putAtt("units","meter");
    varZeta.putAtt("time","ocean_time");
    varZeta.putAtt("coordinates","lon_rho lat_rho ocean_time");
    varZeta.putAtt("field","free-surface, scalar, series");
    varZeta.putAtt("_FillValue",ncFloat, 9.99999993e+36);
    varZeta.putVar(zeta());

    vector<NcDim> oceanTimeSRhoEtaRhoXiRhoDims;
    oceanTimeSRhoEtaRhoXiRhoDims.push_back(oceanTimeDim);
    oceanTimeSRhoEtaRhoXiRhoDims.push_back(sRhoDim);
    oceanTimeSRhoEtaRhoXiRhoDims.push_back(etaRhoDim);
    oceanTimeSRhoEtaRhoXiRhoDims.push_back(xiRhoDim);

    NcVar uVar = dataFile.addVar("u", ncFloat, oceanTimeSRhoEtaRhoXiRhoDims);
    uVar.putAtt("long_name","u-momentum component at RHO-points");
    uVar.putAtt("units","meter second-1");
    uVar.putAtt("grid","grid");
    uVar.putAtt("loction","face");
    uVar.putAtt("coordinates","lon_rho lat_rho s_rho ocean_time");
    uVar.putAtt("field","u-velocity, scalar, series");
    uVar.putAtt("time","ocean_time");
    uVar.putAtt("_FillValue",ncFloat, 9.99999993e+36);
    uVar.putVar(u());

    NcVar vVar = dataFile.addVar("v", ncFloat, oceanTimeSRhoEtaRhoXiRhoDims);
    vVar.putAtt("long_name","v-momentum component at RHO-points");
    vVar.putAtt("units","meter second-1");
    vVar.putAtt("grid","grid");
    vVar.putAtt("loction","face");
    vVar.putAtt("coordinates","lon_rho lat_rho s_rho ocean_time");
    vVar.putAtt("field","v-velocity, scalar, series");
    vVar.putAtt("time","ocean_time");
    vVar.putAtt("_FillValue",ncFloat, 9.99999993e+36);
    vVar.putVar(v());

    vector<NcDim> oceanTimeSWEtaRhoXiRhoDims;
    oceanTimeSWEtaRhoXiRhoDims.push_back(oceanTimeDim);
    oceanTimeSWEtaRhoXiRhoDims.push_back(sWDim);
    oceanTimeSWEtaRhoXiRhoDims.push_back(etaRhoDim);
    oceanTimeSWEtaRhoXiRhoDims.push_back(xiRhoDim);

    NcVar wVar = dataFile.addVar("w", ncFloat, oceanTimeSWEtaRhoXiRhoDims);
    wVar.putAtt("long_name","vertical momentum component");
    wVar.putAtt("units","meter second-1");
    wVar.putAtt("grid","grid");
    wVar.putAtt("loction","face");
    wVar.putAtt("coordinates","lon_rho lat_rho s_w ocean_time");
    wVar.putAtt("field","w-velocity, scalar, series");
    wVar.putAtt("time","ocean_time");
    wVar.putAtt("_FillValue",ncFloat, 9.99999993e+36);
    wVar.putVar(w());

    NcVar aktVar = dataFile.addVar("akt", ncFloat, oceanTimeSWEtaRhoXiRhoDims);
    aktVar.putAtt("long_name","temperature vertical diffusion coefficient");
    aktVar.putAtt("units","meter2 second-1");
    aktVar.putAtt("grid","grid");
    aktVar.putAtt("loction","face");
    aktVar.putAtt("coordinates","lon_rho lat_rho s_w ocean_time");
    aktVar.putAtt("field","AKt, scalar, series");
    aktVar.putAtt("time","ocean_time");
    aktVar.putVar(akt());

}

void OceanModelAdapter::process() {
    LOG4CPLUS_ERROR(logger,"OceanModelAdapter::process must be implemented in a concrete adapter!");
    exit(-1);
}


