//
// Created by Raffaele Montella on 10/12/20.
//

#include "OceanModelAdapter.hpp"
#include <math.h>

OceanModelAdapter::OceanModelAdapter() {
    logger = log4cplus::Logger::getInstance(LOG4CPLUS_TEXT("WaComM"));
}

Array1<double> &OceanModelAdapter::OceanTime() { return _data.oceanTime; }
Array1<double> &OceanModelAdapter::SRho() { return _data.sRho; }
Array1<double> &OceanModelAdapter::SW() { return _data.sW; }
Array2<double> &OceanModelAdapter::Mask() { return _data.mask; }
Array2<double> &OceanModelAdapter::Lon() { return _data.lon; }
Array2<double> &OceanModelAdapter::Lat() { return _data.lat; }
Array2<double> &OceanModelAdapter::LonRad() { return _data.lonRad; }
Array2<double> &OceanModelAdapter::LatRad() { return _data.latRad; }
Array1<double> &OceanModelAdapter::Depth() { return _data.depth; }
Array2<double> &OceanModelAdapter::H() { return _data.h; }
Array3<float> &OceanModelAdapter::Zeta() { return _data.zeta; }
Array4<float> &OceanModelAdapter::U() { return _data.u; }
Array4<float> &OceanModelAdapter::V() { return _data.v; }
Array4<float> &OceanModelAdapter::W() { return _data.w; }
Array4<float> &OceanModelAdapter::AKT() { return _data.akt; }



void OceanModelAdapter::saveAsNetCDF(std::string &fileName) {

    size_t ocean_time=_data.oceanTime.Nx();
    size_t s_rho=_data.sRho.Nx();
    size_t s_w=_data.sW.Nx();
    size_t eta_rho=_data.mask.Nx();
    size_t xi_rho=_data.mask.Ny();

    LOG4CPLUS_DEBUG(logger,"Saving in: " << fileName);

    // Open the file for read access
    netCDF::NcFile dataFile(fileName, NcFile::replace,NcFile::nc4);

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
    oceanTimeVar.putVar(_data.oceanTime());

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
    sRhoVar.putVar(_data.sRho());

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
    sWVar.putVar(_data.sW());

    vector<NcDim> etaRhoXiRhoDims;
    etaRhoXiRhoDims.push_back(etaRhoDim);
    etaRhoXiRhoDims.push_back(xiRhoDim);
    NcVar latRhoVar = dataFile.addVar("lat_rho", ncDouble, etaRhoXiRhoDims);
    latRhoVar.putAtt("long_name","latitude of rho-points");
    latRhoVar.putAtt("unit","degree_north");
    latRhoVar.putAtt("standard_name","latitude");
    latRhoVar.putAtt("field","lat_rho, scalar");
    latRhoVar.putAtt("_coordinateaxistype","lat");
    latRhoVar.putVar(_data.lat());

    NcVar lonRhoVar = dataFile.addVar("lon_rho", ncDouble, etaRhoXiRhoDims);
    lonRhoVar.putAtt("long_name","longitude of rho-points");
    lonRhoVar.putAtt("unit","degree_east");
    lonRhoVar.putAtt("standard_name","longitude");
    lonRhoVar.putAtt("field","lon_rho, scalar");
    lonRhoVar.putAtt("_coordinateaxistype","lon");
    lonRhoVar.putVar(_data.lon());

    NcVar maskRhoVar = dataFile.addVar("mask_rho", ncDouble, etaRhoXiRhoDims);
    maskRhoVar.putAtt("long_name","mask on RHO-points");
    maskRhoVar.putAtt("coordinates","lon_rho lat_rho");
    maskRhoVar.putAtt("units","1");
    maskRhoVar.putVar(_data.mask());

    NcVar hVar = dataFile.addVar("h", ncDouble, etaRhoXiRhoDims);
    hVar.putAtt("long_name","bathymetry at RHO-point");
    hVar.putAtt("units","meter");
    hVar.putAtt("coordinates","lon_rho lat_rho");
    hVar.putAtt("field","bath, scalar");
    hVar.putVar(_data.h());

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
    varZeta.putVar(_data.zeta());

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
    uVar.putVar(_data.u());

    NcVar vVar = dataFile.addVar("v", ncFloat, oceanTimeSRhoEtaRhoXiRhoDims);
    vVar.putAtt("long_name","v-momentum component at RHO-points");
    vVar.putAtt("units","meter second-1");
    vVar.putAtt("grid","grid");
    vVar.putAtt("loction","face");
    vVar.putAtt("coordinates","lon_rho lat_rho s_rho ocean_time");
    vVar.putAtt("field","v-velocity, scalar, series");
    vVar.putAtt("time","ocean_time");
    vVar.putAtt("_FillValue",ncFloat, 9.99999993e+36);
    vVar.putVar(_data.v());

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
    wVar.putVar(_data.w());

    NcVar aktVar = dataFile.addVar("akt", ncFloat, oceanTimeSWEtaRhoXiRhoDims);
    aktVar.putAtt("long_name","temperature vertical diffusion coefficient");
    aktVar.putAtt("units","meter2 second-1");
    aktVar.putAtt("grid","grid");
    aktVar.putAtt("loction","face");
    aktVar.putAtt("coordinates","lon_rho lat_rho s_w ocean_time");
    aktVar.putAtt("field","AKt, scalar, series");
    aktVar.putAtt("time","ocean_time");
    aktVar.putVar(_data.akt());

}

oceanmodel_data *OceanModelAdapter::dataptr() {
    return &_data;
}

void OceanModelAdapter::kji2deplatlon(double k, double j, double i, double &dep, double &lat, double &lon) {
    // Get the integer part and the fraction part of particle k
    auto kI=(int)k; double kF=k-kI;

    // Get the integer part and the fraction part of particle j
    auto jI=(int)j; double jF=j-jI;

    // Get the integer part and the fraction part of particle i
    auto iI=(int)i; double iF=i-iI;

    // Check if the source must be skipped
    if (jI < 0 || iI < 0 || jI>=_data.mask.Nx() || iI>=_data.mask.Ny()) {
        lon = 1e37;
        lat = 1e37;
        dep = 1e37;
        return;
    }

    // Perform the bilinear interpolation (2D) in order to get
    // the lat at the source position.
    double lon1=_data.lon(jI,    iI)    *(1.0-iF)  *(1.0-jF);
    double lon2=_data.lon(jI+1,  iI)    *(1.0-iF)  *     jF;
    double lon3=_data.lon(jI+1,  iI+1)  *     iF   *     jF;
    double lon4=_data.lon(jI  ,  iI+1)  *     iF   *(1.0-jF);

    // The current lon (longitude) at the source position
    lon=lon1+lon2+lon3+lon4;


    // Perform the bilinear interpolation (2D) in order to get
    // the lat at the source position.
    double lat1=_data.lat(jI,    iI)    *(1.0-iF)  *(1.0-jF);
    double lat2=_data.lat(jI+1,  iI)    *(1.0-iF)  *     jF;
    double lat3=_data.lat(jI+1,  iI+1)  *     iF   *     jF;
    double lat4=_data.lat(jI  ,  iI+1)  *     iF   *(1.0-jF);

    // The current lat (latitude) at the source position
    lat=lat1+lat2+lat3+lat4;

    // Perform the bilinear interpolation (2D) in order to get
    // the h (depth) at the particle position.
    double h1=_data.h(jI,    iI)    *(1.0-iF)  *(1.0-jF);
    double h2=_data.h(jI+1,  iI)    *(1.0-iF)  *     jF;
    double h3=_data.h(jI+1,  iI+1)  *     iF   *     jF;
    double h4=_data.h(jI  ,  iI+1)  *     iF   *(1.0-jF);

    // The current h (depth) at the particle position
    double h=h1+h2+h3+h4;

    double aDep = abs(
            abs(h * _data.sW(kI-1)) -
            abs(h * _data.sW(kI))
    );
    dep=-(h*abs(_data.sW(kI))+abs(kF*aDep));
}

void OceanModelAdapter::deplatlon2kji(double dep, double lat, double lon, double &k, double &j, double &i) {
    int minK, minJ, minI;
    double d, d1, d2, dd, minD=1e37;

    double latRad=0.0174533*lat;
    double lonRad=0.0174533*lon;

    size_t eta_rho = _data.mask.Nx();
    size_t xi_rho = _data.mask.Ny();
    size_t s_w = _data.w.Ny();

    for (int j=0; j<eta_rho; j++) {
        for (int i=0; i<eta_rho; i++) {
            // Calculate the distance in radiants of latitude between the grid cell where is
            // currently located the particle and the next one.
            d1=(latRad-_data.latRad(j,i));

            // Calculate the distance in radiants of longitude between the grid cell where is
            // currently located the particle and the next one.
            d2=(lonRad-_data.lonRad(j,i));

            // Calculate the grid cell diagonal horizontal size using the Haversine method
            // https://www.movable-type.co.uk/scripts/latlong.html
            dd=pow(sin(0.5*d1),2) +
               pow(sin(0.5*d2),2)*
               cos(latRad)*
               cos(_data.latRad(j,i));
            d=2.0*atan2(pow(dd,.5),pow(1.0-dd,.5))*6371.0;

            if (d<minD) {
                minD=d;
                minJ=j;
                minI=i;
            }
        }
    }

    double dLat=latRad-_data.latRad(minJ, minI);
    double dLon=lonRad-_data.lonRad(minJ, minI);
    int otherJ=minJ+sgn(dLat);
    int otherI=minI+sgn(dLon);
    if (dLat!=0) {
        double aLat = abs(_data.latRad(minJ, minI)-_data.latRad(otherJ, otherI));
        double jF=abs(dLat)/aLat;
        j=min(minJ,otherJ)+jF;
    } else {
        j=minJ;
    }

    if (dLon!=0) {
        double aLon = abs(_data.lonRad(minJ, minI)-_data.lonRad(otherJ, otherI));
        double iF=abs(dLon)/aLon;
        i=min(minI,otherI)+iF;
    } else {
        i=minI;
    }

    // Get the integer part and the fraction part of particle j
    auto jI=(int)j; double jF=j-jI;

    // Get the integer part and the fraction part of particle i
    auto iI=(int)i; double iF=i-iI;

    // Convert dep to positive down
    dep = abs(dep);

    // Perform the bilinear interpolation (2D) in order to get
    // the h (depth) at the particle position.
    double h1=_data.h(jI,    iI)    *(1.0-iF)  *(1.0-jF);
    double h2=_data.h(jI+1,  iI)    *(1.0-iF)  *     jF;
    double h3=_data.h(jI+1,  iI+1)  *     iF   *     jF;
    double h4=_data.h(jI  ,  iI+1)  *     iF   *(1.0-jF);

    // The current h (depth) at the particle position
    double h=h1+h2+h3+h4;

    // Check if the depth is deeper than h
    if (dep>h) {
        // The position is about at the bottom
        k=-(int) s_w + 1;
    } else
        // CHeck if it is on the surface
        if (dep==0) {
            k=-1;
    } else {
            minD = 1e37;
            double hs;
            for (int k = (-(int) s_w + 1); k <= 0; k++) {
                hs = h * abs(_data.sW(k));
                d = abs(hs - dep);
                if (d < minD) {
                    minD = d;
                    minK = k;
                }
            }
            double dDep = dep - h * abs(_data.sW(minK));
            if (dDep != 0) {
                int otherK = minK + sgn(dDep);
                double aDep = abs(
                        abs(h * _data.sW(minK)) -
                        abs(h * _data.sW(otherK))
                );
                double kF = abs(dDep) / aDep;
                k = min(minK, otherK) + kF;
            } else {
                k = minK;
            }
        }

}

// Returns -1 if a < 0 and 1 if a > 0
double OceanModelAdapter::sgn(double a) { return (a > 0) - (a < 0); }




