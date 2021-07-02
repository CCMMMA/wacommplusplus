//
// Created by Raffaele Montella on 29/01/21.
//

#include "WacommAdapter.hpp"

WacommAdapter::WacommAdapter(string &fileName): fileName(fileName) {
    logger = log4cplus::Logger::getInstance(LOG4CPLUS_TEXT("WaComM"));
}

WacommAdapter::~WacommAdapter() = default;

void WacommAdapter::process()  {
    LOG4CPLUS_DEBUG(logger, "Native file:" + fileName);

    // Open the file for read access
    netCDF::NcFile dataFile(fileName, NcFile::read,NcFile::nc4);

    NcDim oceanTimeDim = dataFile.getDim("ocean_time");
    NcDim sRhoDim = dataFile.getDim("s_rho");
    NcDim sWDim = dataFile.getDim("s_w");
    NcDim etaRhoDim = dataFile.getDim("eta_rho");
    NcDim xiRhoDim = dataFile.getDim("eta_xi");

    size_t ocean_time=oceanTimeDim.getSize();
    size_t s_rho=sRhoDim.getSize();
    size_t s_w=sWDim.getSize();
    size_t eta_rho=etaRhoDim.getSize();
    size_t xi_rho=xiRhoDim.getSize();


    this->OceanTime().Allocate(ocean_time);
    this->SRho().Allocate(s_rho, -(int)s_rho+1);
    this->SW().Allocate(s_w, -(int)s_w+1);
    this->DepthIntervals().Allocate(s_w,-(int)s_w+2);
    this->Mask().Allocate(eta_rho,xi_rho);
    this->Lon().Allocate(eta_rho,xi_rho);
    this->Lat().Allocate(eta_rho,xi_rho);
    this->LonRad().Allocate(eta_rho,xi_rho);
    this->LatRad().Allocate(eta_rho,xi_rho);
    this->H().Allocate(eta_rho,xi_rho);
    this->Zeta().Allocate(ocean_time,eta_rho,xi_rho);
    this->U().Allocate(ocean_time,s_rho,eta_rho,xi_rho,0,-(int)s_rho+1,0,0);
    this->V().Allocate(ocean_time,s_rho,eta_rho,xi_rho,0,-(int)s_rho+1,0,0);
    this->W().Allocate(ocean_time,s_w,eta_rho,xi_rho,0,-(int)s_w+1,0,0);
    this->AKT().Allocate(ocean_time,s_w,eta_rho,xi_rho,0,-(int)s_w+1,0,0);

    NcVar oceanTimeVar = dataFile.getVar("ocean_time");
    oceanTimeVar.getVar(this->OceanTime()());

    //cout << "OceanModelAdapter::loadFromNetCDF : " << _data.oceanTime(0) << endl;


    NcVar sRhoVar = dataFile.getVar("s_rho");
    sRhoVar.getVar(this->SRho()());

    NcVar sWVar = dataFile.getVar("s_w");
    sWVar.getVar(this->SW()());

    NcVar latRhoVar = dataFile.getVar("lat_rho");
    latRhoVar.getVar(this->Lat()());

    NcVar lonRhoVar = dataFile.getVar("lon_rho");
    lonRhoVar.getVar(this->Lon()());

    NcVar maskRhoVar = dataFile.getVar("mask_rho");
    maskRhoVar.getVar(this->Mask()());

    NcVar hVar = dataFile.getVar("h");
    hVar.getVar(this->H()());

    //cout << "OceanModelAdapter::loadFromNetCDF _data.h(650,550):" << _data.h(650,550) << endl;

    NcVar varZeta = dataFile.getVar("zeta");
    varZeta.getVar(this->Zeta()());

    NcVar uVar = dataFile.getVar("u");
    uVar.getVar(this->U()());

    NcVar vVar = dataFile.getVar("v");
    vVar.getVar(this->V()());

    NcVar wVar = dataFile.getVar("w");
    wVar.getVar(this->W()());

    NcVar aktVar = dataFile.getVar("akt");
    aktVar.getVar(this->AKT()());

    for(int k=-(int)s_w+2; k<=0;k++) {
        this->DepthIntervals().operator()(k)=this->SW()(k)-this->SW()(k-1);
    }

    // lon_u, lat_v in radiants
    #pragma omp parallel for collapse(2) default(none) shared(eta_rho, xi_rho )
    for ( int j=0; j<eta_rho;j++) {
        for (int i=0;i<xi_rho;i++) {
            this->LonRad().operator()(j,i)=0.0174533*this->Lon()(j,i);
            this->LatRad().operator()(j,i)=0.0174533*this->Lat()(j,i);
        }
    }
}