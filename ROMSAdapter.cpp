//
// Created by Raffaele Montella on 07/12/20.
//

#include "ROMSAdapter.hpp"

ROMSAdapter::ROMSAdapter(string &fileName): fileName(fileName) {
    logger = log4cplus::Logger::getInstance(LOG4CPLUS_TEXT("WaComM"));
}

void ROMSAdapter::process()
{
    LOG4CPLUS_INFO(logger,"Loading:"+fileName);

    // Open the file for read access
    netCDF::NcFile dataFile(fileName, NcFile::read);

    // Retrieve the variable named "mask_rho"
    NcVar varMaskRho=dataFile.getVar("mask_rho");
    size_t eta_rho = varMaskRho.getDim(0).getSize();
    size_t xi_rho = varMaskRho.getDim(1).getSize();
    Array2<double>mask_rho(eta_rho,xi_rho);
    varMaskRho.getVar(mask_rho());

    // Retrieve the variable named "mask_u"
    NcVar varMaskU=dataFile.getVar("mask_u");
    size_t eta_u = varMaskU.getDim(0).getSize();
    size_t xi_u = varMaskU.getDim(1).getSize();
    Array2<double> mask_u(eta_u,xi_u);
    varMaskU.getVar(mask_u());

    // Retrieve the variable named "mask_v"
    NcVar varMaskV=dataFile.getVar("mask_v");
    size_t eta_v = varMaskV.getDim(0).getSize();
    size_t xi_v = varMaskV.getDim(1).getSize();
    Array2<double> mask_v(eta_v,xi_v);
    varMaskV.getVar(mask_v());

    // Retrieve the variable named "lat_rho"
    NcVar varLatRho=dataFile.getVar("lat_rho");
    Array2<double> lat_rho(eta_rho,xi_rho);
    varLatRho.getVar(lat_rho());

    // Retrieve the variable named "lon_rho"
    NcVar varLonRho=dataFile.getVar("lon_rho");
    Array2<double> lon_rho(eta_rho,xi_rho);
    varLonRho.getVar(lon_rho());

    // Retrieve the variable named "lat_v"
    NcVar varLatV=dataFile.getVar("lat_v");
    Array2<double> lat_v(eta_v,xi_v);
    varLatV.getVar(lat_v());

    // Retrieve the variable named "lon_u"
    NcVar varLonU=dataFile.getVar("lon_u");
    Array2<double> lon_u(eta_u,xi_u);
    varLonU.getVar(lon_u());

    // Retrieve the variable named "h"
    NcVar varH=dataFile.getVar("h");
    Array2<double> h(eta_rho,xi_rho);
    varH.getVar(h());

    // Retrieve the variable named "u"
    NcVar varU=dataFile.getVar("u");
    size_t ocean_time = varU.getDim(0).getSize();
    size_t s_rho = varU.getDim(1).getSize();
    Array4<float> u(ocean_time, s_rho, eta_u, xi_u,0,-(int)s_rho+1,0,0);
    varU.getVar(u());

    // Retrieve the variable named "v"
    NcVar varV=dataFile.getVar("v");
    Array4<float> v(ocean_time, s_rho, eta_v, xi_v,0,-(int)s_rho+1,0,0);
    varV.getVar(v());

    // Retrieve the variable named "w"
    NcVar varW=dataFile.getVar("w");
    size_t s_w = varW.getDim(1).getSize();
    Array4<float> w(ocean_time,s_w,eta_rho,xi_rho,0,-(int)s_w+1,0,0);
    varW.getVar(w());

    // Retrieve the variable named "AKt"
    NcVar varAKt=dataFile.getVar("AKt");
    Array4<float> akt(ocean_time,s_w,eta_rho,xi_rho,0,-(int)s_w+1,0,0);
    varAKt.getVar(akt());


    // Retrieve the variable named "s_w"
    NcVar varSW=dataFile.getVar("s_w");
    varSW.getDim(0).getSize();
    Array1<double> sW(s_w,-(int)s_w+1);
    varSW.getVar(sW());

    // Retrieve the variable named "s_rho"
    NcVar varSRho=dataFile.getVar("s_rho");
    varSRho.getDim(0).getSize();
    Array1<double> sRho(s_rho, -(int)s_rho+1);
    varSRho.getVar(sRho());

    // Retrieve the variable named "ocean_time"
    NcVar varOceanTime=dataFile.getVar("ocean_time");
    varOceanTime.getDim(0).getSize();
    Array1<double> oceanTime(ocean_time);
    varOceanTime.getVar(oceanTime());

    // Retrieve the variable named "zeta"
    NcVar varZeta=dataFile.getVar("zeta");
    Array3<float> zeta(ocean_time,eta_rho,xi_rho);
    varZeta.getVar(zeta());

    LOG4CPLUS_INFO(logger,"Loaded!");

    this->OceanTime().Allocate(ocean_time);
    this->SRho().Allocate(s_rho, -(int)s_rho+1);
    this->SW().Allocate(s_w, -(int)s_w+1);
    this->Depth()->Allocate(s_w,-(int)s_w+2);
    this->Mask().Allocate(eta_rho,xi_rho);
    this->Lon().Allocate(eta_rho,xi_rho);
    this->Lat().Allocate(eta_rho,xi_rho);
    this->LonRad()->Allocate(eta_rho,xi_rho);
    this->LatRad()->Allocate(eta_rho,xi_rho);
    this->H()->Allocate(eta_rho,xi_rho);
    this->Zeta()->Allocate(ocean_time,eta_rho,xi_rho);
    this->U()->Allocate(ocean_time,s_rho,eta_rho,xi_rho,0,-(int)s_rho+1,0,0);
    this->V()->Allocate(ocean_time,s_rho,eta_rho,xi_rho,0,-(int)s_rho+1,0,0);
    this->W()->Allocate(ocean_time,s_w,eta_rho,xi_rho,0,-(int)s_w+1,0,0);
    this->AKT()->Allocate(ocean_time,s_w,eta_rho,xi_rho,0,-(int)s_w+1,0,0);

    LOG4CPLUS_INFO(logger,"Copying 1D ...");
/*
    for(int t=0; t<ocean_time;t++) {
        this->OceanTime()(t)=oceanTime(t);
    }*/
    memcpy( this->OceanTime(), oceanTime(), ocean_time*sizeof(double) );
    /*
    for(int k=-(int)s_rho+1; k<=0;k++) {
        this->SRho()(k)=sRho(k);
    }*/
    memcpy( this->SRho(), sRho(), s_rho*sizeof(double) );
    /*
    for(int k=-(int)s_w+1; k<=0;k++) {
        this->SW()(k)=sW(k);
    }*/
    memcpy( this->SW(), sW(), s_w*sizeof(double) );

    for(int k=-(int)s_w+2; k<=0;k++) {
        this->Depth()->operator()(k)=sW(k)-sW(k-1);
    }

    LOG4CPLUS_INFO(logger,"Copying 2D...");

    memcpy( this->Mask(), mask_rho(), eta_rho*xi_rho*sizeof(double) );
    memcpy( this->Lat(), lat_rho(), eta_rho*xi_rho*sizeof(double) );
    memcpy( this->Lon(), lon_rho(), eta_rho*xi_rho*sizeof(double) );
    memcpy( this->H()->operator double *(), h(), eta_rho*xi_rho*sizeof(double) );

    LOG4CPLUS_INFO(logger,"Copying 3D...");
    memcpy( this->Zeta()->operator float *(), zeta(), ocean_time*eta_rho*xi_rho*sizeof(float) );

    LOG4CPLUS_INFO(logger,"Convert lon & lat in radiants...");
    // lon_u, lat_v in radiants
    #pragma omp for collapse(2)
    for ( int j=0; j<eta_rho;j++) {
        for (int i=0;i<xi_rho;i++) {
            this->LonRad()->operator()(j,i)=0.0174533*lon_rho(j,i);
            this->LatRad()->operator()(j,i)=0.0174533*lat_rho(j,i);
        }
    }

    LOG4CPLUS_INFO(logger,"Interpolation 2D...");
    uv2rho(mask_rho, mask_u, mask_v, u,v);

    LOG4CPLUS_INFO(logger,"Interpolation 3D...");
    wakt2rho(mask_rho,mask_u, mask_v, w, akt);


    LOG4CPLUS_INFO(logger,"...done!");
}

void ROMSAdapter::uv2rho(Array2<double>& mask_rho, Array2<double>& mask_u, Array2<double>& mask_v,
                         Array4<float>& u, Array4<float>& v) {
    float uw1, uw2, vw1, vw2;

    LOG4CPLUS_INFO(logger,"uv2rho");

    size_t ocean_time = u.Nx();
    size_t s_rho = u.Ny();
    size_t eta_v = mask_v.Nx();
    size_t xi_v = mask_v.Ny();
    size_t eta_u = mask_v.Nx();
    size_t xi_u = mask_u.Ny();
    size_t eta_rho = mask_rho.Nx();
    size_t xi_rho = mask_rho.Ny();

    #pragma omp for collapse(4)
    for (int t=0; t < ocean_time; t++)
    {
        for (int k=(-(int)s_rho+1); k <=0; k++)
        {
            for (int j=0; j < eta_v; j++)
            {
                for (int i=0; i< xi_u; i++)
                {

                    if ( j>=0 && i>=0 && j < eta_rho && i < xi_rho && mask_rho(j,i) > 0.0 )
                    {
                        if ( j>=0 && i>=0 && j<eta_u && i <xi_u && mask_u(j,i) > 0.0 )
                        {
                            uw1=u(t,k,j,i);
                        } else {
                            uw1=0.0;
                        }
                        if ( j>0 && i>=0 && j<eta_u && i <xi_u && mask_u(j-1,i) > 0.0 )
                        {
                            uw2=u(t,k,j-1,i);
                        } else
                        {
                            uw2=0.0;
                        }
                        if ( j>=0 && i>=0 && j<eta_v && i <xi_v && mask_v(j,i) > 0.0 )
                        {
                            vw1=v(t,k,j,i);
                        } else
                        {
                            vw1=0.0;
                        }
                        if ( j>=0 && i>0 && j<eta_v && i <xi_v && mask_v(j,i-1) > 0.0 )
                        {
                            vw2=v(t,k,j,i-1);
                        } else
                        {
                            vw2=0.0;
                        }

                        this->U()->operator()(t,k,j,i)=0.5*(uw1+uw2);
                        this->V()->operator()(t,k,j,i)=0.5*(vw1+vw2);
                    } else {
                        this->U()->operator()(t,k,j,i)=0.0;
                        this->V()->operator()(t,k,j,i)=0.0;
                    }
                }
            }
        }
    }
}

void ROMSAdapter::wakt2rho(Array2<double>& mask_rho, Array2<double>& mask_u, Array2<double>& mask_v,
                           Array4<float>& w, Array4<float>& akt ) {
    float ww0, ww1, ww2, ww3, aktw0, aktw1, aktw2, aktw3;

    LOG4CPLUS_INFO(logger,"wakt2rho");

    size_t ocean_time = w.Nx();
    size_t s_w = w.Ny();
    size_t eta_v = mask_v.Nx();
    size_t xi_u = mask_u.Ny();
    size_t eta_rho = mask_rho.Nx();
    size_t xi_rho = mask_rho.Ny();


    #pragma omp for collapse(4)
    for (int t=0; t < ocean_time; t++) {
        for (int k=(-(int)s_w+1); k <= 0; k++) {
            for (int j=0; j < eta_v; j++) {
                for (int i=0; i< xi_u; i++) {

                    if ( j>=0 && i>=0 && j<eta_rho && i <xi_rho && mask_rho(j,i) > 0.0 )
                    {
                        /*
                        ww0 = w(t,k,j,i);
                        aktw0 = akt(t,k,j,i);

                        if (ww0 >= 9.99e36) ww0=0;
                        if (aktw0 >= 9.99e36) aktw0=0;
                        */
                        if ( j>=0 && i>0 && j<eta_rho && i <xi_rho && mask_rho(j,i-1) > 0.0 )
                        {
                            ww1=w(t, k,j,i-1);
                            aktw1=akt(t, k,j,i-1);
                            /*
                            if (ww1 >= 9.99e36) ww1=0;
                            if (aktw1 >= 9.99e36) aktw1=0;
                             */
                        } else {
                            ww1=0.0;
                            aktw1=0.0;
                        }

                        if ( j>0 && i>=0 && j<eta_rho && i <xi_rho && mask_rho(j-1,i) > 0.0 ) {
                            ww2=w(t, k,j-1,i);
                            aktw2=akt(t, k,j-1,i);
                            /*
                            if (ww2 >= 9.99e36) ww2=0;
                            if (aktw2 >= 9.99e36) aktw2=0;
                             */
                        } else {
                            ww2=0.0;
                            aktw2=0.0;
                        }

                        if ( j>0 && i>0 && j<eta_rho && i <xi_rho && mask_rho(j-1,i-1) > 0.0 ) {
                            ww3=w(t, k, j-1,i-1);
                            aktw3=akt(t, k, j-1,i-1);
                            /*
                            if (ww3 >= 9.99e36) ww3=0;
                            if (aktw3 >= 9.99e36) aktw3=0;
                             */
                        } else {
                            ww3=0.0;
                            aktw3=0.0;
                        }

                        this->W()->operator()(t,k,j,i)=0.25*(w(t,k,j,i)+ww1+ww2+ww3);
                        this->AKT()->operator()(t,k,j,i)=0.25*(akt(t,k,j,i)+aktw1+aktw2+aktw3);

                    } else {

                        this->W()->operator()(t,k,j,i)=0.0;
                        this->AKT()->operator()(t,k,j,i)=0.0;
                    }
                }
            }
        }
    }
}



ROMSAdapter::~ROMSAdapter() = default;




