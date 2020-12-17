//
// Created by Raffaele Montella on 12/5/20.
//

#include <cfloat>
#include <utility>
#include "Particle.hpp"
#include "Config.hpp"

Particle::Particle(double k, double j, double i,
                   double health, double tpart):
                   k(k), j(j), i(i),
                   health(health),tpart(tpart)
{
    logger = log4cplus::Logger::getInstance(LOG4CPLUS_TEXT("WaComM"));

}

Particle::Particle(double k, double j, double i, double tpart): k(k), j(j), i(i), tpart(tpart){
    health=health0;
    logger = log4cplus::Logger::getInstance(LOG4CPLUS_TEXT("WaComM"));
}


Particle::~Particle() = default;

bool Particle::isAlive() const {
    return health>0;
}



void Particle::move(const std::shared_ptr<Config>& config, int ocean_time_idx,
                    const std::shared_ptr<OceanModelAdapter>& oceanModelAdapter) {

    // Get the integration time (default 30s)
    double dti=config->Dti();

    // Get the time in seconds between two input ocean data (default 3600s, 1h)
    double deltat=config->Deltat();

    // Ask Angelo Riccio (default 86400)
    double tau0=config->Tau0();

    // Probability to survive (default 1.0e-4)
    double survprob=config->Survprob();

    // Reduction Coefficient (default 1)
    double crid=config->ReductionCoefficient();

    // Sedimentation velocity (m-1. default )
    double sv=config->SedimentationVelocity();

    double idet=0,jdet=0,kdet=0;

    // Number of integration intervals
    double iint=deltat/dti;

    // Get the number of the sigma levels
    size_t s_w = oceanModelAdapter->W()->Ny();

    // Get the domain size in south-north number of rows
    size_t eta_rho = oceanModelAdapter->Mask().Nx();

    // Get the domain size in west-east number of columns
    size_t xi_rho = oceanModelAdapter->Mask().Ny();

    // For each integration interval
    for (int t=0;t<iint;t++) {

        // Check of the particle health is less than its probability to survive
        if (health<survprob) {
            // The particle is dead
            health=-1;
            // No reason to continue, exit the integration loop
            break;
        }

        /*
         * IMPORTANT: the k index is negative
         * in the interval [-s_w, 0]
         */

        // Get the integer part and the fraction part of particle k
        auto kI=(int)k; double kF=k-kI;

        // Get the integer part and the fraction part of particle j
        auto jI=(int)j; double jF=j-jI;

        // Get the integer part and the fraction part of particle i
        auto iI=(int)i; double iF=i-iI;

        // Check if the particle is out of the domain
        if (jI<0 || iI<0 || jI>=eta_rho|| iI>=xi_rho) {

            // Set the particle health
            health=-1;

            // no reason to continue,  exit the integration loop
            break;
        }
        // Check if the particle beached
        if (oceanModelAdapter->Mask()(jI,iI)<=0) {

            // Set the particle health
            health=-1;

            // no reason to continue,  exit the integration loop
            break;
        }

        /*
        // Check if the particle jumped outside the water :-)
        if (k>0) {
            // The particle must stay in the water
            k=0;
        }

        // Check if the particle sunk
        if (k>=s_w) {
            // Set the particle health as dead
            health=-1;

            // no reason to continue,  exit the integration loop
            break;
        }
        */

        LOG4CPLUS_DEBUG(logger, this->to_string());

        // The particle is alive!
        // Perform the bilinear interpolation (2D) in order to get
        // the u component of the current field in the particle position.
        float u1=oceanModelAdapter->U()->operator()(ocean_time_idx,   kI,    jI,    iI)    *(1.0-iF)  *(1.0-jF);
        float u2=oceanModelAdapter->U()->operator()(ocean_time_idx,   kI,    jI+1,  iI)    *(1.0-iF)  *     jF;
        float u3=oceanModelAdapter->U()->operator()(ocean_time_idx,   kI,    jI+1,  iI+1)  *     iF   *     jF;
        float u4=oceanModelAdapter->U()->operator()(ocean_time_idx,   kI,    jI  ,  iI+1)  *     iF   *(1.0-jF);

        // The current u component in the particle position
        float uu=u1+u2+u3+u4;

        LOG4CPLUS_DEBUG(logger, "uu:" << uu);

        // Perform the bilinear interpolation (2D) in order to get
        // the v component of the current field in the particle position.
        float v1=oceanModelAdapter->V()->operator()(ocean_time_idx,    kI, jI,     iI)     *(1.0-iF)   *(1.0-jF);
        float v2=oceanModelAdapter->V()->operator()(ocean_time_idx,    kI, jI+1,   iI)     *(1.0-iF)   *     jF;
        float v3=oceanModelAdapter->V()->operator()(ocean_time_idx,    kI, jI+1,   iI+1)   *     iF    *     jF;
        float v4=oceanModelAdapter->V()->operator()(ocean_time_idx,    kI, jI,     iI+1)   *     iF    *(1.0-jF);

        // The current v component in the particle position
        float vv=v1+v2+v3+v4;

        LOG4CPLUS_DEBUG(logger, "vv:" << vv );


        // 0,0,690,533

        // Perform the bilinear interpolation (3D) in order to get
        // the w component of the current field in the particle position.
        float w1=oceanModelAdapter->W()->operator()(ocean_time_idx,    kI,     jI,     iI)*(1.0-iF)  *(1.0-jF)  *(1.0-kF);
        float w2=oceanModelAdapter->W()->operator()(ocean_time_idx,    kI,     jI+1,   iI)*(1.0-iF)  *     jF   *(1.0-kF);
        float w3=oceanModelAdapter->W()->operator()(ocean_time_idx,    kI,     jI+1,   iI+1)*     iF   *     jF   *(1.0-kF);
        float w4=oceanModelAdapter->W()->operator()(ocean_time_idx,    kI,     jI,     iI+1)*     iF   *(1.0-jF)  *(1.0-kF);
        float w5=oceanModelAdapter->W()->operator()(ocean_time_idx,kI-1, jI, iI)*(1.0-iF)*(1.0-jF)*kF;
        float w6=oceanModelAdapter->W()->operator()(ocean_time_idx,kI-1,jI+1,iI)*(1.0-iF)*jF*kF;
        float w7=oceanModelAdapter->W()->operator()(ocean_time_idx,kI-1,jI+1,iI+1)*iF*jF*kF;
        float w8=oceanModelAdapter->W()->operator()(ocean_time_idx, kI-1, jI, iI+1)*iF*(1.0-jF)*kF;

        // The current w component in the particle position
        float ww=w1+w2+w3+w4+w5+w6+w7+w8;

        LOG4CPLUS_DEBUG(logger, "ww:" << ww);



        // Perform the bilinear interpolation (3D) in order to get
        // the akt in the particle position.
        float a1=oceanModelAdapter->AKT()->operator()(ocean_time_idx,  kI,     jI,     iI)*(1.0-iF)   *(1.0-jF)   *(1.0-kF);
        float a2=oceanModelAdapter->AKT()->operator()(ocean_time_idx,  kI,     jI+1,   iI)*(1.0-iF)   *     jF    *(1.0-kF);
        float a3=oceanModelAdapter->AKT()->operator()(ocean_time_idx,  kI,     jI+1,   iI+1)*     iF    *     jF    *(1.0-kF);
        float a4=oceanModelAdapter->AKT()->operator()(ocean_time_idx,  kI,     jI,     iI+1)*     iF    *(1.0-jF)   *(1.0-kF);
        float a5=oceanModelAdapter->AKT()->operator()(ocean_time_idx,  kI-1,   jI,     iI)*(1.0-iF)*(1.0-jF)*kF;
        float a6=oceanModelAdapter->AKT()->operator()(ocean_time_idx,  kI-1,   jI+1,   iI)*(1.0-iF)*jF*kF;
        float a7=oceanModelAdapter->AKT()->operator()(ocean_time_idx,  kI-1,   jI+1,   iI+1)*iF*jF*kF; ;
        float a8=oceanModelAdapter->AKT()->operator()(ocean_time_idx,kI-1,jI,iI+1)*iF*(1.0-jF)*kF;

        // The AKT at the particle position.
        float aa=a1+a2+a3+a4+a5+a6+a7+a8;

        LOG4CPLUS_DEBUG(logger, "aa:" << aa);

        // Evaluate the particle leap due to the current field (deterministic leap).
        double dileap=uu*dti;
        double djleap=vv*dti;
        double dkleap=(sv+ww)*dti;

        // Calculation of sigma profile
        // sigmaPROF=sigma(Ixx,Iyy)*(1-Zdet/(-H(Ixx,Iyy))) ! Here is H
        //double sigmaProf=oceanModelAdapter->sigma()(jI, iI)*(1-kdet/(-oceanModelAdapter->H()(jI,iI)))
        double sigmaprof=3.46*(1+k/s_w);

        // Extract 3 pseudorandom numbers
        double gi=0,gj=0,gk=0;
        if (config->Random()) {
            for (int a = 0; a < 12; a++) {
                gi = gi + gen() - 0.5;
                gj = gj + gen() - 0.5;
                gk = gk + gen() - 0.5;
            }
        }

        // Random leap
        double rileap=gi*sigmaprof;
        double rjleap=gj*sigmaprof;
        double rkleap=gk*aa*crid;

        // Final leap
        double ileap=dileap+rileap;
        double jleap=djleap+rjleap;
        double kleap=dkleap+rkleap;

        LOG4CPLUS_DEBUG(logger, "kleap:" << kleap << " jleap:" << jleap << " ileap:" << ileap );
        double d1,d2,dd,jidist, kdist;

        // Calculate the distance in radiants of latitude between the grid cell where is
        // currently located the particle and the next one.
        d1=(oceanModelAdapter->LatRad()->operator()(jI+1,iI)-oceanModelAdapter->LatRad()->operator()(jI,iI));

        // Calculate the distance in radiants of longitude between the grid cell where is
        // currently located the particle and the next one.
        d2=(oceanModelAdapter->LonRad()->operator()(jI,iI+1)-oceanModelAdapter->LonRad()->operator()(jI,iI));

        // Calculate the grid cell diagonal horizontal size using the Haversine method
        // https://www.movable-type.co.uk/scripts/latlong.html
        dd=pow(sin(0.5*d1),2) +
                pow(sin(0.5*d2),2)*
                cos(oceanModelAdapter->LatRad()->operator()(jI+1,iI))*
                cos(oceanModelAdapter->LatRad()->operator()(jI,iI));
        jidist=2.0*atan2(pow(dd,.5),pow(1.0-dd,.5))*6371.0;

        // Calculate a vertical distance
        double depth=oceanModelAdapter->Depth()->operator()(kI);
        double hcbz=oceanModelAdapter->HCorrectedByZeta(ocean_time_idx,jI,iI);
        //cout << "jI:" << jI << " iI:" << iI << " depth(" << kI <<"):"<<depth<< " hcbz:"<< hcbz << endl;
        kdist=oceanModelAdapter->Depth()->operator()(kI)*oceanModelAdapter->HCorrectedByZeta(ocean_time_idx,jI,iI);
        if ( abs(kleap) > abs(kdist) ) {
            kleap=sign(kdist,kleap);
        }



        // Calculate the new particle j candidate
        jdet=j+0.001*jleap/jidist;

        // Calculate the new particle i candidate
        idet=i+0.001*ileap/jidist;

        // Calculate the new particle k candidate
        kdet=k+kleap/kdist;

        // Reflect if out-of-column
        // Check if the new k have to be limited by the sealfoor
        if ( kdet < (-(int)s_w+2)) {
            // Limit it on the bottom
            kdet=2.0*(-(int)s_w+2)-kdet;
        }

        // Check if the new k have to be limited by the seafloor
        if ( kdet > 0. ) {
            // Limit it on the surface
            kdet=-kdet;
        }

        // Reflect if crossed the coastline

        // Calculate the integer part of the j and i candidates
        int jdetI=(int)(jdet);
        int idetI=(int)(idet);

        // Check if the candidate new particle position is on land (mask=0)
        if ( oceanModelAdapter->Mask()(jdetI,idetI) <= 0.0 ) {
            // Reflect the particle
            if ( idetI < iI ) {
                idet=(double)iI + abs(i-idet);
            } else if ( idetI > iI ) {
                idet=(double)idetI- mod(idet,1.0);
            }
            if ( jdetI < jdet ) {
                jdet=(double)jdetI+ abs(j-jdet);
            } else if ( jdetI > jI ) {
                jdet=(double)jdetI - mod(jdet,1.0);
            }
        }

        // Assign the new particle position
        i=idet;
        j=jdet;
        k=kdet;

        // Update the paticle age
        tpart=tpart+dti;

        // Decay the particle
        health=health0*exp(-tpart/tau0);
    }
}

double Particle::K() const  { return k; }
double Particle::J() const { return j; }
double Particle::I() const { return i; }

// Returns a single pseudorandom number from the uniform distribution over the range [0,1[.
// https://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fNUMBER.html
double Particle::gen() { return Random::get<double>(0.0, 1.0); }

// Returns -1 if a < 0 and 1 if a > 0
double Particle::sgn(double a) { return (a > 0) - (a < 0); }

// Computes the remainder of the division of a by p.
// https://gcc.gnu.org/onlinedocs/gfortran/MOD.html
double Particle::mod(double a, double p) { return a-p*(int)(a/p); }

// Returns the value of a with the sign of b.
// https://gcc.gnu.org/onlinedocs/gfortran/SIGN.html
double Particle::sign(double a, double b) { return abs(a)*sgn(b); }

std::string Particle::to_string() const {
    std::stringstream ss;
    ss << k << " " << j << " " << i << " " << health << " " << tpart;
    return ss.str();
}

double Particle::TPart() const {
    return tpart;
}

double Particle::Health() const {
    return health;
}





