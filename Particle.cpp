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

Particle::Particle(double k, double j, double i): k(k), j(j), i(i){
    health=health0;
    tpart=0;
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
    size_t s_w = oceanModelAdapter->W().Ny();

    // Get the domain size in south-north number of rows
    size_t eta_rho = oceanModelAdapter->Mask().Nx();

    // Get the domain size in west-east number of columns
    size_t xi_rho = oceanModelAdapter->Mask().Ny();

    // For each integration interval
    for (int t=0;t<iint;t++) {

        // Check of the particle health is less than its probability to survive
        if (health<survprob) {
            // The particle is dead
            pstatus=0;
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
        auto jI=(int)i; double jF=j-jI;

        // Get the integer part and the fraction part of particle i
        auto iI=(int)i; double iF=i-iI;

        // Checl if the particle is out of the domain
        if (jI<0 || iI<0 || jI>=eta_rho|| iI>=xi_rho) {

            // Set the particle health
            health=-1;

            // The particle is dead
            pstatus=0;

            // no reason to continue,  exit the integration loop
            break;
        }
        // Check if the particle beached
        if (oceanModelAdapter->Mask()(jI,iI)<1) {

            // Set the particle health
            health=-1;

            // The particle is dead
            pstatus=0;

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
            // Set the particle health
            health=-13;

            // The particle is dead
            pstatus=0;

            // no reason to continue,  exit the integration loop
            break;
        }
*/
        LOG4CPLUS_DEBUG(logger, this->to_string());



        // The particle is alive!
        // Perform the bilinear interpolation (2D) in order to get
        // the u component of the current field in the particle position.
        double u1=oceanModelAdapter->U()(ocean_time_idx,   kI,    jI,    iI)    *(1.0-iF)  *(1.0-jF);
        double u2=oceanModelAdapter->U()(ocean_time_idx,   kI,    jI+1,  iI)    *(1.0-iF)  *     jF;
        double u3=oceanModelAdapter->U()(ocean_time_idx,   kI,    jI+1,  iI+1)  *     iF   *     jF;
        double u4=oceanModelAdapter->U()(ocean_time_idx,   kI,    jI  ,  iI+1)  *     iF   *(1.0-jF);

        // The current u component in the particle position
        double uu=u1+u2+u3+u4;

        LOG4CPLUS_DEBUG(logger, "uu:" << uu);

        // Perform the bilinear interpolation (2D) in order to get
        // the v component of the current field in the particle position.
        double v1=oceanModelAdapter->V()(ocean_time_idx,    kI, jI,     iI)     *(1.0-iF)   *(1.0-jF);
        double v2=oceanModelAdapter->V()(ocean_time_idx,    kI, jI+1,   iI)     *(1.0-iF)   *     jF;
        double v3=oceanModelAdapter->V()(ocean_time_idx,    kI, jI+1,   iI+1)   *     iF    *     jF;
        double v4=oceanModelAdapter->V()(ocean_time_idx,    kI, jI,     iI+1)   *     iF    *(1.0-jF);

        // The current v component in the particle position
        double vv=v1+v2+v3+v4;

        LOG4CPLUS_DEBUG(logger, "vv:" << vv );



        // Perform the bilinear interpolation (3D) in order to get
        // the w component of the current field in the particle position.
        double w1=oceanModelAdapter->W()(ocean_time_idx,    kI,     jI,     iI)     *(1.0-iF)  *(1.0-jF)  *(1.0-kF);
        double w2=oceanModelAdapter->W()(ocean_time_idx,    kI,     jI+1,   iI)     *(1.0-iF)  *     jF   *(1.0-kF);
        double w3=oceanModelAdapter->W()(ocean_time_idx,    kI,     jI+1,   iI+1)   *     iF   *     jF   *(1.0-kF);
        double w4=oceanModelAdapter->W()(ocean_time_idx,    kI,     jI,     iI+1)   *     iF   *(1.0-jF)  *(1.0-kF);
        double w5=oceanModelAdapter->W()(ocean_time_idx,    kI-1,   jI,     iI)     *(1.0-iF)  *(1.0-jF)  *     kF;
        double w6=oceanModelAdapter->W()(ocean_time_idx,    kI-1,   jI+1,   iI)     *(1.0-iF)  *     jF   *     kF;
        double w7=oceanModelAdapter->W()(ocean_time_idx,    kI-1,   jI+1,   iI+1)   *     iF   *     jF   *     kF;
        double w8=oceanModelAdapter->W()(ocean_time_idx,    kI-1,   jI,     iI+1)   *     iF   *(1.0-jF)  *     kF;

        // The current w component in the particle position
        double ww=w1+w2+w3+w4+w5+w6+w7+w8;

        LOG4CPLUS_DEBUG(logger, "ww:" << ww);



        // Perform the bilinear interpolation (3D) in order to get
        // the akt in the particle position.
        double a1=oceanModelAdapter->AKT()(ocean_time_idx,  kI,     jI,     iI)     *(1.0-iF)   *(1.0-jF)   *(1.0-kF);
        double a2=oceanModelAdapter->AKT()(ocean_time_idx,  kI,     jI+1,   iI)     *(1.0-iF)   *     jF    *(1.0-kF);
        double a3=oceanModelAdapter->AKT()(ocean_time_idx,  kI,     jI+1,   iI+1)   *     iF    *     jF    *(1.0-kF);
        double a4=oceanModelAdapter->AKT()(ocean_time_idx,  kI,     jI,     iI+1)   *     iF    *(1.0-jF)   *(1.0-kF);
        double a5=oceanModelAdapter->AKT()(ocean_time_idx,  kI-1,   jI,     iI)     *(1.0-iF)   *(1.0-jF)   *kF;
        double a6=oceanModelAdapter->AKT()(ocean_time_idx,  kI-1,   jI+1,   iI)     *(1.0-iF)   *     jF    *kF;
        double a7=oceanModelAdapter->AKT()(ocean_time_idx,  kI-1,   jI+1,   iI+1)   *iF         *     jF    *kF;
        double a8=oceanModelAdapter->AKT()(ocean_time_idx,  kI-1,   jI,     iI+1)   *iF         *(1.0-jF)   *kF;


        // The AKT at the particle position.
        double aa=a1+a2+a3+a4+a5+a6+a7+a8;

        LOG4CPLUS_DEBUG(logger, "aa:" << aa);

        // Evaluate the particle leap due to the current field (deterministic leap).
        double dileap=uu*dti;
        double djleap=vv*dti;
        double dkleap=(sv+ww)*dti;

        // Extract 3 pseudorandom numbers
        double gi=0,gj=0,gk=0;
        for (int a=0;a<12;a++) {
            gi=gi+gen()-0.5;
            gj=gj+gen()-0.5;
            gk=gk+gen()-0.5;
        }

        // Calculation of sigma control
        // sigmaPROF=sigma(Ixx,Iyy)*(1-Zdet/(-H(Ixx,Iyy))) ! Here is H
        //double sigmaProf=oceanModelAdapter->sigma()(jI, iI)*(1-kdet/(-oceanModelAdapter->H()(jI,iI)))
        double sigmaprof=3.46*(1+kdet/s_w);

        // Random leap
        double rileap=gi*sigmaprof;
        double rjleap=gj*sigmaprof;
        double rkleap=gk*aa*crid;

        // Final leap
        double ileap=dileap+rileap;
        double jleap=djleap+rjleap;
        double kleap=dkleap+rkleap;

        LOG4CPLUS_DEBUG(logger, "kleap:" << kleap << " jleap:" << jleap << " ileap:" << ileap );
        double d1,d2,dist;

        d1=(oceanModelAdapter->LatRad()(jI+1,iI)-oceanModelAdapter->LatRad()(jI,iI));
        d2=(oceanModelAdapter->LonRad()(jI,iI+1)-oceanModelAdapter->LonRad()(jI,iI));
        d1=pow(sin(0.5*d1),2) +
                pow(sin(0.5*d2),2)*
                cos(oceanModelAdapter->LatRad()(jI+1,iI))*
                cos(oceanModelAdapter->LatRad()(jI,iI));
        dist=2.0*atan2(pow(d1,.5),pow(1.0-d1,.5))*6371.0;
        idet=i+0.001*ileap/dist;

        d1=(oceanModelAdapter->LatRad()(jI+1,iI)-oceanModelAdapter->LatRad()(jI,iI));
        d2=(oceanModelAdapter->LonRad()(jI,iI+1)-oceanModelAdapter->LonRad()(jI,iI));
        d1=pow(sin(0.5*d1),2) +
                pow(sin(0.5*d2),2)*
                cos(oceanModelAdapter->LatRad()(jI+1,iI))*
                cos(oceanModelAdapter->LatRad()(jI,iI));
        dist=2.0*atan2(sqrt(d1),pow(1.0-d1,.5))*6371.0;
        jdet=j+0.001*jleap/dist;

        //printf("depth:%f, zeta:%f\n",depth[kI],zeta[jI][iI]);
        dist=oceanModelAdapter->Depth()(kI)*oceanModelAdapter->Zeta()(ocean_time_idx,jI,iI);
        if (dist<=0) dist=1e-4;
        //printf("dist:%f, kleap:%f\n",dist,kleap);
        if ( abs(kleap) > abs(dist) ) {
            kleap=sign(dist,kleap);
        }
        kdet=k+kleap/dist;

        // Check if the new k have to be limited by the sealfoor
        if ( kdet >= s_w ) {
            // Limit it on the bottom
            kdet=s_w-1;
        }

        // Check if the new k have to be limited by the surface
        if ( kdet < 0. ) {
            // Limit it on the surface
            kdet=0;
        }

        //printf("kdet:%f jdet:%f idet:%f\n",kdet,jdet,idet);

        int jdetI=(int)(jdet);
        int idetI=(int)(idet);
        if ( oceanModelAdapter->Mask()(jdetI,idetI) <= 0.0 ) {
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

        // Add the integration time to the time
        time=time+dti;

        // Decay the particle
        health=health0*exp(-time/tau0);
    }
}

double Particle::K() const  { return k; }
double Particle::J() const { return j; }
double Particle::I() const { return i; }

int Particle::KasInt() const  { return (int)(round(k));}
int Particle::JasInt() const { return (int)(round(j));}
int Particle::IasInt() const { return (int)(round(i));}

/* The random generator should be initialized only one time */
double Particle::gen()
{
    std::random_device randomDevice;
    auto randomGenerator=std::mt19937 (randomDevice());
    std::uniform_real_distribution<double> randomDistribution=std::uniform_real_distribution<double>(0, std::nextafter(1, DBL_MAX));
    return randomDistribution(randomGenerator);
}

double Particle::sgn(double a) { return (a > 0) - (a < 0); }
double Particle::mod(double a, double p) { return a-p*(int)(a/p); }
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





