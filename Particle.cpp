//
// Created by Raffaele Montella on 12/5/20.
//

#include <cfloat>
#include <utility>
#include "Particle.hpp"
#include "Config.hpp"


Particle::Particle(unsigned long id, double k, double j, double i,
                   double health, double age, double time)
{
#ifdef DEBUG
    logger = log4cplus::Logger::getInstance(LOG4CPLUS_TEXT("WaComM"));
#endif

    _data.id=id;
    _data.k=k;
    _data.j=j;
    _data.i=i;
    _data.health=health;
    _data.age=age;
    _data.time=time;

}

Particle::Particle(unsigned long id, double k, double j, double i, double time) {

#ifdef DEBUG
    logger = log4cplus::Logger::getInstance(LOG4CPLUS_TEXT("WaComM"));
#endif
    _data.id=id;
    _data.k=k;
    _data.j=j;
    _data.i=i;
    _data.health=health0;
    _data.age=0;
    _data.time=time;
}

Particle::Particle(particle_data data) {
#ifdef DEBUG
    logger = log4cplus::Logger::getInstance(LOG4CPLUS_TEXT("WaComM"));
#endif
    _data = data;
}

Particle::~Particle() = default;

particle_data Particle::data() {
    return _data;
}

void Particle::data(particle_data data) {
    _data = data;
}



bool Particle::isAlive() const {
    return _data.health>0;
}



void Particle::move(config_data *configData, int ocean_time_idx, Array1<double> &oceanTime, Array2<double> &mask,
                    Array2<double> &lonRad, Array2<double> &latRad, Array1<double> &depth,
                    Array2<double> &h, Array3<float> &zeta, Array4<float> &u, Array4<float> &v, Array4<float> &w,
                    Array4<float> &akt) {

    unsigned seed = Random::get<unsigned>(0, UINT32_MAX);

    particle_data localParticleData;
    memcpy(&localParticleData, &_data, sizeof(particle_data));

    // Get the random flag
    bool random=configData->random;

    // Get the integration time (default 30s)
    double dti=configData->dti;

    // Get the time in seconds between two input ocean data (default 3600s, 1h)
    double deltat=configData->deltat;

    // Ask Angelo Riccio (default 86400)
    double tau0=configData->tau0;

    // Probability to survive (default 1.0e-4)
    double survprob=configData->survprob;

    // Reduction Coefficient (default 1)
    double crid=configData->crid;

    // Sedimentation velocity (m-1. default )
    double sv=configData->sv;

    // Shore limit (positive depth: default 0.25)
    double shoreLimit=configData->shoreLimit;

    double idet=0,jdet=0,kdet=0;

    // Number of integration intervals
    double iint=deltat/dti;

    // Get the number of the sigma levels
    size_t s_w = w.Ny();

    // Get the domain size in south-north number of rows
    size_t eta_rho = mask.Nx();

    // Get the domain size in west-east number of columns
    size_t xi_rho = mask.Ny();

    // For each integration interval
    for (int t=0;t<iint;t++) {

        // Check if the paticle is not yet active
        if (localParticleData.time>(oceanTime(ocean_time_idx)+(t*dti))) {
            // The particle is not already active (already emitted, but not active)
            break;
        }

        // Check of the particle health is less than its probability to survive
        if (localParticleData.health<survprob) {
            // The particle is dead
            localParticleData.health=-1;
            // No reason to continue, exit the integration loop
            break;
        }

        /*
         * IMPORTANT: the k index is negative
         * in the interval [-s_w, 0]
         */

        // Get the integer part and the fraction part of particle k
        auto kI=(int)localParticleData.k; double kF=localParticleData.k-kI;

        // Get the integer part and the fraction part of particle j
        auto jI=(int)localParticleData.j; double jF=localParticleData.j-jI;

        // Get the integer part and the fraction part of particle i
        auto iI=(int)localParticleData.i; double iF=localParticleData.i-iI;

        // Check if the particle is out of the domain
        if (jI<0 || iI<0 || jI>=eta_rho|| iI>=xi_rho) {

            // Set the particle health
            localParticleData.health=-1;

            // no reason to continue,  exit the integration loop
            break;
        }

        /*
        // Check if the particle beached
        if (mask(jI,iI)<=0) {

            // Set the particle health
            localParticleData.health=-1;

            // no reason to continue,  exit the integration loop
            break;
        }
        */
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

#ifdef DEBUG
        LOG4CPLUS_DEBUG(logger, this->to_string());
#endif

        // The particle is alive!

        // Perform the bilinear interpolation (2D) in order to get
        // the zeta at the particle position.
        float z1 = zeta(ocean_time_idx, jI, iI) * (1.0 - iF) * (1.0 - jF);
        float z2 = zeta(ocean_time_idx, jI + 1, iI) * (1.0 - iF) * jF;
        float z3 = zeta(ocean_time_idx, jI + 1, iI + 1) * iF * jF;
        float z4 = zeta(ocean_time_idx, jI, iI + 1) * iF * (1.0 - jF);

        // The current zeta at the particle position
        float zz = z1 + z2 + z3 + z4;

        // Perform the bilinear interpolation (2D) in order to get
        // the h (depth) at the particle position.
        double h1=h(jI,    iI)    *(1.0-iF)  *(1.0-jF);
        double h2=h(jI+1,  iI)    *(1.0-iF)  *     jF;
        double h3=h(jI+1,  iI+1)  *     iF   *     jF;
        double h4=h(jI  ,  iI+1)  *     iF   *(1.0-jF);

        // The current h (depth) at the particle position
        double hh=h1+h2+h3+h4;

        // Calulate the corrected depth
        double hc=hh+zz;

        // Check if the partiche is still floating...
        if (hc>shoreLimit) {
            // The particle is still floating

            // Perform the bilinear interpolation (2D) in order to get
            // the u component of the current field in the particle position.
            float u1 = u(ocean_time_idx, kI, jI, iI) * (1.0 - iF) * (1.0 - jF);
            float u2 = u(ocean_time_idx, kI, jI + 1, iI) * (1.0 - iF) * jF;
            float u3 = u(ocean_time_idx, kI, jI + 1, iI + 1) * iF * jF;
            float u4 = u(ocean_time_idx, kI, jI, iI + 1) * iF * (1.0 - jF);

            // The current u component in the particle position
            float uu = u1 + u2 + u3 + u4;

#ifdef DEBUG
            LOG4CPLUS_DEBUG(logger, "uu:" << uu);
#endif

            // Perform the bilinear interpolation (2D) in order to get
            // the v component of the current field in the particle position.
            float v1 = v(ocean_time_idx, kI, jI, iI) * (1.0 - iF) * (1.0 - jF);
            float v2 = v(ocean_time_idx, kI, jI + 1, iI) * (1.0 - iF) * jF;
            float v3 = v(ocean_time_idx, kI, jI + 1, iI + 1) * iF * jF;
            float v4 = v(ocean_time_idx, kI, jI, iI + 1) * iF * (1.0 - jF);

            // The current v component in the particle position
            float vv = v1 + v2 + v3 + v4;

#ifdef DEBUG
            LOG4CPLUS_DEBUG(logger, "vv:" << vv );
#endif

            // 0,0,690,533

            // Perform the bilinear interpolation (3D) in order to get
            // the w component of the current field in the particle position.
            float w1 = w(ocean_time_idx, kI, jI, iI) * (1.0 - iF) * (1.0 - jF) * (1.0 - kF);
            float w2 = w(ocean_time_idx, kI, jI + 1, iI) * (1.0 - iF) * jF * (1.0 - kF);
            float w3 = w(ocean_time_idx, kI, jI + 1, iI + 1) * iF * jF * (1.0 - kF);
            float w4 = w(ocean_time_idx, kI, jI, iI + 1) * iF * (1.0 - jF) * (1.0 - kF);
            float w5 = w(ocean_time_idx, kI - 1, jI, iI) * (1.0 - iF) * (1.0 - jF) * kF;
            float w6 = w(ocean_time_idx, kI - 1, jI + 1, iI) * (1.0 - iF) * jF * kF;
            float w7 = w(ocean_time_idx, kI - 1, jI + 1, iI + 1) * iF * jF * kF;
            float w8 = w(ocean_time_idx, kI - 1, jI, iI + 1) * iF * (1.0 - jF) * kF;

            // The current w component in the particle position
            float ww = w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8;

#ifdef DEBUG
            LOG4CPLUS_DEBUG(logger, "ww:" << ww);
#endif


            // Perform the bilinear interpolation (3D) in order to get
            // the akt in the particle position.
            float a1 = akt(ocean_time_idx, kI, jI, iI) * (1.0 - iF) * (1.0 - jF) * (1.0 - kF);
            float a2 = akt(ocean_time_idx, kI, jI + 1, iI) * (1.0 - iF) * jF * (1.0 - kF);
            float a3 = akt(ocean_time_idx, kI, jI + 1, iI + 1) * iF * jF * (1.0 - kF);
            float a4 = akt(ocean_time_idx, kI, jI, iI + 1) * iF * (1.0 - jF) * (1.0 - kF);
            float a5 = akt(ocean_time_idx, kI - 1, jI, iI) * (1.0 - iF) * (1.0 - jF) * kF;
            float a6 = akt(ocean_time_idx, kI - 1, jI + 1, iI) * (1.0 - iF) * jF * kF;
            float a7 = akt(ocean_time_idx, kI - 1, jI + 1, iI + 1) * iF * jF * kF;;
            float a8 = akt(ocean_time_idx, kI - 1, jI, iI + 1) * iF * (1.0 - jF) * kF;

            // The AKT at the particle position.
            float aa = a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8;

#ifdef DEBUG
            LOG4CPLUS_DEBUG(logger, "aa:" << aa);
#endif

            // Evaluate the particle leap due to the current field (deterministic leap).
            double dileap = uu * dti;
            double djleap = vv * dti;
            double dkleap = (sv + ww) * dti;

            // Calculation of sigma profile
            // sigmaPROF=sigma(Ixx,Iyy)*(1-Zdet/(-H(Ixx,Iyy))) ! Here is H
            //double sigmaProf=oceanModelAdapter->sigma()(jI, iI)*(1-kdet/(-oceanModelAdapter->H()(jI,iI)))
            double sigmaprof = 3.46 * (1 + localParticleData.k / s_w);

            // Extract 3 pseudorandom numbers
            double gi = 0, gj = 0, gk = 0;

            if (random) {
                for (int a = 0; a < 12; a++) {
                    gi = gi + gen(&seed) - 0.5;
                    gj = gj + gen(&seed) - 0.5;
                    gk = gk + gen(&seed) - 0.5;
                }
            }


            // Random leap
            double rileap = gi * sigmaprof;
            double rjleap = gj * sigmaprof;
            double rkleap = gk * aa * crid;

            // Final leap
            double ileap = dileap + rileap;
            double jleap = djleap + rjleap;
            double kleap = dkleap + rkleap;

#ifdef DEBUG
            LOG4CPLUS_DEBUG(logger, "kleap:" << kleap << " jleap:" << jleap << " ileap:" << ileap );
#endif

            double d1, d2, dd, jidist, kdist;

            // Calculate the distance in radiants of latitude between the grid cell where is
            // currently located the particle and the next one.
            d1 = (latRad(jI + 1, iI) - latRad(jI, iI));

            // Calculate the distance in radiants of longitude between the grid cell where is
            // currently located the particle and the next one.
            d2 = (lonRad(jI, iI + 1) - lonRad(jI, iI));

            // Calculate the grid cell diagonal horizontal size using the Haversine method
            // https://www.movable-type.co.uk/scripts/latlong.html
            dd = pow(sin(0.5 * d1), 2) +
                 pow(sin(0.5 * d2), 2) *
                 cos(latRad(jI + 1, iI)) *
                 cos(latRad(jI, iI));
            jidist = 2.0 * atan2(pow(dd, .5), pow(1.0 - dd, .5)) * 6371.0;



            //cout << "jI:" << jI << " iI:" << iI << " depth(" << kI <<"):"<<depth<< " hcbz:"<< hcbz << endl;
            kdist = depth(kI) * (hh + zz);
            if (abs(kleap) > abs(kdist)) {
                kleap = sign(kdist, kleap);
            }

            // Calculate the new particle j candidate
            jdet = localParticleData.j + 0.001 * jleap / jidist;

            // Calculate the new particle i candidate
            idet = localParticleData.i + 0.001 * ileap / jidist;

            // Calculate the new particle k candidate
            kdet = localParticleData.k + kleap / kdist;

            /*
            // Reflect if out-of-column
            // Check if the new k have to be limited by the sealfoor
            if (kdet < (-(int) s_w + 2)) {
                // Limit it on the bottom
                kdet = 2.0 * (-(int) s_w + 2) - kdet;
            }

            // Check if the new k have to be limited by the seafloor
            if (kdet > 0.) {
                // Limit it on the surface
                kdet = -kdet;
            }
            */

            // Check if the new k have to be limited by the sealfoor
            if (kdet < (-(int) s_w + 2)) {
                // Limit it on the bottom
                kdet = (-(int) s_w + 2) ;
            }

            // Check if the new k have to be limited by the seafloor
            if (kdet > 0.) {
                // Limit it on the surface
                kdet = 0;
            }

            // Reflect if crossed the coastline

            // Calculate the integer part of the j and i candidates
            int jdetI = (int) (jdet);
            int idetI = (int) (idet);

            // Check if the candidate position is within the domain
            if (jdetI >= 0 && idetI >= 0 && jdetI < eta_rho && idetI < xi_rho) {
                /*
                // Check if the candidate new particle position is on land (mask=0)
                if (mask(jdetI, idetI) <= 0.0) {
                    // Reflect the particle
                    if (idetI < iI) {
                        idet = (double) iI + abs(localParticleData.i - idet);
                    } else if (idetI > iI) {
                        idet = (double) idetI - mod(idet, 1.0);
                    }
                    if (jdetI < jdet) {
                        jdet = (double) jdetI + abs(localParticleData.j - jdet);
                    } else if (jdetI > jI) {
                        jdet = (double) jdetI - mod(jdet, 1.0);
                    }
                }
                */

                // Assign the new particle position
                localParticleData.i = idet;
                localParticleData.j = jdet;
                localParticleData.k = kdet;
            }


        }

        // Update the paticle age
        localParticleData.age=localParticleData.age+dti;

        // Decay the particle
        localParticleData.health=health0*exp(-localParticleData.age/tau0);
    }
    memcpy(&_data, &localParticleData, sizeof(particle_data));
}

double Particle::K() const  { return _data.k; }
double Particle::J() const { return _data.j; }
double Particle::I() const { return _data.i; }

// Returns a single pseudorandom number from the uniform distribution over the range [0,1[.
// https://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fNUMBER.html
//double Particle::gen() { return Random::get<double>(0.0, 1.0); }
double Particle::gen(unsigned int *seed) { return (double) rand_r(seed)/RAND_MAX; }

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
    ss << _data.id << " " << _data.k << " " << _data.j << " " << _data.i << " " << _data.health << " " << _data.age << " " << _data.time;
    return ss.str();
}

double Particle::Age() const {
    return _data.age;
}

double Particle::Health() const {
    return _data.health;
}

double Particle::Time() const {
    return _data.time;
}

unsigned long Particle::Id() const {
    return _data.id;
}





