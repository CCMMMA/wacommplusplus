//
// Created by Raffaele Montella on 12/5/20.
//

#include "Particle.hpp"
#include "Config.hpp"

Particle::Particle(double j, double k, double i,
                   double health, double tpart):
                   j(j), k(k), i(i),
                   health(health),tpart(tpart)
{
}

Particle::~Particle() {
}

bool Particle::isAlive() {
    return health>0;
}



void Particle::move(Config &config,
                    int ocean_time_idx,
                    Array::Array1<double> &depth, Array::Array2<double> &zeta,
                    Array::Array2<double> &lon, Array::Array2<double> &lat,
                    Array::Array2<double> &mask,
                    Array::Array4<double> &u, Array::Array4<double> &v, Array::Array4<double> &w,
                    Array::Array4<double> &akt) {

    double dti=config.Dti();
    double deltat=config.Deltat();
    double tau0=config.Tau0();
    double survprob=config.Survprob();

    double idet=0,jdet=0,kdet=0;
    double crid=1;
    double vs=0;

    double iint=deltat/dti;
    double t=0;

    size_t s_w = w.Nx();
    size_t eta_rho = mask.Nx();
    size_t xi_rho = mask.Ny();

    for (unsigned t=0;t<iint;t++) {
        if (health<survprob) {
            pstatus=0;
            return;
        }
        unsigned int iI=(unsigned int)i; double iF=i-iI;
        unsigned int jI=(unsigned int)j; double jF=j-jI;
        unsigned int kI=(unsigned int)k; double kF=k-kI;

        if (mask(jI,iI)<=0.0 || jI>=eta_rho|| iI>=xi_rho) {
            health=-1;
            pstatus=0;
            return;
        }

        if (k>=s_w) {
            health=-13;
            pstatus=0;
            return;
        }

        //printf("Alive:%u kI:%u jI:%u iI:%u kF:%f jF:%f iF:%f\n",id,kI,jI,iI,kF,jF,iF);

        double u1=u(ocean_time_idx,kF,jF,iF)*(1-iF)*(1-jF);
        double u2=u(ocean_time_idx,kF,jF+1,iF)*(1-iF)*jF;
        double u3=u(ocean_time_idx,kF,jF+1,iF+1)*iF*jF;
        double u4=u(ocean_time_idx,kF,jF,iF+1)*iF*(1-jF);
        double uu=u1+u2+u3+u4;

        //printf("uu:%f\n",uu);

        double v1=v(ocean_time_idx,kF,jF,iF)*(1-iF)*(1-jF);
        double v2=v(ocean_time_idx,kF,jF+1,iF)*(1-iF)*jF;
        double v3=v(ocean_time_idx,kF,jF+1,iF+1)*iF*jF;
        double v4=v(ocean_time_idx,kF,jF,iF+1)*iF*(1-jF);
        double vv=v1+v2+v3+v4;

        //printf("vv:%f\n",vv);

        double w1=w(ocean_time_idx,kF,jF,iF)*(1-iF)*(1-jF)*(1-kF);
        double w2=w(ocean_time_idx,kF,jF+1,iF)*(1-iF)*jF*(1-kF);
        double w3=w(ocean_time_idx,kF,jF+1,iF+1)*iF*jF*(1-kF);
        double w4=w(ocean_time_idx,kF,jF,iF+1)*iF*(1-jF)*(1-kF);
        double w5=w(ocean_time_idx,kF+1,jF,iF)*(1-iF)*(1-jF)*kF;
        double w6=w(ocean_time_idx,kF+1,jF+1,iF)*(1-iF)*jF*kF;
        double w7=w(ocean_time_idx,kF+1,jF+1,iF+1)*iF*jF*kF;
        double w8=w(ocean_time_idx,kF+1,jF,iF+1)*iF*(1-jF)*kF;
        double ww=w1+w2+w3+w4+w5+w6+w7+w8;

        //printf("ww:%f\n",ww);

        double ileap=uu*dti;
        double jleap=vv*dti;
        double kleap=(vs+ww)*dti;

        //printf ("kleap:%f jleap:%f ileap:%f\n",kleap,jleap,ileap); 

        double sigmaprof=3.46*(1+kdet/s_w);
        double gi=0,gj=0,gk=0;
        for (int a=0;a<12;a++) {
            gi=gi+gen()-0.5;
            gj=gj+gen()-0.5;
            gk=gk+gen()-0.5;
        }

        double a1=akt(ocean_time_idx,kF,jF,iF)*(1-iF)*(1-jF)*(1-kF);
        double a2=akt(ocean_time_idx,kF,jF+1,iF)*(1-iF)*jF*(1-kF);
        double a3=akt(ocean_time_idx,kF,jF+1,iF+1)*iF*jF*(1-kF);
        double a4=akt(ocean_time_idx,kF,jF,iF+1)*iF*(1-jF)*(1-kF);
        double a5=akt(ocean_time_idx,kF+1,jF,iF)*(1-iF)*(1-jF)*kF;
        double a6=akt(ocean_time_idx,kF+1,jF+1,iF)*(1-iF)*jF*kF;
        double a7=akt(ocean_time_idx,kF+1,jF+1,iF+1)*iF*jF*kF;
        double a8=akt(ocean_time_idx,kF+1,jF,iF+1)*iF*(1-jF)*kF;
        double aa=a1+a2+a3+a4+a5+a6+a7+a8;
        //printf("aa:%f\n",aa);

        double rileap=gi*sigmaprof;
        double rjleap=gj*sigmaprof;
        double rkleap=gk*aa*crid;

        //printf ("rkleap:%f rjleap:%f rileap:%f\n",rkleap,rjleap,rileap);

        ileap=ileap+rileap;
        jleap=jleap+rjleap;
        kleap=kleap+rkleap;

        double d1,d2,dist;

        d1=(lat(jI+1,iI)-lat(jI,iI));
        d2=(lon(jI,iI+1)-lon(jI,iI));
        d1=pow(sin(0.5*d1),2) + pow(sin(0.5*d2),2)* cos(lat(jI+1,iI))*cos(lat(jI,iI));
        dist=2.0*atan2(pow(d1,.5),pow(1.0-d1,.5))*6371.0;
        idet=i+0.001*ileap/dist;

        d1=(lat(jI+1,iI)-lat(jI,iI));
        d2=(lon(jI,iI+1)-lon(jI,iI));
        d1=pow(sin(0.5*d1),2) + pow(sin(0.5*d2),2)* cos(lat(jI+1,iI))*cos(lat(jI,iI));
        dist=2.0*atan2(sqrt(d1),pow(1.0-d1,.5))*6371.0;
        jdet=j+0.001*jleap/dist;

        //printf("depth:%f, zeta:%f\n",depth[kI],zeta[jI][iI]);
        dist=depth(kI)*zeta(jI,iI);
        if (dist<=0) dist=1e-4;
        //printf("dist:%f, kleap:%f\n",dist,kleap);
        if ( abs(kleap) > abs(dist) ) {
            kleap=sign(dist,kleap);
        }
        kdet=k+kleap/dist;

        //printf("kdet:%f jdet:%f idet:%f\n",kdet,jdet,idet);

        if ( kdet >= s_w ) {
            kdet=s_w-1;
        }
        if ( kdet < 0. ) {
            kdet=0;
        }

        //printf("kdet:%f jdet:%f idet:%f\n",kdet,jdet,idet);

        unsigned int jdetI=(unsigned int)jdet;
        unsigned int idetI=(unsigned int)idet;
        if ( mask(jdetI,idetI) <= 0.0 ) {
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

        i=idet;
        j=jdet;
        k=kdet;

        time=time+dti;
        health=health0*exp(-time/tau0);
    }
}

double Particle::K() { return k; }
double Particle::J() { return j; }
double Particle::I() { return i; }

double Particle::gen() { return 0; }


