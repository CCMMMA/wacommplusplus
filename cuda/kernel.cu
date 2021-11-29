//
// Created by Ciro Giuseppe De Vita and Gennaro Mellone on 24/12/20.
//

#include <cuda.h>
#include <cuda_runtime_api.h>
#include <curand_kernel.h>

#include "../Config.hpp"
#include "../Particle.hpp"
#include "kernel.h"

// Returns -1 if a < 0 and 1 if a > 0
__host__ __device__ double sgn(double a) { return (a > 0) - (a < 0); }

// Computes the remainder of the division of a by p.
// https://gcc.gnu.org/onlinedocs/gfortran/MOD.html
__host__ __device__ double mod(double a, double p) { return a-p*(int)(a/p); }

// Returns the value of a with the sign of b.
// https://gcc.gnu.org/onlinedocs/gfortran/SIGN.html
__host__ __device__ double sign(double a, double b) { return abs(a)*sgn(b); }

__global__ void move(config_data *pConfigData, particle_data *pParticleData, int ocean_time_idx, int ocean_time, int s_w, int s_rho, int eta_rho, int xi_rho, double *pOceanTime, double *pMask, double *pLonRad, double *pLatRad, double *pDepthIntervals, double *pH, float *pZeta, float *pU, float *pV, float *pW, float *pAkt, int sizeSectionParticles, int numThread){
	int idx = threadIdx.x + blockIdx.x * blockDim.x;

	if (idx < sizeSectionParticles){
		double health0=1;

    		// Get the random flag
    		bool random=pConfigData->random;
		curandState state;
		if (random){
                	curand_init (clock64(), idx, 0, &state);
		}

    		// Get the integration time (default 30s)
    		double dti=pConfigData->dti;

    		// Get the time in seconds between two input ocean data (default 3600s, 1h)
   		double deltat=pConfigData->deltat;

    		// Ask Angelo Riccio (default 86400)
    		double tau0=pConfigData->tau0;

    		// Probability to survive (default 1.0e-4)
    		double survprob=pConfigData->survprob;

    		// Reduction Coefficient (default 1)
    		double crid=pConfigData->crid;

    		// Sedimentation velocity (m-1. default )
    		double sv=pConfigData->sv;

    		double idet=0,jdet=0,kdet=0;

    		// Number of integration intervals
    		double iint=deltat/dti;

		// For each integration interval
        	for (int t = 0; t < iint; t++) {
			// Check if the paticle is not yet active
            		if (pParticleData[idx].time > (pOceanTime[ocean_time_idx] + (t * dti))) {
                		// The particle is not already active (already emitted, but not active)
                		break;
            		}

            		// Check of the particle health is less than its probability to survive
            		if (pParticleData[idx].health < survprob) {
                		// The particle is dead
                		pParticleData[idx].health = -1;
                		// No reason to continue, exit the integration loop
                		break;
            		}

            		// Get the integer part and the fraction part of particle k
            		auto kI = (int) pParticleData[idx].k;
            		double kF = pParticleData[idx].k - kI;

            		// Get the integer part and the fraction part of particle j
            		auto jI = (int) pParticleData[idx].j;
            		double jF = pParticleData[idx].j - jI;

            		// Get the integer part and the fraction part of particle i
            		auto iI = (int) pParticleData[idx].i;
            		double iF = pParticleData[idx].i - iI;

            		// Check if the particle is out of the domain
            		if (jI < 0 || iI < 0 || jI >= eta_rho || iI >= xi_rho) {
                		// Set the particle health
                		pParticleData[idx].health = -1;

                		// no reason to continue,  exit the integration loop
                		break;
            		}

            		// Check if the particle beached
            		if (pMask[jI * xi_rho + iI] <= 0) {
                		// Set the particle health
                		pParticleData[idx].health = -1;

                		// no reason to continue,  exit the integration loop
                		break;
            		}

            		// The particle is alive!

            		// Perform the bilinear interpolation (2D) in order to get
            		// the zeta at the particle position.
            		float z1 = pZeta[ocean_time_idx * (eta_rho * xi_rho) + (jI) * xi_rho + (iI)] * (1.0 - iF) * (1.0 - jF);
	    		float z2 = pZeta[ocean_time_idx * (eta_rho * xi_rho) + (jI + 1) * xi_rho + (iI)] * (1.0 - iF) * jF;
	    		float z3 = pZeta[ocean_time_idx * (eta_rho * xi_rho) + (jI + 1) * xi_rho + (iI + 1)] * iF * jF;
	    		float z4 = pZeta[ocean_time_idx * (eta_rho * xi_rho) + (jI) * xi_rho + (iI + 1)] * iF * (1.0 - jF);

            		// The current zeta at the particle position
            		float zeta = z1 + z2 + z3 + z4;

            		// Perform the bilinear interpolation (2D) in order to get
            		// the h (depth) at the particle position.
	    		double h1 = pH[jI * xi_rho + iI] * (1.0 - iF) * (1.0 - jF);
            		double h2 = pH[(jI + 1) * xi_rho + iI] * (1.0 - iF) * jF;
            		double h3 = pH[(jI + 1) * xi_rho + (iI + 1)] * iF * jF;
            		double h4 = pH[jI * xi_rho + (iI + 1)] * iF * (1.0 - jF);

            		// The current h (depth) at the particle position
            		double h = h1 + h2 + h3 + h4;

            		// Perform the bilinear interpolation (2D) in order to get
            		// the u component of the current field in the particle position.

	    		int oy_u = (-(int)s_rho+1)*(eta_rho * xi_rho);
	    		//int oy_u = 0;
	    		float u1 = pU[(ocean_time_idx * (s_rho * eta_rho * xi_rho) + kI * (eta_rho *  xi_rho) + jI * xi_rho + iI) - oy_u] * (1.0 - iF) * (1.0 - jF);
            		float u2 = pU[(ocean_time_idx * (s_rho * eta_rho * xi_rho) + kI * (eta_rho *  xi_rho) + (jI + 1) * xi_rho + iI) - oy_u] * (1.0 - iF) * jF;
            		float u3 = pU[(ocean_time_idx * (s_rho * eta_rho * xi_rho) + kI * (eta_rho *  xi_rho) + (jI + 1) * xi_rho + (iI + 1)) - oy_u] * iF * jF;
            		float u4 = pU[(ocean_time_idx * (s_rho * eta_rho * xi_rho) + kI * (eta_rho *  xi_rho) + jI * xi_rho + (iI + 1)) - oy_u] * iF * (1.0 - jF);
	
            		// The current u component in the particle position
            		float uu = u1 + u2 + u3 + u4;

            		// Perform the bilinear interpolation (2D) in order to get
            		// the v component of the current field in the particle position.
	    		int oy_v = (-(int)s_rho+1)*(eta_rho * xi_rho);
	    		float v1 = pV[(ocean_time_idx * (s_rho * eta_rho * xi_rho) + kI * (eta_rho *  xi_rho) + jI * xi_rho + iI) - oy_v] * (1.0 - iF) * (1.0 - jF);
            		float v2 = pV[(ocean_time_idx * (s_rho * eta_rho * xi_rho) + kI * (eta_rho *  xi_rho) + (jI + 1) * xi_rho + iI) - oy_v] * (1.0 - iF) * jF;
            		float v3 = pV[(ocean_time_idx * (s_rho * eta_rho * xi_rho) + kI * (eta_rho *  xi_rho) + (jI + 1) * xi_rho + (iI + 1)) - oy_v] * iF * jF;
            		float v4 = pV[(ocean_time_idx * (s_rho * eta_rho * xi_rho) + kI * (eta_rho *  xi_rho) + jI * xi_rho + (iI + 1)) - oy_v] * iF * (1.0 - jF);

            		// The current v component in the particle position
            		float vv = v1 + v2 + v3 + v4;

            		// Perform the bilinear interpolation (3D) in order to get
            		// the w component of the current field in the particle position.
	    		int oy_w = (-(int)s_w+1)*(eta_rho * xi_rho);
	    		float w1 = pW[(ocean_time_idx * (s_w * eta_rho * xi_rho) + kI * (eta_rho *  xi_rho) + jI * xi_rho + iI) - oy_w] * (1.0 - iF) * (1.0 - jF) * (1.0 - kF);
            		float w2 = pW[(ocean_time_idx * (s_w * eta_rho * xi_rho) + kI * (eta_rho *  xi_rho) + (jI + 1) * xi_rho + iI) - oy_w] * (1.0 - iF) * jF * (1.0 - kF);
            		float w3 = pW[(ocean_time_idx * (s_w * eta_rho * xi_rho) + kI * (eta_rho *  xi_rho) + (jI + 1) * xi_rho + (iI + 1)) - oy_w] * iF * jF * (1.0 - kF);
            		float w4 = pW[(ocean_time_idx * (s_w * eta_rho * xi_rho) + kI * (eta_rho *  xi_rho) + jI * xi_rho + (iI + 1)) - oy_w] * iF * (1.0 - jF) * (1.0 - kF);
	    		float w5 = pW[(ocean_time_idx * (s_w * eta_rho * xi_rho) + (kI - 1) * (eta_rho *  xi_rho) + jI * xi_rho + iI) - oy_w] * (1.0 - iF) * (1.0 - jF) * kF;
            		float w6 = pW[(ocean_time_idx * (s_w * eta_rho * xi_rho) + (kI - 1) * (eta_rho *  xi_rho) + (jI + 1) * xi_rho + iI) - oy_w] * (1.0 - iF) * jF * kF;
            		float w7 = pW[(ocean_time_idx * (s_w * eta_rho * xi_rho) + (kI - 1) * (eta_rho *  xi_rho) + (jI + 1) * xi_rho + (iI + 1)) - oy_w] * iF * jF * kF;
            		float w8 = pW[(ocean_time_idx * (s_w * eta_rho * xi_rho) + (kI - 1) * (eta_rho *  xi_rho) + jI * xi_rho + (iI + 1)) - oy_w] * iF * (1.0 - jF) * kF;

            		// The current w component in the particle position
            		float ww = w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8;

            		// Perform the bilinear interpolation (3D) in order to get
            		// the akt in the particle position.
	    		int oy_a = (-(int)s_w+1)*(eta_rho * xi_rho);
	    		float a1 = pAkt[(ocean_time_idx * (s_rho * eta_rho * xi_rho) + kI * (eta_rho *  xi_rho) + jI * xi_rho + iI) - oy_a] * (1.0 - iF) * (1.0 - jF) * (1.0 - kF);
            		float a2 = pAkt[(ocean_time_idx * (s_rho * eta_rho * xi_rho) + kI * (eta_rho *  xi_rho) + (jI + 1) * xi_rho + iI) - oy_a] * (1.0 - iF) * jF * (1.0 - kF);
            		float a3 = pAkt[(ocean_time_idx * (s_rho * eta_rho * xi_rho) + kI * (eta_rho *  xi_rho) + (jI + 1) * xi_rho + (iI + 1)) - oy_a] * iF * jF * (1.0 - kF);
            		float a4 = pAkt[(ocean_time_idx * (s_rho * eta_rho * xi_rho) + kI * (eta_rho *  xi_rho) + jI * xi_rho + (iI + 1)) - oy_a] * iF * (1.0 - jF) * (1.0 - kF);
            		float a5 = pAkt[(ocean_time_idx * (s_rho * eta_rho * xi_rho) + (kI - 1) * (eta_rho *  xi_rho) + jI * xi_rho + iI) - oy_a] * (1.0 - iF) * (1.0 - jF) * kF;
            		float a6 = pAkt[(ocean_time_idx * (s_rho * eta_rho * xi_rho) + (kI - 1) * (eta_rho *  xi_rho) + (jI + 1) * xi_rho + iI) - oy_a] * (1.0 - iF) * jF * kF;
            		float a7 = pAkt[(ocean_time_idx * (s_rho * eta_rho * xi_rho) + (kI - 1) * (eta_rho *  xi_rho) + (jI + 1) * xi_rho + (iI + 1)) - oy_a] * iF * jF * kF;
            		float a8 = pAkt[(ocean_time_idx * (s_rho * eta_rho * xi_rho) + (kI - 1) * (eta_rho *  xi_rho) + jI * xi_rho + (iI + 1)) - oy_a] * iF * (1.0 - jF) * kF;

            		// The AKT at the particle position.
            		float aa = a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8;

            		// Evaluate the particle leap due to the current field (deterministic leap).
            		double dileap = uu * dti;
            		double djleap = vv * dti;
            		double dkleap = (sv + ww) * dti;

            		// Calculation of sigma profile
            		// sigmaPROF=sigma(Ixx,Iyy)*(1-Zdet/(-H(Ixx,Iyy))) ! Here is H
            		//double sigmaProf=oceanModelAdapter->sigma()(jI, iI)*(1-kdet/(-oceanModelAdapter->H()(jI,iI)))
            		double sigmaprof = 3.46 * (1 + pParticleData[idx].k / s_w);

            		// Extract 3 pseudorandom numbers
            		double gi = 0, gj = 0, gk = 0;
			if (random) {
            			for (int a = 0; a < 12; a++) {
                			gi = gi + curand_uniform( &state ) - 0.5;
                			gj = gj + curand_uniform( &state ) - 0.5;
                			gk = gk + curand_uniform( &state ) - 0.5;
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

            		double d1,d2,dd,jidist, kdist;

            		// Calculate the distance in radiants of latitude between the grid cell where is
            		// currently located the particle and the next one.
            		d1=pLatRad[xi_rho*(jI+1)+iI]-pLatRad[xi_rho*jI+iI];

            		// Calculate the distance in radiants of longitude between the grid cell where is
            		// currently located the particle and the next one.
            		d2=pLonRad[xi_rho*jI+iI+1]-pLonRad[xi_rho*jI+iI];

            		// Calculate the grid cell diagonal horizontal size using the Haversine method
            		// https://www.movable-type.co.uk/scripts/latlong.html
            		dd=pow(sin(0.5*d1),2) +
               		   pow(sin(0.5*d2),2)*
               		   cos(pLatRad[xi_rho*(jI+1)+iI])*
               		   cos(pLatRad[xi_rho*(jI  )+iI]);
            		jidist=2.0*atan2(pow(dd,.5),pow(1.0-dd,.5))*6371.0;

            		kdist=pDepthIntervals[kI-(-(int)s_w+2)]*(h+zeta);
            		if ( abs(kleap) > abs(kdist) ) {
                		kleap=sign(kdist,kleap);
            		}

            		// Calculate the new particle j candidate
            		jdet=pParticleData[idx].j+0.001*jleap/jidist;

            		// Calculate the new particle i candidate
            		idet=pParticleData[idx].i+0.001*ileap/jidist;

            		// Calculate the new particle k candidate
            		kdet=pParticleData[idx].k+kleap/kdist;

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

            		// Check if the candidate position is within the domain
            		if (jdetI>= 0 && idetI >= 0 && jdetI<eta_rho && idetI <xi_rho) {
                		// Check if the candidate new particle position is on land (mask=0)
                		if (pMask[xi_rho*jdetI+idetI] <= 0.0) {
                    			// Reflect the particle
                    			if (idetI < iI) {
                        			idet = (double) iI + abs(pParticleData[idx].i - idet);
                    			} else if (idetI > iI) {
                        			idet = (double) idetI - mod(idet, 1.0);
                    			}
                    			if (jdetI < jdet) {
                       	 			jdet = (double) jdetI + abs(pParticleData[idx].j - jdet);
                    			} else if (jdetI > jI) {
                        			jdet = (double) jdetI - mod(jdet, 1.0);
                    			}
                		}
            		}


            		// Assign the new particle position
            		pParticleData[idx].i=idet;
            		pParticleData[idx].j=jdet;
            		pParticleData[idx].k=kdet;

            		// Update the paticle age
            		pParticleData[idx].age=pParticleData[idx].age+dti;
            		// Decay the particle
            		pParticleData[idx].health=health0*exp(-pParticleData[idx].age/tau0);
        	}
	}
}

cudaError_t cudaMoveParticle(config_data *pConfigData, particle_data *pParticleData, int ocean_time_idx, int ocean_time, int s_w, int s_rho, int eta_rho, int xi_rho, double *pOceanTime, double *pMask, double *pLonRad, double *pLatRad, double *pDepthIntervals, double *pH, float *pZeta, float *pU, float *pV, float *pW, float *pAkt, int sizeSectionParticles, int numThread, int numGPU){

    dim3 nBlocks, nThreadPerBlock = 512;

    nBlocks = sizeSectionParticles/nThreadPerBlock.x + ((sizeSectionParticles%nThreadPerBlock.x) == 0?0:1);

    move<<<nBlocks, nThreadPerBlock>>>(pConfigData, pParticleData, ocean_time_idx, ocean_time, s_w, s_rho, eta_rho, xi_rho, pOceanTime, pMask, pLonRad, pLatRad, pDepthIntervals, pH, pZeta, pU, pV, pW, pAkt, sizeSectionParticles, numThread);

    cudaDeviceSynchronize();

	return cudaGetLastError();
}
