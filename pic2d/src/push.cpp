/**********************************************************************

  2D3V PIC-MCC CODE '2D Arc-PIC'

  Copyright 2010-2015 CERN and Helsinki Institute of Physics.
  This software is distributed under the terms of the
  GNU General Public License version 3 (GPL Version 3),
  copied verbatim in the file LICENCE.md. In applying this
  license, CERN does not waive the privileges and immunities granted to it
  by virtue of its status as an Intergovernmental Organization
  or submit itself to any jurisdiction.

  Project website: http://arcpic.web.cern.ch/
  Developers: Helga Timko, Kyrre Sjobak

  push.cpp:
  Particle mover; eletrostatic treatment only (optional B_ext)

***********************************************************************/

#include <stdio.h>
#include <math.h>

#include  "pic.h"
#include  "dim.h"
#include  "mydef.h"

#include "push.h"

void  push_2D( ParticleSpecies* pa, const double Eg_r[], const double Eg_z[], int NZ ) {
  
  for(size_t n=0; n < pa->GetN(); n++) {
    double hr  = pa->r[n];
    int j      = (int)hr;
    hr -= j;
    
    double hz  = pa->z[n];
    int k   = (int)hz;
    hz -= k;
    
    // Field Interpolation
    double Ez  = (1-hr)*(1-hz)*Eg_z[j*NZ+k] + (1-hr)*hz*Eg_z[j*NZ+(k+1)] + hr*(1-hz)*Eg_z[(j+1)*NZ+k] + hr*hz*Eg_z[(j+1)*NZ+(k+1)];
    double Er  = (1-hr)*(1-hz)*Eg_r[j*NZ+k] + (1-hr)*hz*Eg_r[j*NZ+(k+1)] + hr*(1-hz)*Eg_r[(j+1)*NZ+k] + hr*hz*Eg_r[(j+1)*NZ+(k+1)];
    
    // Acceleration
    pa->vz[n] -= 2.*Ez;
    pa->vr[n] -= 2.*Er;
    
    // Move particle
    double r0 = pa->r[n];
    pa->z[n] += pa->vz[n];
    pa->r[n]  = sqrt( SQU(r0 + pa->vr[n]) + SQU(pa->vt[n]) );
    
    // Rotate coordinate system
    double sa, ca;
    if ( pa->r[n] > 1.e-8 ) {
      sa = pa->vt[n] / pa->r[n];
      ca = (r0 + pa->vr[n]) / pa->r[n];
    }
    else {
      if ( pa->vr[n] < -1.e-8 ) {
	sa = 0.;
	ca = -1.;
      }
      else {
	sa = 0.;
	ca = 1.;
      }
    }
    
    double vr = pa->vr[n];
    pa->vr[n] =  ca * vr + sa * pa->vt[n];
    pa->vt[n] = -sa * vr + ca * pa->vt[n];
  }
}

void  push_magnetic_onlyBz_2D( ParticleSpecies* pa, const double Eg_r[], const double Eg_z[], 
			       const double Bz, double factor, int NZ ) {
  
  double B = factor*Bz;
  double B_alt = 2*B/(1+B*B);
  
  for(size_t n=0; n < pa->GetN(); n++) {
    double hr  = pa->r[n];
    int    j   = (int)hr;
    hr -= j;
    
    double hz  = pa->z[n];
    int k   = (int)hz;
    hz -= k;
    
    // Field Interpolation
    double Ez  = (1-hr)*(1-hz)*Eg_z[j*NZ+k] + (1-hr)*hz*Eg_z[j*NZ+(k+1)] + hr*(1-hz)*Eg_z[(j+1)*NZ+k] + hr*hz*Eg_z[(j+1)*NZ+(k+1)];
    double Er  = (1-hr)*(1-hz)*Eg_r[j*NZ+k] + (1-hr)*hz*Eg_r[j*NZ+(k+1)] + hr*(1-hz)*Eg_r[(j+1)*NZ+k] + hr*hz*Eg_r[(j+1)*NZ+(k+1)];
    
    // Acceleration: Boris
    double vza = pa->vz[n] - Ez;
    double vra = pa->vr[n] - Er;
    
    double vrb = vra - pa->vt[n]*B;
    double vtb = pa->vt[n] + vra*B;
    
    double vrc = vra - vtb*B_alt;
    
    pa->vt[n] += vrb*B_alt;
    pa->vz[n]  = vza - Ez;
    pa->vr[n]  = vrc - Er;
    
    // Move particle
    double r0 = pa->r[n];
    pa->z[n] += pa->vz[n];
    pa->r[n]  = sqrt( SQU(r0 + pa->vr[n]) + SQU(pa->vt[n]) );
    
    // Rotate coordinate system
    double sa, ca;
    if ( pa->r[n] > 1.e-8 ) {
      sa = pa->vt[n]/pa->r[n];
      ca = (r0 + pa->vr[n])/pa->r[n];
    }
    else {
      if ( pa->vr[n] < -1.e-8 ) {
	sa = 0.;
	ca = -1.;
      }
      else {
	sa = 0.;
	ca = 1.;
      }
    }
    double vr = pa->vr[n];
    pa->vr[n]  = ca * vr + sa * pa->vt[n];
    pa->vt[n] = -sa * vr + ca * pa->vt[n];
  }

}

void  push_magnetic_2D( ParticleSpecies* pa, const double Eg_r[], const double Eg_z[],
			const double Bextz, const double Bextt, double factor, int NZ ) {
  
  double Bz = factor*Bextz;
  double Bt = factor*Bextt;
  double B_alt = 2./(1.+Bz*Bz+Bt*Bt);
  
  for(size_t n=0; n < pa->GetN(); n++) {
    double hr  = pa->r[n];
    int j   = (int)hr;
    hr -= j;
    
    double hz  = pa->z[n];
    int k   = (int)hz;
    hz -= k;
    
    // Field Interpolation
    double Ez  = (1-hr)*(1-hz)*Eg_z[j*NZ+k] + (1-hr)*hz*Eg_z[j*NZ+(k+1)] + hr*(1-hz)*Eg_z[(j+1)*NZ+k] + hr*hz*Eg_z[(j+1)*NZ+(k+1)];
    double Er  = (1-hr)*(1-hz)*Eg_r[j*NZ+k] + (1-hr)*hz*Eg_r[j*NZ+(k+1)] + hr*(1-hz)*Eg_r[(j+1)*NZ+k] + hr*hz*Eg_r[(j+1)*NZ+(k+1)];
    
    // Acceleration: Boris
    double vza = pa->vz[n] - Ez;
    double vra = pa->vr[n] - Er;
    
    double vzb = vza - vra*Bt;
    double vrb = vra - pa->vt[n]*Bz + vza*Bt;
    double vtb = pa->vt[n] + vra*Bz;
    
    double vzc = vza - vrb*Bt*B_alt;
    double vrc = vra - (vtb*Bz - vzb*Bt)*B_alt;
    pa->vt[n] += vrb*Bz*B_alt;
    
    pa->vz[n]  = vzc - Ez;
    pa->vr[n]  = vrc - Er;
    
    // Move particle
    double r0 = pa->r[n];
    pa->z[n] += pa->vz[n];
    pa->r[n] = sqrt( SQU(r0 + pa->vr[n]) + SQU(pa->vt[n]) );
    
    // Rotate coordinate system
    double sa, ca;
    if ( pa->r[n] > 1.e-8 ) {
      sa = pa->vt[n] / pa->r[n];
      ca = (r0 + pa->vr[n])/pa->r[n];
    }
    else {
      if ( pa->vr[n] < -1.e-8 ) {
	sa = 0.;
	ca = -1.;
      }
      else {
	sa = 0.;
	ca = 1.;
      }
    }
    double vr = pa->vr[n];
    pa->vr[n] = ca * vr + sa * pa->vt[n];
    pa->vt[n] = -sa * vr + ca * pa->vt[n];
  }
}

void  push_neutrals_2D( ParticleSpecies* pa ) {
  
  for(size_t n=0; n < pa->GetN(); n++) {
    // No acceleration
    
    // Move particle
    double r0 = pa->r[n];
    pa->z[n] += pa->vz[n];
    pa->r[n] = sqrt( SQU(r0 + pa->vr[n]) + SQU(pa->vt[n]) );
    
    // Rotate coordinate system
    double sa, ca;
    if ( pa->r[n] > 1.e-8 ) {
      sa = pa->vt[n] / pa->r[n];
      ca = (r0 + pa->vr[n])/pa->r[n];
    }
    else {
      if ( pa->vr[n] < -1.e-8 ) {
	sa = 0.;
	ca = -1.;
      }
      else {
	sa = 0.;
	ca = 1.;
      }
    }
    double vr = pa->vr[n];
    pa->vr[n] =  ca * vr + sa * pa->vt[n];
    pa->vt[n] = -sa * vr + ca * pa->vt[n];
  }

}
