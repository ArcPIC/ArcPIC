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

  arcbound_flexFN2.cpp:
  Particle boundary conditions,
  Fowler-Nordheim where beta and alpha are functions of the r-coordinate

***********************************************************************/

#include "arcbound_flexFN2.h"
#include "mydef.h"
#include "random.h"

#define   XTRN extern
#include "var.h"
#include "arrays1.h"
#undef XTRN

#include "ArcPicConfig.h"

#include <cmath>

using namespace std;

// ******** Implementation of FlexFN2 ******************
void FlexFN2::injectFNring(ParticleSpecies* pa, double alpha, double beta, double const Ez[], double r1, double r2) {

  double v_inj_e = picConfig.v_te*0.01;

  for (size_t i = (size_t) r1; i <= (size_t) (r2+0.6); i++) { //Loop over meshpoints
    // r2+0.6 in case nr-0.5 < r2 < nr, needs to activate outer ring.
    // Use 0.6 to avoid roundoff problems (even if 0.5 *should* be sufficient).

    //Emit between max((i-0.5), r1), min(i+0.5,r2)
    double R1 = r1 < i-0.5 ? i-0.5 : r1;
    double R2 = i+0.5 < r2 ? i+0.5 : r2;
    if (R2 <= R1) continue; //0-area rings may be generated due to r2+0.6

    //Sanity check
    if (R1 < 0 || R2 > nr) {
      printf("Error in FlexFN2::injectFNring(): i=%zu, nr=%u, r1=%g, r2=%g\n", i, nr, r1, r2);
      printf("Check your emitter boundaries!\n");
      exit(1);
    }
    
    
    // *** CURRENT CALCULATION *** //
    
    double field = Ez[i*NZ];
    if (field > 0) continue;
        
    // Fowler-Nordheim with Wang-Loew approximation
    // W-L: v(y) ~ 0.956 - 1.062*(3.7947e-5)^2*Eloc/(4.5*4.5)
    // work fct=4.5eV, j in A/cm^2, E in V/m
    
    // rescale the field to GV/m, multiply with beta(r_i)
    field = - 2.69036254e-10*picConfig.dz/SQU(picConfig.Omega_pe)*sqrt(picConfig.T_ref*picConfig.n_ref)*field*beta;
    // Protect against numerical fluctuations
    // (this caps j at the field where FN becomes invalid;
    //  Murphy&Good eq. 57, T=0K, phi=4.5eV).
    if ( field > 12. ) field = 12.;      
    double I_FN = 4.7133e9 * SQU(field) * exp(-62.338/field); // in A/cm^2
    
    //Rescale to units (#Superparticles / omega_pe^-1) / lambda_Db^2
    I_FN *= picConfig.Ndb/(6.7192539e-12*picConfig.n_ref*sqrt(picConfig.T_ref));
    //Area factor
    I_FN *= alpha;
    // [#superparticles]
    I_FN *= PI*(SQU(R2*picConfig.dz)-SQU(R1*picConfig.dz))*(e2inj_step*picConfig.Omega_pe);

    // *** INJECTION *** //
    
    //Number to inject
    size_t Ninj = (size_t) I_FN;
    if ( RAND <= I_FN-Ninj ) Ninj++;
    pa->ExpandBy(Ninj);
    
    //Inject Ninj superparticles
    double r1; //Temp variable
    for (size_t k = 0; k < Ninj; k++) {
      
      //Velocity
      do { r1 = sqrt(-2.*log(RAND+1.e-20)); } while( r1 > 5. );
      double r2 = RAND * TWOPI;
      pa->vr.push_back( r1*cos(r2)*v_inj_e );
      pa->vt.push_back( r1*sin(r2)*v_inj_e );
      do { r1 = sqrt(-2.*log(RAND+1.e-20)); } while( r1 > 5. );
      pa->vz.push_back( r1*v_inj_e );
      
      //Position
      // uniform random position on annulus (R1, R2) as found above
      pa->r.push_back( sqrt( (SQU(R2)-SQU(R1))*RAND + SQU(R1) ) );
      pa->z.push_back( zmin );
      
      //Fractional timestep push (euler method)
      double Rp = RAND; //1-R; how (fractionally) far into the timestep
      // are we at z=0?
      pa->r.back()  += Rp*e2inj_step*pa->vr.back();
      pa->z.back()  += Rp*e2inj_step*pa->vz.back();
      //Assume no magnetic field, E = Ez on surface
      pa->vz.back() -= (Rp*e2inj_step-0.5)*2*Ez[i*NZ];
      
      if ( pa->r.back() < 0 )      pa->r.back() = -pa->r.back(); //Reflect on axis
      else if (pa->r.back() > nr ) pa->r.back() = 2*nr - pa->r.back();
      
      pa->m.push_back( 1 );
      
      current_e[0] += 2*( int(pa->r.back()) ) + 1;
      current_cathode[ int(pa->r.back()) ] += 1;
    }

    injected_e[0] += Ninj;
    
  } // END loop over meshpoints
}

// ******** Implementation of FlexFN2_ring ******************
FlexFN2_ring::FlexFN2_ring(std::vector<char*>& options) {
  if (options.size() != 5) {
    cout << "Error in FlexFN2_ring(): Expected 5 options, got "
	 << options.size() << endl;
    exit(1);
  }

  sscanf(options[0], "%*[^:]%*[:] %lg", &(this->ring_alpha));
  sscanf(options[1], "%*[^:]%*[:] %lg", &(this->ring_beta));
  sscanf(options[2], "%*[^:]%*[:] %lg", &(this->ring_r1));
  sscanf(options[3], "%*[^:]%*[:] %lg", &(this->ring_r2));
  sscanf(options[4], "%*[^:]%*[:] %u",  &(this->file_timestep));
}
void FlexFN2_ring::print_par() const {
  printf( " - ring_alpha           %g \n", ring_alpha);
  printf( " - ring_beta            %g \n", ring_beta);
  printf( " - ring_r1              %g \n", ring_r1);
  printf( " - ring_r2              %g \n", ring_r2);
  printf( " - file_timestep        %u \n", file_timestep );
}

void FlexFN2_ring::init(unsigned int nr, double zmin, double zmax, double rmax) {
  ArcBounds::init(nr,zmin,zmax,rmax);
}
void FlexFN2_ring::re_init(unsigned int nr, double zmin, double zmax, double rmax) {
  FlexFN2_ring::init(nr, zmin, zmax, rmax);
}

void FlexFN2_ring::inject_e(ParticleSpecies* pa, double const Ez[]) {
  injectFNring(pa, ring_alpha, ring_beta, Ez, ring_r1, ring_r2);
}

