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
  
  arcbounds_sumTipFNField.h:  
  Particle boundary conditions,
  Fowler-Nordheim getting the field at r=0
  by summing over all electrons directly,
  not by using the grid

***********************************************************************/

#include "arcbound_sumTipFNField.h"
#include "mydef.h"

#include "random.h"

#define   XTRN extern
#include "var.h"
#include "arrays1.h"
#undef XTRN

#include "ArcPicConfig.h"

#include <vector>

#include <cmath>

using namespace std;

SumTipFNField::SumTipFNField(vector<char*>& options) : ofile_tip(NULL) {
  if (options.size() != 4) {
    cout << "Error in SumTipFNField(): Expected 4 options, got " << options.size() << endl;
    exit(1);
  }

  sscanf(options[0], "%*[^:]%*[:] %lg", &(this->beta_tip));
  sscanf(options[1], "%*[^:]%*[:] %lg", &(this->Remission));
  sscanf(options[2], "%*[^:]%*[:] %lg", &(this->Remission_theor));
  sscanf(options[3], "%*[^:]%*[:] %u", &(this->file_timestep));
}
SumTipFNField::~SumTipFNField() {
  if (ofile_tip == NULL) fclose(ofile_tip);
}
void SumTipFNField::print_par() const {
  printf( " - beta_tip             %g \n", this->beta_tip );
  printf( " - Remission            %g \n", this->Remission );
  printf( " - Remission_theor      %g \n", this->Remission_theor );
  printf( " - file_timestep        %u \n", this->file_timestep );
}

void SumTipFNField::inject_e(ParticleSpecies* pa, double const Ez[]) {
  //Calculate field
  double fieldConst = N_sp*1.602176565e-19/(4*PI*8.854187e-12);               // -N_sp*e/(4*pi*eps0) [V*m]
  double unitConst0 = picConfig.dz*Ldb*1e-2;                                  //dimless length (dz) -> m
  double unitConst1 = picConfig.T_ref * SQU(picConfig.dz/picConfig.Omega_pe); //dimless potential -> V
  double unitConst2 = unitConst1 * 2.0 / (Ldb*1e-2) / picConfig.dz;           //dimless field -> V/m

  double field = 0.0;
  for(size_t i = 0; i < pa->GetN(); i++) {
    double R3   = pow( SQU(pa->z[i]) + SQU(pa->r[i]), 3.0/2.0);
    double fLoc = 2*fieldConst * ( pa->z[i] / (R3 * SQU(unitConst0)) );
    //printf("%g %g %g\n", pa[i].p.z, pa[i].p.r, fLoc);
    field += fLoc;
  }
  double field2 = -(circuit->getUNz()-circuit->getU0())*unitConst1/(nz*unitConst0) + field;

  //For output
  field_PIC = Ez[0]*unitConst2;
  field_direct = field2;
  field_direct_particlesOnly = field;
  
  // printf("FIELD = %g %g %g %g %g \n\n", 
  // 	 field, 
  // 	 field2, 
  // 	 Ez[0]*unitConst2, 
  // 	 (circuit->getUNz()-circuit->getU0())*unitConst1,
  // 	 (circuit->getUNz()-circuit->getU0())*unitConst1/(nz*unitConst0) );

  //Calculate number of particles to inject (Fowler-Nordheim)
  double Eloc = -field2*beta_tip*1e-9; //GV/m
  //printf("Eloc = %g\n", Eloc);
  if (Eloc > 12.0) Eloc = 12.0;
  else if (Eloc < 0.0) return;
  
  double j = 4.7133e9 * SQU(Eloc) * exp(-62.338/Eloc); // in A/cm^2
  j *= PI*(picConfig.Ndb*e2inj_step*picConfig.Omega_pe*SQU(Remission_theor*picConfig.dz))/(6.7193e-12*picConfig.n_ref*sqrt(picConfig.T_ref)); //dimless, 2D: rem^2

  // FN INJECTION FROM FIELD EMITTER
  // FE: inject with flat distribution over the emission radius
  double v_inj_e = picConfig.v_te*0.01;
  double r1, r2; // temp variables from Gaussian RNG

  //Inject 'tmp' num particles
  size_t tmp = size_t(j);
  j -= tmp;
  if ( RAND <= j ) tmp++;
  //printf("tmp = %zu\n", tmp);
  tip_emitted = tmp; //for output
  pa->ExpandBy(tmp);
  
  for  (size_t k=0; k<tmp; k++ ) {
    // Gaussian scheme
    do { r1 = sqrt(-2.*log(RAND+1.e-20)); } while( r1 > 5. );
    r2 = RAND * TWOPI;
    pa->vr.push_back( r1*cos(r2)*v_inj_e );
    pa->vt.push_back( r1*sin(r2)*v_inj_e );
    
    // Gaussian scheme
    do { r1 = sqrt(-2.*log(RAND+1.e-20)); } while( r1 > 5. );
    pa->vz.push_back( r1*v_inj_e );
    
    pa->z.push_back( zmin + pa->vz.back());
    pa->r.push_back( Remission * sqrt(RAND) + pa->vr.back() );
    
    if ( pa->r.back() < 0 )     pa->r.back() = 1.e-20;
    else if (pa->r.back() > nr) pa->r.back() = 2*nr - pa->r.back();

    pa->m.push_back( 1 );
    
    current_e[0] += 2*( int(pa->r.back()) ) + 1;
    current_cathode[ int(pa->r.back()) ] += 1;
  }
  
  injected_e[0] += tmp;
}

void SumTipFNField::writeFile_extras(unsigned int nstep) {
  if (ofile_tip == NULL) {
    ofile_tip = fopen("arcbounds_sumTipFNField.dat", "w");
    fprintf(ofile_tip, "## timestep  field_PIC field_direct field_direct_particlesOnly tip_emitted\n");
    fflush(ofile_tip);
  }
  if (nstep % file_timestep != 0) return;

  fprintf(ofile_tip, "%11d %20f %20f %20f %10zu\n", nstep, field_PIC, field_direct, field_direct_particlesOnly, tip_emitted);
  fflush(ofile_tip);
}
