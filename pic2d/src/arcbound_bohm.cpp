/**********************************************************************

  2D3V PIC-MCC CODE '2D Arc-PIC'

  Copyright 2010-2019 CERN and Helsinki Institute of Physics.
  This software is distributed under the terms of the
  GNU General Public License version 3 (GPL Version 3),
  copied verbatim in the file LICENCE.md. In applying this
  license, CERN does not waive the privileges and immunities granted to it
  by virtue of its status as an Intergovernmental Organization
  or submit itself to any jurisdiction.

  Project website: http://arcpic.web.cern.ch/
  Developers: Kyrre Sjobak, Riccardo Cassegrande

  arcbound_bohm.cpp:
  Particle boundary conditions, implements a Bohm sheet condition

***********************************************************************/

#include "arcbound_bohm.h"
#include "mydef.h"

#include "random.h"

#define   XTRN extern
#include "var.h"
#include "arrays1.h"
#undef XTRN

#include "ArcPicConfig.h"

#include <cmath>

using namespace std;

ArcBohm::ArcBohm(vector<char*>& options) {

  if (options.size() != 17) {
    cout << "Error in ArcOriginalNewHS(): Expected 16 options, got "
	 << options.size() << endl;
    exit(1);
  }
/*
  //Interpret lines
  char foo;

  sscanf(options[0],  "%*[^:]%*[:] %lg", &(this->beta_tip));
  sscanf(options[1],  "%*[^:]%*[:] %lg", &(this->beta_f));
  sscanf(options[2],  "%*[^:]%*[:] %lg", &(this->j_melt));

  sscanf(options[3],  "%*[^:]%*[:] %c",  &(foo));
  if      (foo == 'y') this->do_erosion = true;
  else if (foo == 'n') this->do_erosion = false;
  else {
    printf("Error in ArcOriginalNewHS::ArcOriginalNewHS(): do_erosion has to be either 'y' or 'n' \n");
    exit(1);
  }  
  if (this->do_erosion and (beta_tip < beta_f) ) {
    printf("Error in ArcOriginalNewHS::ArcOriginalNewHS(): when do_erosion = true, must have beta_tip > beta_flat\n"); 
    printf(" (else beta_tip gets set to beta_flat when doing the erosion)\n");
    exit(1);
  }

  sscanf(options[4],  "%*[^:]%*[:] %lg", &(this->Remission));
  sscanf(options[5],  "%*[^:]%*[:] %lg", &(this->Remission_theor));
  sscanf(options[6],  "%*[^:]%*[:] %lg", &(this->Rborder));

  sscanf(options[7],  "%*[^:]%*[:] %lg", &(this->SEY));
  sscanf(options[8],  "%*[^:]%*[:] %lg", &(this->r_Cu_e));
  sscanf(options[9],  "%*[^:]%*[:] %lg", &(this->r_Cu_e_flat));


  sscanf(options[10],  "%*[^:]%*[:] %c",  &(foo));
  if      (foo == 'y') this->evap_flat_center = true;
  else if (foo == 'n') this->evap_flat_center = false;
  else {
    printf("Error in ArcOriginalNewHS::ArcOriginalNewHS(): evap_flat_center has to be either 'y' or 'n' \n");
    exit(1);
  }
  
  sscanf(options[11],  "%*[^:]%*[:] %c",  &(foo));
  if      (foo == 'y') this->fracInjectStep = true;
  else if (foo == 'n') this->fracInjectStep = false;
  else {
    printf("Error in ArcOriginalNewHS::ArcOriginalNewHS(): fracInjectStep has to be either 'y' or 'n' \n");
    exit(1);
  }

  sscanf(options[12],  "%*[^:]%*[:] %lg", &(this->alpha_flat));

  sscanf(options[13],  "%*[^:]%*[:] %lg", &(this->heatspike_threshold));
  sscanf(options[14],  "%*[^:]%*[:] %lg", &(this->heatspike_yield_ratio));

  sscanf(options[15],  "%*[^:]%*[:] %c",  &(foo));
  if      (foo == 'a' or foo == 'b' or foo == 'c') {
    this->heatspike_model = foo;
  }
  else {
    printf("Error in ArcOriginalNewHS::ArcOriginalNewHS(): heatspike_model has to be either 'a', 'b' or 'c' \n");
    exit(1);
  }

  sscanf(options[16], "%*[^:]%*[:] %u",  &(this->file_timestep));
*/
}

ArcBohm::~ArcBohm() {
  
}

void ArcBohm::print_par() const {
    /*
  printf( " - beta_tip:            %g \n",          this->beta_tip         );
  printf( " - beta_f:              %g \n",          this->beta_f           );
  printf( " - j_melt:              %g [A/cm^2]\n",  this->j_melt           );
  double j_max = 4.7133e9 * SQU(12.0) * exp(-62.338/12.0);
  if (j_melt > j_max) {
    printf( "   (NOTE: j_melt > j_max = %g [A/cm^2] -- No melting possible!)\n", j_max);
  }
  printf( " - do_erosion:          %c \n", this->do_erosion ? 'y' : 'n'    );
  printf( " - Remission:           %g [dz]\n",      this->Remission        );
  printf( " - Remission_theor:     %g [dz]\n",      this->Remission_theor  );
  printf( " - evap_flat_center:    %c \n", this->evap_flat_center ? 'y' : 'n');
  printf( " - Rborder:             %g [dz]\n",      this->Rborder          );
  printf( " - SEY:                 %g \n",          this->SEY              );
  printf( " - r_Cu_e:              %g \n",          this->r_Cu_e           );
  printf( " - r_Cu_e_flat:         %g \n",          this->r_Cu_e_flat      );
  printf( " - evap_flat_center:    %c \n", this->evap_flat_center ? 'y' : 'n');
  printf( " - alpha_flat:          %g \n",          this->alpha_flat       );
  printf( " - fracInjectStep:      %c \n", this->fracInjectStep ? 'y' : 'n');
  printf( " - heatspike_threshold: %g [particles/cm^2/s]\n",  this->heatspike_threshold);
  printf( "   = %g [particles/Ldb^2/injection_timestep]\n",
          this->heatspike_threshold*(picConfig.Ndb*n2inj_step*picConfig.Omega_pe) /
          (4.1938226e7*picConfig.n_ref*sqrt(picConfig.T_ref)) );
  printf( "   = %g [particles/inner cell/injection_timestep]\n",
          this->heatspike_threshold*(picConfig.Ndb*n2inj_step*picConfig.Omega_pe) /
          (4.1938226e7*picConfig.n_ref*sqrt(picConfig.T_ref)) * PI*SQU(picConfig.dr) );
  printf( "   = %g [particles/2nd cell/injection_timestep]\n",
          this->heatspike_threshold*(picConfig.Ndb*n2inj_step*picConfig.Omega_pe) /
          (4.1938226e7*picConfig.n_ref*sqrt(picConfig.T_ref)) * 3*PI*SQU(picConfig.dr) );
  printf( " - heatspike_yieldRatio: %g [outgoing/incomming]\n", this->heatspike_yield_ratio);
  printf( " - heatspike_model:     %c \n",          this->heatspike_model  );
  printf( " - file_timestep:       %u \n",          this->file_timestep    );
  */
}

void ArcBohm::init(unsigned int nr, double zmin, double zmax, double rmax) {
  ArcBounds::init(nr,zmin,zmax,rmax);
  
  /*
  v_inj_e = picConfig.v_te*0.01;
  v_inj_i = vt_ions[2];

  sput_cathode_current = new unsigned int[nr];
  for (unsigned int i = 0; i < nr; i++) {
    sput_cathode_current[i] = 0;
  }

  //Initialization
  this->has_melted          = false;
  this->erosion             = 0;
  this->emitted_tip_melting = 0;
  
  this->emitted_tip_evap    = 0;
  this->emitted_flat_evap   = new int[nr];
  for (unsigned int i = 0; i < nr; i++) {
    emitted_flat_evap[i] = 0;
  }

  this->emitted_tip_output         = 0;
  this->emitted_flat_output        = 0;
  this->emitted_SEY_output         = 0;
  this->emitted_evap_output        = 0;
  this->emitted_sputter_cat_output = 0;
  this->emitted_sputter_ano_output = 0;
  this->emitted_heatspike_output   = 0;

  this->heatspike_sigma    = 0;
  this->heatspike_incident = 0;
*/
}

void ArcBohm::re_init(unsigned int nr, double zmin, double zmax, double rmax) {
  //ArcBounds::re_init(nr, zmin, zmax, rmax);
  ArcBohm::init(nr, zmin, zmax, rmax); //ArcBounds::re_init just calls init()
}



void ArcBohm::timestep(unsigned int nstep, bool isOutputTimestep) {
  //Write files and reset counting arrays
  ArcBounds::timestep(nstep, isOutputTimestep);
}

void ArcBohm::inject_i(ParticleSpecies* pa, double const Ez[], unsigned int sort) {
    //TODO
};