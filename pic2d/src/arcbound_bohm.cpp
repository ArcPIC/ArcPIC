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

  if (options.size() != 5) {
    cout << "Error in ArcBohm(): Expected 5 options, got "
	 << options.size() << endl;
    exit(1);
  }

  //Interpret lines
  char foo;

  sscanf(options[0],  "%*[^:]%*[:] %lg", &(this->v_injAnode_ion));
  sscanf(options[1],  "%*[^:]%*[:] %lg", &(this->v_injAnode_temp));
  sscanf(options[2],  "%*[^:]%*[:] %lg", &(this->rate_injAnode));
  sscanf(options[3],  "%*[^:]%*[:] %u",  &(this->maxR_injAnode));

  sscanf(options[4], "%*[^:]%*[:] %u",  &(this->file_timestep));
}

ArcBohm::~ArcBohm() {
}

void ArcBohm::print_par() const {
  printf( " - v_injAnode_ion:     %g \n",          this->v_injAnode_ion         );
  printf( " - v_injAnode_temp:    %g \n",          this->v_injAnode_temp        );
  printf( " - rate_injAnode:      %g [UNIT]\n",    this->rate_injAnode          );
  printf( " - maxR_injAnode:      %g [dz]\n",      this->maxR_injAnode          );

  printf( " - file_timestep:      %u \n",          this->file_timestep          );
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