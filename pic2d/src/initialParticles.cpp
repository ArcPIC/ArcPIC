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

  initialParticles.cpp:
  Pre-fill the system with some particles
   (neutral background gas etc.)

***********************************************************************/

#include "initialParticles.h"
#include "input.h"
#include "random.h"

#define XTRN extern
#include "var.h"
#include "arrays1.h"
#include "dim.h"
#include "mydef.h"
#undef XTRN

#include "ArcPicConfig.h"

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cmath>

using namespace std;

InitialParticles* InitialParticles::LoadInitialParticles(FILE* in_file) {
  vector<char*> options;
  char* type = readInputSection(in_file, options, true);

  InitialParticles* ret = NULL;
    
  if (! strncmp(type, "None", NAME_MAXLEN) ) {
    cout << "Not using an initial particle distribution, iParts==NULL" << endl;
  }
  else if (! strncmp(type, "UniformRestricted", NAME_MAXLEN) ) {
    ret = new UniformRestricted(options);
  }
  else {
    cout << "Unknown InitialParticle type '" << type << "'" << endl;
    exit(1);
  }

  //Cleanup
  delete[] type;
  for (size_t i = 0; i < options.size(); i++) {
    delete[] options[i];
  }
  options.clear();
  
  return ret;
}

UniformRestricted::UniformRestricted (std::vector<char*> options) {
  if (options.size() != 9) {
    cout << "Error in UniformRestricted::UniformRestricted(): Expected 9 options, got "
	 << options.size() << endl;
    exit(1);
  }

  sscanf(options[0], "%*[^:]%*[:] %lg", &density);
  sscanf(options[1], "%*[^:]%*[:] %lg", &maxR);
  sscanf(options[2], "%*[^:]%*[:] %lg", &maxZ);
  
  sscanf(options[3], "%*[^:]%*[:] %lg", &minR);
  sscanf(options[4], "%*[^:]%*[:] %lg", &minZ);

  char foo;
  sscanf(options[5],  "%*[^:]%*[:] %c",  &(foo));
  if      (foo == 'y') doInject_e = true;
  else if (foo == 'n') doInject_e = false;
  else {
    printf("Error in UniformRestricted::UniformRestricted(): doInject_e has to be either 'y' or 'n' \n");
    exit(1);
  }
  sscanf(options[6],  "%*[^:]%*[:] %c",  &(foo));
  if      (foo == 'y') doInject_i = true;
  else if (foo == 'n') doInject_i = false;
  else {
    printf("Error in UniformRestricted::UniformRestricted(): doInject_i has to be either 'y' or 'n' \n");
    exit(1);
  }
  sscanf(options[7],  "%*[^:]%*[:] %c",  &(foo));
  if      (foo == 'y') doInject_n = true;
  else if (foo == 'n') doInject_n = false;
  else {
    printf("Error in UniformRestricted::UniformRestricted(): doInject_n has to be either 'y' or 'n' \n");
    exit(1);
  }
  
  sscanf(options[8], "%*[^:]%*[:] %lg", &Tinj);

}

void UniformRestricted::init() {
  //Sanity checks
  if (this->maxR > picConfig.Rmax) {
    printf("Error in UniformRestricted::init(): maxR=%g greater than the grid size Rmax %g\n", this->maxR, picConfig.Rmax);
    exit(1);
  }
  if (this->maxZ > picConfig.Zmax) {
    printf("Error in UniformRestricted::init(): maxZ=%g greater than the grid size Zmax %g\n", this->maxZ, picConfig.Zmax);
    exit(1);
  }
  if (this->minR < 0.0) {
    printf("Error in UniformRestricted::init(): minR=%g < 0.0 \n", this->minR);
    exit(1);
  }
  if (this->minZ < 0.0) {
    printf("Error in UniformRestricted::init(): minZ=%g < 0.0 \n", this->minZ);
    exit(1);
  }
  if (this->maxR <= this->minR) {
    printf("Error in UniformRestricted::init(): maxR=%g <= minR=%g \n", this->maxR, this->minR);
    exit(1);
  }
  if (this->maxZ <= this->minZ) {
    printf("Error in UniformRestricted::init(): maxZ=%g <= minZ=%g \n", this->maxZ, this->minZ);
    exit(1);
  }
  
  //Volume to inject in:
  vol = PI*(SQU(maxR)-SQU(minR))*(maxZ-minZ); //[dZ^3]
  //Debye length (also calculated in init_reactions(), but init_reactions() runs later)
  this->Ldb = 7.43e2*sqrt(picConfig.T_ref/picConfig.n_ref);
  vol *= (this->Ldb)*(this->Ldb)*(this->Ldb)*picConfig.dz*picConfig.dz*picConfig.dz; //[cm^3]
  
  //Superparticle ratio (also in init_reactions(), but init_reactions() runs later)
  this->N_sp = picConfig.n_ref*(this->Ldb)*(this->Ldb)*(this->Ldb)/picConfig.Ndb;
  
  //Number of particles to inject pr. species
  num_inject = (size_t) (density*vol/(this->N_sp));
}

void UniformRestricted::print_par() {
  printf( " - density         %g [cm^-3] \n", density );
  printf( " - maxR            %g [dZ]\n",     maxR );
  printf( " - maxZ            %g [dZ]\n",     maxZ );
  printf( " - minR            %g [dZ]\n",     minR );
  printf( " - minZ            %g [dZ]\n",     minZ );
  printf( " - doInject_e      %c \n",         doInject_e ? 'y' : 'n' );
  printf( " - doInject_i      %c \n",         doInject_i ? 'y' : 'n' );
  printf( " - doInject_n      %c \n",         doInject_n ? 'y' : 'n' );
  printf( " - Tinj            %g [eV] \n",    Tinj );
  printf( " - Ldb             %g [cm]\n",     this->Ldb );
  printf( " - Vol             %g [cm^3]\n",   vol );
  printf( " - N_sp            %g \n",         this->N_sp );
  printf( " - num_inject      %zu \n",        num_inject );
}
void UniformRestricted::inject_e(ParticleSpecies* pa) {
  double vInj = sqrt(Tinj/picConfig.T_ref)*picConfig.v_te;
  if (doInject_e) injector(pa, vInj);
}
void UniformRestricted::inject_n(ParticleSpecies* pa) {
  if ( pa->charge != 0. ) {
    printf ( "Error in UniformRestricted::inject_n(): Species '%s' is not a neutral!", pa->name );
    exit(1);
  }

  double vInj = dt_ion * sqrt(1/M_ions[Lastion]) * sqrt(Tinj/picConfig.T_ref)*picConfig.v_te;
  if (doInject_n) injector(pa, vInj);
}
void UniformRestricted::inject_i(ParticleSpecies* pa) {
  if ( pa->charge == 0. ) {
    printf ( "Error in UniformRestricted::inject_i(): Species '%s' is not an ion!", pa->name );
    exit(1);
  }
  
  double vInj = dt_ion * sqrt(1/pa->mass) * sqrt(Tinj/picConfig.T_ref)*picConfig.v_te;
  if (doInject_i) injector(pa, vInj);
}

void UniformRestricted::injector(ParticleSpecies* pa, double vInj) {
    pa->ReserveSpace(num_inject);
  
  for (size_t i = 0; i < num_inject; i++) {
    //Position: Uniform within requested area (swept into volume)
    pa->z.push_back( minZ + (maxZ-minZ)*RAND );
    pa->r.push_back( sqrt((maxR*maxR-minR*minR)*RAND + minR*minR) );
    
    //Velocity: Inject thermal gas
    pa->vz.push_back( GausRandom(0, vInj) );
    pa->vr.push_back( GausRandom(0, vInj) );
    pa->vt.push_back( GausRandom(0, vInj) );
    
    //Stuff
    pa->m.push_back( 1 );
  }
}
