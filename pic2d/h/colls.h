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

  colls.h:
  Header file for collisions.cpp

***********************************************************************/

#ifndef COLLISIONS_H
#define COLLISIONS_H

#include "ParticleSpecies.h"

void coll_el_knm_2D( ParticleSpecies* pa,
		     int nr,int nz,int NZ,
		     double Mpa_over_me, int kind,
		     Vec3d *momcheck, double *engcheck, int ncoll);

void coll_ion_neutral_noSP_2D( ParticleSpecies* neutrals, double M_n,
			       ParticleSpecies* ions, double M_i,
			       int nr, int nz, int NZ, Reaction React,
			       Vec3d *momcheck, double *engcheck   );

#define COLL_N_N_2D_OMP_MINPARTICLES 1000 //Lower bound for where paralellization kicks in
void coll_n_n_2D( ParticleSpecies* neutrals,
		  int nr, int nz, int NZ, Reaction React,
		  Vec3d *momcheck, double *engcheck );

void coll_el_all_fake_2D( ParticleSpecies* molecules, double M_m,   // molecules
			  ParticleSpecies* electrons,              // electrons
			  int nr, int nz, int NZ, Reaction React);

void coll_el_neutrals_2D( ParticleSpecies* neutrals, double M_n,
			  ParticleSpecies* electrons,
			  ParticleSpecies* ions, int nr, int nz, int NZ, Reaction React,
			  Vec3d *momcheck, double *engcheck );

#endif
