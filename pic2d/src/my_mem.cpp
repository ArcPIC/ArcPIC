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

  my_mem.cpp:
  Memory (de-)allocation

***********************************************************************/

#include  <slu_ddefs.h>
#include  <iostream>

#include  "pic.h"
#include  "dim.h"

#define XTRN extern
#include  "var.h"
#include  "arrays1.h"
#undef XTRN

void allocate_arrays( int nr, int nz, int** perm_c, int** perm_r, double** rhs ) {
  
  // Allocate particle arrays
  try {
    elec =  new Particle[NPART];
    ions =  new Particle[NSpecies*NPART];
  }
  catch (const std::bad_alloc&) {
    printf("Error in allocate_arrays():\n");
    printf("Could not allocate the requested memory for particle arrays\n");
    printf("Most likely, NPART is too high.\n");
    exit(1);
  }

  //Only allocate VDF arrays if they are needed
  if ( OUT_VDF == 0 ) {
    vdf_ez   = new double[Nvdst*(nr+1)/2*(nz+1)/2];
    vdf_er   = new double[Nvdst*(nr+1)/2*(nz+1)/2];
    vdf_eabs = new double[Nvdst*(nr+1)/2*(nz+1)/2];
    vdf_iz   = new double[Nvdst*(nr+1)/2*(nz+1)/2];
    vdf_ir   = new double[Nvdst*(nr+1)/2*(nz+1)/2];
    vdf_iabs = new double[Nvdst*(nr+1)/2*(nz+1)/2];
    vdf_nz   = new double[Nvdst*(nr+1)/2*(nz+1)/2];
    vdf_nr   = new double[Nvdst*(nr+1)/2*(nz+1)/2];
    vdf_nabs = new double[Nvdst*(nr+1)/2*(nz+1)/2];
  }
  else {
    vdf_ez   = NULL;
    vdf_er   = NULL;
    vdf_eabs = NULL;
    vdf_iz   = NULL;
    vdf_ir   = NULL;
    vdf_iabs = NULL;
    vdf_nz   = NULL;
    vdf_nr   = NULL;
    vdf_nabs = NULL;
  }
  
  //Various matrices
  mom_el  = new Moments[(nr+1)*(nz+1)];
  mom_ion = new Moments[NSpecies*(nr+1)*(nz+1)];

  Vcell = new double[(nr+1)];

  // Densities
  n_e = new double[(nr+1)*(nz+1)];
  n_i = new double[NSpecies*(nr+1)*(nz+1)];
  n_e_av = new double[(nr+1)*(nz+1)];
  n_i_av = new double[NSpecies*(nr+1)*(nz+1)];

  // EM-fields
  phi = new double[(nr+1)*(nz+1)];
  phi_av = new double[(nr+1)*(nz+1)];
  E_grid_r = new double[(nr+1)*(nz+1)];
  E_grid_z = new double[(nr+1)*(nz+1)];
  E_av_r = new double[(nr+1)*(nz+1)];
  E_av_z = new double[(nr+1)*(nz+1)];
  E_ion_r = new double[Lastion*(nr+1)*(nz+1)];
  E_ion_z = new double[Lastion*(nr+1)*(nz+1)];

  diagn_Te = new double[(nr+1)*(nz+1)];
  diagn_ve = new double[(nr+1)*(nz+1)];
  diagn_ne = new double[(nr+1)*(nz+1)];
  diagn_dens = new double[(nr+1)*(nz+1)];

  // Ordering arrays for collisions
  e_order = new size_t[(nr+1)*(nz+1)];
  i_order = new size_t[NSpecies*(nr+1)*(nz+1)];
  
  // Energy outputting
  En_i = new double[NSpecies];
  
  //Field solver matrices for SuperLU
  int ni = (nr+1)*(nz+1);
  int nnz = 0;
  if      ( BC == 0 ) nnz = 2*(nr+1) + (nz-1) + 4*(nz-1) + 5*(nr-1)*(nz-1);
  else if ( BC == 1 ) nnz = 2*(nr+1) + 3*(nz-1) + 4*(nz-1) + 5*(nr-1)*(nz-1);
  else if ( BC == 2 ) nnz = 3*(nr-1) + 3*(nr-1) + (nz+1) + (nz+1) + 5*(nr-1)*(nz-1);
  else if ( BC == 3 ) nnz = 5*(nr+1)*(nz+1);
  else if ( BC == 4 ) nnz = 2*(nr+1) + 8*(nz-1) + 5*(nr-1)*(nz-1);
  else                printf("Error with solver B.C.s \n");
  
  printf ("In allocate_arrays, ni= %d, nnz= %d \n", ni, nnz);
  
  if ( !(*perm_r = intMalloc(ni)) ) ABORT("Malloc fails for perm_r[].");
  if ( !(*perm_c = intMalloc(ni)) ) ABORT("Malloc fails for perm_c[].");
  if ( !(*rhs = doubleMalloc(nnz)) ) ABORT("Malloc fails for rhs[].");

  int i;
  for (i=0; i<ni; i++) {
      (*perm_r)[i]=0;
      (*perm_c)[i]=0;
  }
  printf("In allocate, perm_c %d \n", (*perm_c)[0]);
}


void delete_arrays( int** perm_c, int** perm_r, double** rhs ) {
  
  delete[] elec;
  delete[] ions;
  
  SUPERLU_FREE (*perm_r);
  SUPERLU_FREE (*perm_c);
  SUPERLU_FREE (*rhs);
  
}


