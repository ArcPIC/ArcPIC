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

  arrays1.h:
  Define external arrays

***********************************************************************/

#ifndef ARRAYS1_H
#define ARRAYS1_H

XTRN Particle *elec;
XTRN Particle *ions;
XTRN Particle *temp;

XTRN const char *Names[NSpecies];   // Name of ions
XTRN size_t nr_i[NSpecies];         // Nr of ions
XTRN double q_ions[NSpecies];       // Charge of ions
XTRN double M_ions[NSpecies];       // Mass of ions
XTRN double cs_ions[NSpecies];      // Sound velocity
XTRN double vt_ions[NSpecies];      // Thermal velocity

XTRN Moments *mom_el, *mom_ion;

XTRN double* Vcell;

// Densities
XTRN double *n_e, *n_i, *n_e_av, *n_i_av;

// EM-fields
XTRN double *phi, *phi_av;
XTRN double *E_grid_r, *E_grid_z;
XTRN double *E_av_r, *E_av_z;
XTRN double *E_ion_r, *E_ion_z;

// Stability diagnosis
XTRN double *diagn_Te, *diagn_ve,
  *diagn_ne, *diagn_dens;

// Ordering arrays for collisions
XTRN size_t *e_order, *i_order;

// Energy outputting
XTRN double *En_i;

// VDF arrays (only allocated if VDF is activated)
XTRN double* vdf_ez;
XTRN double* vdf_er;
XTRN double* vdf_eabs;
XTRN double* vdf_iz;
XTRN double* vdf_ir;
XTRN double* vdf_iabs;
XTRN double* vdf_nz;
XTRN double* vdf_nr;
XTRN double* vdf_nabs;

#endif
