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

  var.h:
  Global variables

***********************************************************************/

#ifndef VAR_H
#define VAR_H

#include "circuit.h"
#include "arcbounds.h"
#include "initialParticles.h"
#include "ArcPicConfig.h"

XTRN ArcPicConfig picConfig;

// New-style classes for handling different parts of the physics
XTRN Circuit* circuit;
XTRN ArcBounds* pbounds;
XTRN InitialParticles* iParts;

XTRN int diagn_start,
  dt_diagn,
  dt_ion,
  n_aver,
  n_aver_diagn,
  n_aver_ion,
  nav_start,
  nav_time,
  nav_dt,
  ncoll_el, ncoll_ion,
  nr, nz, //Number of cells
  NR, NZ, //Number of grids (nr+=1)
  nsteps, nstepsmax, nstepsmin;


XTRN int CONTINUATION, BC;

XTRN bool OUT_COORD;     // Outputting particle coordinates?
XTRN bool OUT_EFIELD;    // Outputting electric field?
XTRN bool OUT_VDF;       // Outputting VDF?
XTRN bool MAGNETIC;      // Magnetic push?
XTRN bool BINARY_OUTPUT; // Write binary output files (HDF5)?
XTRN bool DOCOLL;        // Enable collissions?
XTRN bool DODEBUG;       // Enable extra debugging checks?

XTRN int e2inj_step, i2inj_step, n2inj_step;

XTRN double En_e, En_f, En_p, En_tot; //En_i, En_n,

XTRN Reaction  React_Cu_el, React_Cu_ion, React_Cup_Cu_el, React_Cu_Cu;
XTRN double Ldb, N_sp;

//RNG and paralellization setup
XTRN unsigned long int RNGbaseSeed;
XTRN int numParaThreads;

#endif
