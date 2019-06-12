/**********************************************************************

  2D3V PIC-MCC CODE '2D Arc-PIC'

  Copyright 2010-2019 CERN, University of Oslo, and Helsinki Institute of Physics.
  This software is distributed under the terms of the
  GNU General Public License version 3 (GPL Version 3),
  copied verbatim in the file LICENCE.md. In applying this
  license, CERN does not waive the privileges and immunities granted to it
  by virtue of its status as an Intergovernmental Organization
  or submit itself to any jurisdiction.

  Project website: http://arcpic.web.cern.ch/
  Developers: Helga Timko, Kyrre Sjobak

  ArcPicConfig.h:
    The purpose of the ArcPicConfig class is to have one consistent place to store
    the formerly global configuration variables controlling things like scaling,
    number of time steps, flags, etc.

***********************************************************************/

#ifndef ARCPICCONFIG_H
#define ARCPICCONFIG_H

#include "H5Cpp.h"

class ArcPicConfig {

 public:
    double Ampl,
    av_start, av_time,
    Bt_ext, Bz_ext,
    cs,
    dr, dz,
    dt_out,
    lambda_De,
    me_over_mi,
    mi_over_me,
    Ndb,
    n_ref, T_ref,       /*  Ref. dens. and temp. for rescaling  */
    Omega_pe,
    qe, qi,
    Ti_over_Te,
    v_te, v_ti,
    vi_0,
    Zmin, Zmax, Rmin, Rmax;

    void print_parameters_2D( void );

    void outputfile_addParameterMetadata(H5::H5File* outputFile, const int nsteps);

    //from init.cpp
    void calc_parameters_2D( void );
    void re_init( void );
    void init_reactions( void );

};

#endif
