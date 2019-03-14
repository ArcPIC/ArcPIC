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

  push.h:
  Header file for push.cpp

***********************************************************************/

#include "ParticleSpecies.h"

void  push_2D( ParticleSpecies* pa, const double Eg_r[], const double Eg_z[],
               int NZ );

void  push_magnetic_onlyBz_2D( ParticleSpecies* pa, const double Eg_r[], const double Eg_z[],
                               const double Bz, double factor, int NZ );

void  push_magnetic_2D( ParticleSpecies* pa, const double Eg_r[], const double Eg_z[],
                        const double Bextz, const double Bextt, double factor, int NZ );

void  push_neutrals_2D( ParticleSpecies* pa );


