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

  checkbounds.cpp:
  checking read-in data for 2D boundary conditions

***********************************************************************/

#include  "pic.h"

#include  <stdio.h>

#include "checkbounds.h"

bool checkbounds_2D( ParticleSpecies* pa, double rmin, double rmax, double zmin, double zmax ) {

  bool allInside = true;

  for (size_t n=0; n<pa->GetN(); n++) {

    if ( pa->r[n] >= rmax ) {
      printf("Erroneous data: n=%zu  r=%g  rmax=%g\n", n, pa->r[n], rmax);
      pa->r[n] = rmax - 1.e-10;
      allInside = false;
    }
    else if ( pa->r[n] < rmin ) {
      printf("Erroneous data: n=%zu  r=%g  rmin=%g\n", n, pa->r[n], rmin);
      pa->r[n] = rmin + 1.e-10;
      allInside = false;
    }

    if ( pa->z[n] >= zmax ) {
      printf("Erroneous data: n=%zu  z=%g  zmax=%g\n", n, pa->z[n], zmax);
      pa->z[n] = zmax - 1.e-10;
      allInside = false;
    }
    else if ( pa->z[n] < zmin ) {
      printf("Erroneous data: n=%zu  z=%g  zmin=%g\n", n, pa->z[n], zmin);
      pa->z[n] = zmin + 1.e-10;
      allInside = false;
    }

  }

  fflush(stdout);

  return allInside;
}

