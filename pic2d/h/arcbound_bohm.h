/**********************************************************************

  2D3V PIC-MCC CODE '2D Arc-PIC'

  Copyright 2010-2018 CERN and Helsinki Institute of Physics.
  This software is distributed under the terms of the
  GNU General Public License version 3 (GPL Version 3),
  copied verbatim in the file LICENCE.md. In applying this
  license, CERN does not waive the privileges and immunities granted to it
  by virtue of its status as an Intergovernmental Organization
  or submit itself to any jurisdiction.

  Project website: http://arcpic.web.cern.ch/
  Developers: Helga Timko, Kyrre Sjobak, Riccardo Casagrande

  arcbound_bohm.h:
  Particle boundary conditions, implements a Bohm sheet condition

***********************************************************************/

#ifndef ARCBOUND_BOHM_H
#define ARCBOUND_BOHM_H

#include "arcbounds.h"

#include <vector>

class ArcBohm : public ArcRemover {
  //Implements a Bohm sheet condition on the anode.
  // Can be used stand-alone or as a mother class for a more specialized model
 public:
  ArcBohm(std::vector<char*>& options);
  virtual ~ArcBohm();

  void init(unsigned int nr, double zmin, double zmax, double rmax);
  void re_init(unsigned int nr, double zmin, double zmax, double rmax);

  virtual void print_par() const;

  virtual void inject_e(ParticleSpecies* pa, double const Ez[]) {}; //NOP
  virtual void inject_i(ParticleSpecies* pa, double const Ez[], unsigned int sort);
  virtual void inject_n(ParticleSpecies* pa, double const Ez[]) {}; //NOP

  //taken care of by ArcRemover
  //virtual void remove_i(ParticleSpecies* pa, unsigned int sort);
  //virtual void remove_n(ParticleSpecies* pa);
  //virtual void remove_e();

  virtual void timestep(unsigned int nstep, bool isOutputTimestep);

  //Write to the arcbounds_original.dat output file
  virtual void writeFile_extras(unsigned int nstep) {
    //writeFile_arcboundsOriginalDat(nstep);
  }

  //Save and restore backup data
  virtual void backup(FILE* file) {};
  virtual void restoreBackup(FILE* file) {}; 

  virtual const char* getName() const { return "ArcOriginalNewHS"; };

 private:
  //Settings from input file
  double v_injAnode_ion;  // Injection velocity from the anode
  double v_injAnode_temp; // Injection temperature from the anode
  double rate_injAnode;   // Number of ions to inject per timstep per area
  unsigned int maxR_injAnode; // Maximum radius to inject

};


#endif
