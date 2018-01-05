/**********************************************************************

  2D3V PIC-MCC CODE '2D Arc-PIC'

  Copyright 2010-2018 CERN, University of Oslo, and Helsinki Institute of Physics.
  This software is distributed under the terms of the
  GNU General Public License version 3 (GPL Version 3),
  copied verbatim in the file LICENCE.md. In applying this
  license, CERN does not waive the privileges and immunities granted to it
  by virtue of its status as an Intergovernmental Organization
  or submit itself to any jurisdiction.

  Project website: http://arcpic.web.cern.ch/
  Developers: Helga Timko, Kyrre Sjobak

  ParticleSpecies.h:
  Classes for handling the particle arrays

***********************************************************************/

#ifndef ParticleSpecies_h
#define ParticleSpecies_h

#include <vector>

class ParticleSpecies {
  
 public:
  //Constructor
 ParticleSpecies(size_t nr, size_t nz);
  //Destructor
  ~ParticleSpecies();
  
  //Particle arrays, laid out so that they can be vectorized
  std::vector<double> z;
  std::vector<double> r;
  std::vector<double> vz;
  std::vector<double> vr;
  std::vector<double> vt;
  std::vector<int> m;

  //Ordering vector, updated by calling order_2D
  size_t* ordcount;

  //Check the current number of particles
  const int getN() const {
    return z.size();
  };
  

  //Order the particle arrays by cell
  void order_2D();
 private:
  // Local copy of the nr and nz settings
  size_t nr, nz;
  // Temporary array of vectors used for sorting;
  // use a field so we don't destroy&reallocate between each use.
  std::vector<size_t> *temP = NULL;
  //Shuffle an array (i.e. z, r, vz, vr, or vt) according to temP;
  std::vector<double> shuffleTmpArr_dbl;
  std::vector<int>    shuffleTmpArr_int;
  void ShuffleArray(std::vector<double> &arr);
  void ShuffleArray(std::vector<int>    &arr);
};

#endif
