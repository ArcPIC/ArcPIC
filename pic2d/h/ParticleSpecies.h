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

#include "pic.h"

#include <vector>
#include <string>

class ParticleSpecies {
  
 public:
  //Constructor
  ParticleSpecies(size_t nr, size_t nz, std::string name, double mass, double charge);
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
  const size_t GetN() const {
    return z.size();
  };
  //Change the ammount of storage allocated
  const int ExpandBy(const size_t n_expand) {
    size_t N = GetN();
    this->ReserveSpace( N + n_expand );
    return GetN();
  }
  void ReserveSpace(const size_t n) {
    z.reserve(n);
    r.reserve(n);
    vz.reserve(n);
    vr.reserve(n);
    vt.reserve(n);
    m.reserve(n);
  }

  //Copy a particle from position i->j,
  // overwriting the data in position j.
  void CopyParticle(const size_t i, const size_t j) {
    z[j]  = z[i];
    r[j]  = r[i];
    vz[j] = vz[i];
    vr[j] = vr[i];
    vt[j] = vt[i];
    m[j]  = m[i];
  }
  //Resize the particle arrays to something smaller,
  // deleting the extra particles
  void ResizeDelete(const size_t newSize);
  void ResizeDeleteBy(const size_t nDelete);
  
  // Name of the particle species, used for printing, tables, etc.
  const std::string name;

  // Mass of the particle species, relative to the electron mass
  const double mass;
  // Dimensionless particle charge
  const double charge;
  
  //Order the particle arrays by cell
  void Order2D();

  //Update the density map
  void UpdateDensityMap( double V_cell[] );
  void ZeroDensityMap( );
  double* densMap;

  // Extract one particle object
  Particle GetOneParticle(size_t n);
  
 private:
  // Local copy of the nr and nz settings
  const size_t nr, nz;
  
  // Temporary array of vectors used in order_2D to sort the particles;
  // use a member variable so we don't destroy&reallocate the temp array between each use.
  std::vector<size_t> *temP = NULL;
  //Shuffle an array (i.e. z, r, vz, vr, or vt) according to temP;
  std::vector<double> shuffleTmpArr_dbl;
  std::vector<int>    shuffleTmpArr_int;
  void ShuffleArray(std::vector<double> &arr);
  void ShuffleArray(std::vector<int>    &arr);
};

#endif
