#include "ParticleSpecies.h"

#include <stdio.h>
#include <stdlib.h>

#include <vector>

ParticleSpecies::ParticleSpecies(size_t nr, size_t nz, std::string name, double m_over_me, double q) :
  nr(nr), nz(nz), name(name), m_over_me(m_over_me), q(q) {
  ordcount = new size_t[nr*nz];
  temP = new std::vector<size_t>[nr*nz];
}

ParticleSpecies::~ParticleSpecies(){
  delete[] temP;
  delete[] ordcount;
}

const int ParticleSpecies::ExpandBy(size_t n_expand) {
  size_t N = GetN();
  z.resize(N+n_expand);
  r.resize(N+n_expand);
  vz.resize(N+n_expand);
  vr.resize(N+n_expand);
  vt.resize(N+n_expand);
  m.resize(N+n_expand);
  return GetN();
}

void ParticleSpecies::ReserveSpace(size_t n){
  z.reserve(n);
  r.reserve(n);
  vz.reserve(n);
  vr.reserve(n);
  vt.reserve(n);
  m.reserve(n);
}

void ParticleSpecies::order_2D() {
  // Zero ordcount
  for (size_t j=0; j<nr; j++) {
    for (size_t k=0; k<nz; k++){
      ordcount[j*(nz+1) + k]= 0;
    }
  }
  
  //Sort particles
  for (size_t n=0; n<GetN(); n++) {
    int ir = (int)r[n];
    int iz = (int)z[n];

    //Sanity check -- we should be inside the grid
    if ( (ir<0) || (ir>=nr) || (iz<0) || (iz>=nz) ) {
      //Edge cases
      if (r[n] == (double) nr) {
	ir--;
      }
      else if (z[n] == (double) nz) {
	iz--;
      }
      else {
	//Over the edge; it's just wrong!
	fprintf(stderr," !!!ERROR!!! %d %d\n", ir, iz);
	fflush(stderr);

	printf("Error in ParticleSpecies::order_2D(): particle position out-of-range, r=%g, z=%g, (ir,iz)= (%d,%d)\n", r[n], z[n], ir, iz);
	
	exit(1);
      }
    }

    temP[ir*nz+iz].push_back(n);
    ordcount[ir*(nz+1) + iz]++;
  }
  
  //Stuff them back into the particle arrays, but now in the correct order
  ShuffleArray(z);
  ShuffleArray(r);
  ShuffleArray(vz);
  ShuffleArray(vr);
  ShuffleArray(vt);
  ShuffleArray(m);
  for (int cidx = 0; cidx < nr*nz; cidx++) {
    temP[cidx].clear();
  }
  
}

void ParticleSpecies::ShuffleArray(std::vector<double> &arr) {
  shuffleTmpArr_dbl.resize(0);
  shuffleTmpArr_dbl.reserve(this->GetN());
  
  for (int cidx = 0; cidx < nr*nz; cidx++) {
    for (size_t l = 0; l < temP[cidx].size(); l++) {
      shuffleTmpArr_dbl.push_back(arr[temP[cidx][l]]);
    }
  }
  
  arr.swap(shuffleTmpArr_dbl);
}
void ParticleSpecies::ShuffleArray(std::vector<int> &arr) {
  shuffleTmpArr_int.resize(0);
  shuffleTmpArr_int.reserve(this->GetN());
  
  for (int cidx = 0; cidx < nr*nz; cidx++) {
    for (size_t l = 0; l < temP[cidx].size(); l++) {
      shuffleTmpArr_int.push_back(arr[temP[cidx][l]]);
    }
  }
  
  arr.swap(shuffleTmpArr_int);
}
