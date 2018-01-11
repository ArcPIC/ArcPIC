#include "ParticleSpecies.h"

#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <vector>

ParticleSpecies::ParticleSpecies(size_t nr, size_t nz, std::string name, double m_over_me, double q) :
  nr(nr), nz(nz), name(name), m_over_me(m_over_me), q(q) {
  ordcount = new size_t[(nr+1)*(nz+1)];
  temP = new std::vector<size_t>[(nr+1)*(nz+1)];
  densMap = new double[(nr+1)*(nz+1)];
  this->ZeroDensityMap();
}

ParticleSpecies::~ParticleSpecies(){
  delete[] temP;     temP     = NULL;
  delete[] ordcount; ordcount = NULL;
  delete[] temP;     temP     = NULL;
  delete[] densMap;  densMap  = NULL;
}

void ParticleSpecies::ResizeDelete(const size_t newSize) {
  if (newSize > GetN()) {
    std::cout << "Error when resizing ParticleSpecies; newSize > oldSize!" << std::endl;
    exit(1);
  }
  z.resize(newSize);
  r.resize(newSize);
  vz.resize(newSize);
  vr.resize(newSize);
  vt.resize(newSize);
  m.resize(newSize);
}
void ParticleSpecies::ResizeDeleteBy(const size_t nDelete){
  if (nDelete > GetN()){
    std::cout << "Error when downsizing ParticleSpecies; nDelete > oldSize!" << std::endl;
    exit(1);
  }
  ResizeDelete(GetN()-nDelete);
}

void ParticleSpecies::Order2D() {
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

void ParticleSpecies::UpdateDensityMap( double V_cell[] ) {
  // Initialise
  this->ZeroDensityMap();
  
  // Positive ions: only this
  for (size_t n=0; n<this->GetN(); n++) {
    double hr = r[n];
    int j  = (int)hr;
    hr -= j;
    
    double hz = z[n];
    int k  = (int)hz;
    hz -= k;
    
    densMap[j*(nz+1) + k]         += (1-hr)*(1-hz)/V_cell[j];
    densMap[(j+1)*(nz+1) + k]     += hr*(1-hz)/V_cell[j+1];
    densMap[j*(nz+1) + (k+1)]     += (1-hr)*hz/V_cell[j];
    densMap[(j+1)*(nz+1) + (k+1)] += hr*hz/V_cell[j+1];
  }
  
  /*  Multipying by particles charge	   */
  if (this->q != 0.0 ) {
    for(size_t j=0; j<(nr+1)*(nz+1); j++) {
      densMap[j] *= this->q;
    }
  }
}

void ParticleSpecies::ZeroDensityMap() {
  for (size_t j=0; j<(nr+1); j++) {
    for (size_t k=0; k<(nz+1); k++) {
      densMap[j*(nz+1)+k] = 0.;
    }
  }
}

Particle ParticleSpecies::GetOneParticle(size_t n) {
  Particle pa;
  pa.p.z  = this->z[n];
  pa.p.r  = this->r[n];
  pa.p.vz = this->vz[n];
  pa.p.vr = this->vr[n];
  pa.p.vt = this->vt[n];
  pa.p.m  = this->m[n];
  return pa;
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
