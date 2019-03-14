/**********************************************************************

  2D3V PIC-MCC CODE '2D Arc-PIC'

  Copyright 2010-2018 CERN, Helsinki Institute of Physics, and University of Oslo.
  This software is distributed under the terms of the
  GNU General Public License version 3 (GPL Version 3),
  copied verbatim in the file LICENCE.md. In applying this
  license, CERN does not waive the privileges and immunities granted to it
  by virtue of its status as an Intergovernmental Organization
  or submit itself to any jurisdiction.

  Project website: http://arcpic.web.cern.ch/
  Developers: Helga Timko, Kyrre Sjobak

  outputz.h:
  Header file for outputz.cpp

***********************************************************************/

#include "ParticleSpecies.h"

void out_fv_along_2D( double fvz[], double fvr[], double fvabs[], int nr, int nz,
                      const char *dat1, const char *dat2, const char *dat3 );
void out_fv_along_2D_h5( double fvz[], double fvr[], double fvabs[], int nr, int nz,
                         const char* const tablename_z, const char* const tablename_r, const char* const tablename_abs,
                         H5::Group& group_veldist );

void out_dens_2D( double dens_av[], int n_aver, double sign,
                  int nr, int nz, int NZ, double Omega_pe,
                  double dr, double dz,
                  const char *dat_p );
void out_dens_2D_h5( double dens_av[], int n_aver, double sign,
                     int nr, int nz, int NZ, double Omega_pe,
                     double dr, double dz,
                     const char* const tablename, H5::Group& group_dens );

void out_phi_2D( double phi[], int n_aver, int nr, int nz, int NZ,
                 double Omega_pe, double dr, double dz,
                 const char *dat1 );
void out_phi_2D_h5( double phi[], int n_aver, int nr, int nz, int NZ,
                    double Omega_pe, double dr, double dz,
                    const char* const tablename, H5::Group& group_emfield );


void out_efield_2D( double efield_z[], double efield_r[], int n_aver, int nr, int nz, int NZ,
                    double Omega_pe, double dr, double dz, const char *dat1, const char *dat2 );
void out_efield_2D_h5( double efield_z[], double efield_r[], int n_aver, int nr, int nz, int NZ,
                       double Omega_pe, double dr, double dz,
                       const char* const tablename, H5::Group& group_emfield );

void out_vels_2D( Moments  mom[], int nr, int nz, int NZ, double u0,
                  double dr, double dz,
                  const char *dat1, const char *dat2, const char *dat3 );
void out_vels_2D_h5( Moments  mom[], int nr, int nz, int NZ, double u0,
                     double dr, double dz,
                     const char* const tablename1, const char* const tablename2, const char* const tablename3,
                     H5::Group& group_velocities );


void out_temps_2D( Moments  mom[], double u0, double fnorm, int nr, int nz, int NZ,
                   double dr, double dz, const char *dat1, const char *dat2, const char *dat3 );
void out_temps_2D_h5( Moments  mom[], double u0, double fnorm, int nr, int nz, int NZ,
                      double dr, double dz,
                      const char* const tablename1, const char* const tablename2, const char* const tablename3,
                      H5::Group& group_temp );


void out_coords_2D( ParticleSpecies* pa, int fnorm,
                    double omega_pe, double dz,
                    const char *dat );
void out_coords_2D_h5( ParticleSpecies* pa, int fnorm,
                       double omega_pe, double dz,
                       const char* const tablename, H5::Group& group_coords );

H5::H5File* createH5File_timestep( const int nsteps, std::string basename = "output" );

void diagn_stability_2D( const double dens[], Particle pa[], double diag_Te[], double diag_ve[], double diag_ne[],
                         double sign, int np, double u0, int nr, int nz, int NR, int NZ, double Omega_pe,
                         double dr, double dz, int steps, int& check, const char *dat_err );

void diagn_av_stability( const double dens[], double diag_Te[], double diag_ne[], int n_av,
                         double sign, double u0, int nr, int nz, int NR, int NZ, double Omega_pe,
                         double dr, double dz, int steps, int *check, const char *dat_err );
