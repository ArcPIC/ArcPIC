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

  outputz.cpp:
  Output and diagnostic routines

***********************************************************************/

#include <cstdio>

#include <iostream>

#include <string>
#include <sstream>
#include <iomanip>

#include "H5Cpp.h"

#include "dim.h"
#include "pic.h"
#include "mydef.h"
#include "outputz.h"

/***************************
 Outputting routines
***************************/

void out_fv_along_2D( double fvz[], double fvr[], double fvabs[], int nr, int nz,
		      const char *dat1, const char *dat2, const char *dat3 ) {
  FILE *file1, *file2, *file3;
  int j, k, n;
  int halfnr = nr/2;
  int halfnz = nz/2;
  
  file1 = fopen(dat1, "w");
  file2 = fopen(dat2, "w");
  file3 = fopen(dat3, "w");

  for (j=0; j<halfnr; j++) {
    for (k=0; k<halfnz; k++) {
      
      for (n=0; n < Nvdst; n++) {
	fprintf(file1, "% .5e ", fvz[n + Nvdst*(j*(nz/2)+k)]);
	fprintf(file2, "% .5e ", fvr[n + Nvdst*(j*(nz/2)+k)]);
	fprintf(file3, "% .5e ", fvabs[n + Nvdst*(j*(nz/2)+k)]);
      }
      
      fprintf(file1, "\n");
      fprintf(file2, "\n");
      fprintf(file3, "\n");
    }
  }
  
  fclose(file1);
  fclose(file2);
  fclose(file3);
  
}

void out_fv_along_2D_h5( double fvz[], double fvr[], double fvabs[], int nr, int nz,
			 const char* const tablename_z, const char* const tablename_r, const char* const tablename_abs, H5::Group& group_veldist ) {
  int halfnr = nr/2;
  int halfnz = nz/2;
  
  //Create the dataspaces for the file
  hsize_t dims_file[] = {halfnr,halfnz,Nvdst};
  H5::DataSpace dataspace_z  (3, dims_file);
  H5::DataSpace dataspace_r  (3, dims_file);
  H5::DataSpace dataspace_abs(3, dims_file);
  
  H5::DataSet dataset_z   = group_veldist.createDataSet(tablename_z, H5::PredType::NATIVE_DOUBLE, dataspace_z);
  H5::DataSet dataset_r   = group_veldist.createDataSet(tablename_r, H5::PredType::NATIVE_DOUBLE, dataspace_r);
  H5::DataSet dataset_abs = group_veldist.createDataSet(tablename_abs, H5::PredType::NATIVE_DOUBLE, dataspace_abs);

  //Create a window into the file dataspace (1 row)
  hsize_t dims_mem[] = {1,1,Nvdst};
  H5::DataSpace dataspace_mem(3,dims_mem,dims_file);

  const hsize_t offset_memInFile[] = {0,0,0};
  dataspace_z.selectHyperslab  (H5S_SELECT_SET, dims_mem, offset_memInFile);
  dataspace_r.selectHyperslab  (H5S_SELECT_SET, dims_mem, offset_memInFile);
  dataspace_abs.selectHyperslab(H5S_SELECT_SET, dims_mem, offset_memInFile);
  
  for (ssize_t j=0; j<halfnr; j++) {
    for (ssize_t k=0; k<halfnz; k++) {
      double databuffer_z   [Nvdst];
      double databuffer_r   [Nvdst];
      double databuffer_abs [Nvdst];
      for (ssize_t n=0; n < Nvdst; n++) {
	databuffer_z[n]   = fvz[n + Nvdst*(j*(nz/2)+k)];
	databuffer_r[n]   = fvr[n + Nvdst*(j*(nz/2)+k)];
	databuffer_abs[n] = fvabs[n + Nvdst*(j*(nz/2)+k)];
      }
      const hssize_t offset_memInFile_shift[] = {j,k,0};
      double databuffer[] = {1,2,3};
      dataspace_z.offsetSimple(offset_memInFile_shift);
      dataset_z.write(databuffer_z,H5::PredType::NATIVE_DOUBLE, dataspace_mem, dataspace_z);
      dataset_z.write(databuffer,H5::PredType::NATIVE_DOUBLE, dataspace_mem, dataspace_z);
      
      dataspace_r.offsetSimple(offset_memInFile_shift);
      dataset_r.write(databuffer_r,H5::PredType::NATIVE_DOUBLE, dataspace_mem, dataspace_r);
      dataspace_abs.offsetSimple(offset_memInFile_shift);
      dataset_abs.write(databuffer_abs,H5::PredType::NATIVE_DOUBLE, dataspace_mem, dataspace_abs);
    }
  }
}


void out_dens_2D( double dens_av[], int n_aver, double sign, int nr, int nz, int NZ, double Omega_pe,
		  double dr, double dz, const char *dat_p ) {
  FILE *file_p;
  double fp;
  int j, k;
  
  if (n_aver < 1)  n_aver = 1;

  fp = sign*SQU(1./Omega_pe)/n_aver;

  file_p = fopen(dat_p, "w");

  for (j=0; j<=nr; j++) {
    for (k=0; k<=nz; k++) {
      fprintf(file_p, "% .5e % .5e % .5e \n", dr*j, dz*k, fp*dens_av[j*NZ+k]);
    }
  }
  
  fclose(file_p);
  
}

void out_dens_2D_h5( double dens_av[], int n_aver, double sign, int nr, int nz, int NZ, double Omega_pe,
		     double dr, double dz, const char* const tablename, H5::Group& group_dens ) {
  
  if (n_aver < 1)  n_aver = 1;

  double fp = sign*SQU(1./Omega_pe)/n_aver;

  // HDF5 is internally row-major (i.e. C convention)
  hsize_t dims_file[] = {nr+1,nz+1,3};
  H5::DataSpace dataspace_file(3, dims_file);
  H5::DataSet dataset = group_dens.createDataSet(tablename, H5::PredType::NATIVE_DOUBLE, dataspace_file);
  
  //Create a window into the file dataspace (1 row)
  hsize_t dims_mem[] = {1,1,3};
  H5::DataSpace dataspace_mem(3,dims_mem,dims_file);
  
  const hsize_t offset_memInFile[] = {0,0,0};
  dataspace_file.selectHyperslab(H5S_SELECT_SET, dims_mem, offset_memInFile);
  
  for (ssize_t j=0; j<=nr; j++) {
    for (ssize_t k=0; k<=nz; k++) {
      double databuffer[] = { dr*j, dz*k, fp*dens_av[j*NZ+k] };
      
      const hssize_t offset_memInFile_shift[] = {j,k,0};
      dataspace_file.offsetSimple(offset_memInFile_shift);
      dataset.write(databuffer,H5::PredType::NATIVE_DOUBLE, dataspace_mem, dataspace_file);
    }
  }
}


void out_phi_2D( double phi[], int n_aver, int nr, int nz, int NZ,
		 double Omega_pe, double dr, double dz,
		 const char *dat1 ) {

  FILE *file1;
  double f1;
  int j, k;
  
  if( n_aver < 1 )  n_aver = 1;
  
  f1 = SQU(dz/Omega_pe)/n_aver;

  file1 = fopen(dat1, "w");
  for (j=0; j<=nr; j++) {
    for (k=0; k<=nz; k++) {
      fprintf(file1, "% .5e % .5e % .5e \n", dr*j, dz*k, f1*phi[j*NZ+k]);
    }
  }
  fclose(file1);
  
}

void out_phi_2D_h5( double phi[], int n_aver, int nr, int nz, int NZ,
		    double Omega_pe, double dr, double dz,
		    const char* const tablename, H5::Group& group_emfield ) {
  
  if( n_aver < 1 )  n_aver = 1;
  
  double f1 = SQU(dz/Omega_pe)/n_aver;

  //Create the dataspace for the file
  hsize_t dims_file[] = {nr+1,nz+1,3};
  H5::DataSpace dataspace_file(3, dims_file);
  H5::DataSet dataset = group_emfield.createDataSet(tablename, H5::PredType::NATIVE_DOUBLE, dataspace_file);
  //Create a window into the file dataspace (1 row)
  hsize_t dims_mem[] = {1,1,3};
  H5::DataSpace dataspace_mem(3,dims_mem,dims_file);

  const hsize_t offset_memInFile[] = {0,0,0};
  dataspace_file.selectHyperslab(H5S_SELECT_SET, dims_mem, offset_memInFile);
  
  for (ssize_t j=0; j<=nr; j++) {
    for (ssize_t k=0; k<=nz; k++) {
      double databuffer[] = { dr*j, dz*k, f1*phi[j*NZ+k] };

      const hssize_t offset_memInFile_shift[] = {j,k,0};
      dataspace_file.offsetSimple(offset_memInFile_shift);
      dataset.write(databuffer,H5::PredType::NATIVE_DOUBLE, dataspace_mem, dataspace_file);
    }
  }
}


void out_efield_2D( double efield_z[], double efield_r[], int n_aver, int nr, int nz, int NZ,
		    double Omega_pe, double dr, double dz, const char *dat1, const char *dat2 ) {

  FILE *file1, *file2;
  double f1;
  int j, k;
  
  if( n_aver < 1 )  n_aver = 1;
  
  f1 = 2*dz/SQU(Omega_pe)/n_aver;

  file1 = fopen(dat1, "w");
  file2 = fopen(dat2, "w");

  for (j=0; j<=nr; j++) {
    for (k=0; k<=nz; k++) {
      fprintf(file1, "% .5e % .5e % .5e \n", dr*j, dz*k, f1*efield_z[j*NZ+k]);
      fprintf(file2, "% .5e % .5e % .5e \n", dr*j, dz*k, f1*efield_r[j*NZ+k]);
    }
  }
  fclose(file1);
  fclose(file2);
  
}

void out_efield_2D_h5( double efield_z[], double efield_r[], int n_aver, int nr, int nz, int NZ,
		       double Omega_pe, double dr, double dz,
		       const char* const tablename, H5::Group& group_emfield ) {
  
  if( n_aver < 1 )  n_aver = 1;
  
  double f1 = 2*dz/SQU(Omega_pe)/n_aver;
  
  //Create the dataspaces for the file
  hsize_t dims_file[] = {nr+1,nz+1,4};
  H5::DataSpace dataspace(3, dims_file);

  H5::DataSet dataset = group_emfield.createDataSet(tablename, H5::PredType::NATIVE_DOUBLE, dataspace);
  
  //Create a window into the file dataspace (1 row)
  hsize_t dims_mem[] = {1,1,4};
  H5::DataSpace dataspace_mem(3,dims_mem,dims_file);

  const hsize_t offset_memInFile[] = {0,0,0};
  dataspace.selectHyperslab(H5S_SELECT_SET, dims_mem, offset_memInFile);
  
  for (ssize_t j=0; j<=nr; j++) {
    for (ssize_t k=0; k<=nz; k++) {
      double databuffer[] = { dr*j, dz*k, f1*efield_z[j*NZ+k], f1*efield_r[j*NZ+k] };

      const hssize_t offset_memInFile_shift[] = {j,k,0};

      dataspace.offsetSimple(offset_memInFile_shift);
      dataset.write(databuffer,H5::PredType::NATIVE_DOUBLE, dataspace_mem, dataspace);
    }
  }
}

void out_vels_2D( Moments  mom[], int nr, int nz, int NZ, double u0,
		  double dr, double dz, const char *dat1, const char *dat2, const char *dat3 ) {

  FILE *file1, *file2, *file3;
  int j, k;
  
  u0 = MAX(u0,1.e-10);

  for (j=0; j<nr; j++) {
    for (k=0; k<nz; k++) {
      
      if( mom[j*NZ+k].n > 1 ) {
	mom[j*NZ+k].uz /= (u0*mom[j*NZ+k].n);
	mom[j*NZ+k].ur /= (u0*mom[j*NZ+k].n);
	mom[j*NZ+k].ut /= (u0*mom[j*NZ+k].n);
      }
      else {
	mom[j*NZ+k].uz = mom[j*NZ+k].ur = mom[j*NZ+k].ut = 0.;
      }
    }
  }
  
  file1 = fopen(dat1, "w");
  file2 = fopen(dat2, "w");
  file3 = fopen(dat3, "w");
  
  for (j=0; j<nr; j++) {
    for (k=0; k<nz; k++) {
      fprintf(file1, "% .5e % .5e % .5e \n", dr*j, dz*k, mom[j*NZ+k].uz);
      fprintf(file2, "% .5e % .5e % .5e \n", dr*j, dz*k, mom[j*NZ+k].ur);
      fprintf(file3, "% .5e % .5e % .5e \n", dr*j, dz*k, mom[j*NZ+k].ut);
    }
  }
  
  fclose(file1);
  fclose(file2);
  fclose(file3);
  
}

void out_vels_2D_h5( Moments  mom[], int nr, int nz, int NZ, double u0,
		     double dr, double dz,
		     const char* const tablename1, const char* const tablename2, const char* const tablename3, H5::Group& group_velavg ) {
  
  u0 = MAX(u0,1.e-10);

  for (size_t j=0; j<nr; j++) {
    for (size_t k=0; k<nz; k++) {
      if( mom[j*NZ+k].n > 1 ) {
	mom[j*NZ+k].uz /= (u0*mom[j*NZ+k].n);
	mom[j*NZ+k].ur /= (u0*mom[j*NZ+k].n);
	mom[j*NZ+k].ut /= (u0*mom[j*NZ+k].n);
      }
      else {
	mom[j*NZ+k].uz = mom[j*NZ+k].ur = mom[j*NZ+k].ut = 0.;
      }
    }
  }

  //Create the dataspaces for the file
  hsize_t dims_file[] = {nr,nz,3};
  H5::DataSpace dataspace_file1(3, dims_file);
  H5::DataSpace dataspace_file2(3, dims_file);
  H5::DataSpace dataspace_file3(3, dims_file);
  
  H5::DataSet dataset1 = group_velavg.createDataSet(tablename1, H5::PredType::NATIVE_DOUBLE, dataspace_file1);
  H5::DataSet dataset2 = group_velavg.createDataSet(tablename2, H5::PredType::NATIVE_DOUBLE, dataspace_file2);
  H5::DataSet dataset3 = group_velavg.createDataSet(tablename3, H5::PredType::NATIVE_DOUBLE, dataspace_file3);
  
  //Create a window into the file dataspace (1 row)
  hsize_t dims_mem[] = {1,1,3};
  H5::DataSpace dataspace_mem(3,dims_mem,dims_file);

  const hsize_t offset_memInFile[] = {0,0,0};
  dataspace_file1.selectHyperslab(H5S_SELECT_SET, dims_mem, offset_memInFile);
  dataspace_file2.selectHyperslab(H5S_SELECT_SET, dims_mem, offset_memInFile);
  dataspace_file3.selectHyperslab(H5S_SELECT_SET, dims_mem, offset_memInFile);
  
  for (ssize_t j=0; j<nr; j++) {
    for (ssize_t k=0; k<nz; k++) {
      double databuffer1[] = { dr*j, dz*k, mom[j*NZ+k].uz };
      double databuffer2[] = { dr*j, dz*k, mom[j*NZ+k].ur };
      double databuffer3[] = { dr*j, dz*k, mom[j*NZ+k].ut };

      const hssize_t offset_memInFile_shift[] = {j,k,0};

      dataspace_file1.offsetSimple(offset_memInFile_shift);
      dataset1.write(databuffer1,H5::PredType::NATIVE_DOUBLE, dataspace_mem, dataspace_file1);
      
      dataspace_file2.offsetSimple(offset_memInFile_shift);
      dataset2.write(databuffer2,H5::PredType::NATIVE_DOUBLE, dataspace_mem, dataspace_file2);
      
      dataspace_file3.offsetSimple(offset_memInFile_shift);
      dataset3.write(databuffer3,H5::PredType::NATIVE_DOUBLE, dataspace_mem, dataspace_file3);
      
    }
  }
}


void out_temps_2D( Moments  mom[], double u0, double fnorm, int nr, int nz, int NZ,
		   double dr, double dz,
		   const char *dat1, const char *dat2, const char *dat3 ) {

  FILE *file1, *file2, *file3;
  int j, k;
  
  for (j=0; j<nr; j++) {
    for (k=0; k<nz; k++) {
      
      if( mom[j*NZ+k].n > 1 ) {
	mom[j*NZ+k].tz = ( mom[j*NZ+k].tz/(u0*u0*mom[j*NZ+k].n) - SQU(mom[j*NZ+k].uz) )*fnorm;
	mom[j*NZ+k].tr = ( mom[j*NZ+k].tr/(u0*u0*mom[j*NZ+k].n) - SQU(mom[j*NZ+k].ur) )*fnorm;
	mom[j*NZ+k].tt = ( mom[j*NZ+k].tt/(u0*u0*mom[j*NZ+k].n) - SQU(mom[j*NZ+k].ut) )*fnorm;
      }
      else {
	mom[j*NZ+k].tz = mom[j*NZ+k].tr = mom[j*NZ+k].tt = 0.;
      }
      
    }
  }

  file1 = fopen(dat1, "w");
  file2 = fopen(dat2, "w");
  file3 = fopen(dat3, "w");


  for (j=0; j<nr; j++) {
    for (k=0; k<nz; k++) {
      fprintf(file1, "% .5e % .5e % .5e \n", dr*j, dz*k, mom[j*NZ+k].tz);
      fprintf(file2, "% .5e % .5e % .5e \n", dr*j, dz*k, mom[j*NZ+k].tr);
      fprintf(file3, "% .5e % .5e % .5e \n", dr*j, dz*k, mom[j*NZ+k].tt);
    }
  }
  
  fclose(file1);
  fclose(file2);
  fclose(file3);
}

void out_temps_2D_h5( Moments  mom[], double u0, double fnorm, int nr, int nz, int NZ,
		      double dr, double dz,
		      const char* const tablename1, const char* const tablename2, const char* const tablename3, H5::Group& group_temp ) {
  
  for (size_t j=0; j<nr; j++) {
    for (size_t k=0; k<nz; k++) {
      
      if( mom[j*NZ+k].n > 1 ) {
	mom[j*NZ+k].tz = ( mom[j*NZ+k].tz/(u0*u0*mom[j*NZ+k].n) - SQU(mom[j*NZ+k].uz) )*fnorm;
	mom[j*NZ+k].tr = ( mom[j*NZ+k].tr/(u0*u0*mom[j*NZ+k].n) - SQU(mom[j*NZ+k].ur) )*fnorm;
	mom[j*NZ+k].tt = ( mom[j*NZ+k].tt/(u0*u0*mom[j*NZ+k].n) - SQU(mom[j*NZ+k].ut) )*fnorm;
      }
      else {
	mom[j*NZ+k].tz = mom[j*NZ+k].tr = mom[j*NZ+k].tt = 0.;
      }
      
    }
  }
  
  //Create the dataspaces for the file
  hsize_t dims_file[] = {nr,nz,3};
  H5::DataSpace dataspace_file1(3, dims_file);
  H5::DataSpace dataspace_file2(3, dims_file);
  H5::DataSpace dataspace_file3(3, dims_file);
  
  H5::DataSet dataset1 = group_temp.createDataSet(tablename1, H5::PredType::NATIVE_DOUBLE, dataspace_file1);
  H5::DataSet dataset2 = group_temp.createDataSet(tablename2, H5::PredType::NATIVE_DOUBLE, dataspace_file2);
  H5::DataSet dataset3 = group_temp.createDataSet(tablename3, H5::PredType::NATIVE_DOUBLE, dataspace_file3);
  
  //Create a window into the file dataspace (1 row)
  hsize_t dims_mem[] = {1,1,3};
  H5::DataSpace dataspace_mem(3,dims_mem,dims_file);

  const hsize_t offset_memInFile[] = {0,0,0};
  dataspace_file1.selectHyperslab(H5S_SELECT_SET, dims_mem, offset_memInFile);
  dataspace_file2.selectHyperslab(H5S_SELECT_SET, dims_mem, offset_memInFile);
  dataspace_file3.selectHyperslab(H5S_SELECT_SET, dims_mem, offset_memInFile);

  for (ssize_t j=0; j<nr; j++) {
    for (ssize_t k=0; k<nz; k++) {

      double databuffer1[] = { dr*j, dz*k, mom[j*NZ+k].tz };
      double databuffer2[] = { dr*j, dz*k, mom[j*NZ+k].tr };
      double databuffer3[] = { dr*j, dz*k, mom[j*NZ+k].tt };

      const hssize_t offset_memInFile_shift[] = {j,k,0};

      dataspace_file1.offsetSimple(offset_memInFile_shift);
      dataset1.write(databuffer1,H5::PredType::NATIVE_DOUBLE, dataspace_mem, dataspace_file1);
      dataspace_file2.offsetSimple(offset_memInFile_shift);
      dataset2.write(databuffer2,H5::PredType::NATIVE_DOUBLE, dataspace_mem, dataspace_file2);

      dataspace_file3.offsetSimple(offset_memInFile_shift);
      dataset3.write(databuffer3,H5::PredType::NATIVE_DOUBLE, dataspace_mem, dataspace_file3);
    }
  }
}

void out_coords_2D( ParticleSpecies* pa, int fnorm,
		    double omega_pe, double dz, const char *dat ) {
  FILE *file;
  double vnorm = dz/omega_pe/fnorm;
  
  file = fopen(dat, "w");
  
  for (size_t i=0; i < pa->GetN(); i++)  {
    fprintf(file, "% .5e % .5e % .5e % .5e % .5e \n", (pa->z[i])*dz, (pa->r[i])*dz,
	    (pa->vz[i])*vnorm, (pa->vr[i])*vnorm, (pa->vt[i])*vnorm);
  }
  
  fclose(file);
  
}

void out_coords_2D_h5( ParticleSpecies* pa, int fnorm,
		       double omega_pe, double dz,
		       const char* const tablename, H5::Group& group_coords ) {
  
  double vnorm = dz/omega_pe/fnorm;

  //std::cout << "out_coords_2D_h5 tablename="<<tablename<<std::endl;
  
  // HDF5 is internally row-major (i.e. C convention)
  hsize_t dims_file[2] = {pa->GetN(),5};
  H5::DataSpace dataspace_file(2, dims_file);
  H5::DataSet dataset = group_coords.createDataSet(tablename, H5::PredType::NATIVE_DOUBLE, dataspace_file);

  if (pa->GetN()==0) {
    //Don't try to write to an empty dataset
    return;
  }
  
  //Create a window into the file dataspace (1 row)
  hsize_t dims_mem[] = {1,5};
  H5::DataSpace dataspace_mem(2,dims_mem,dims_file);

  const hsize_t offset_memInFile[] = {0,0};
  dataspace_file.selectHyperslab(H5S_SELECT_SET, dims_mem, offset_memInFile);

  for (ssize_t i=0; i<pa->GetN(); i++)  {
    double databuffer[] = { (pa->z[i])*dz, (pa->r[i])*dz, (pa->vz[i])*vnorm, (pa->vr[i])*vnorm, (pa->vt[i])*vnorm};
    
    const hssize_t offset_memInFile_shift[] = {i,0};
    dataspace_file.offsetSimple(offset_memInFile_shift);
    dataset.write(databuffer,H5::PredType::NATIVE_DOUBLE, dataspace_mem, dataspace_file);
  }
}

// Create the HDF5 file for the current ouput timestep
H5::H5File* createH5File_timestep(const int nsteps, std::string basename){
  
  std::ostringstream timestep_string;
  timestep_string << std::setw(8) << std::setfill('0') << nsteps;
  
  std::string fname = "out/" + basename + "_" + timestep_string.str() + ".h5";
  //H5std_string fname_h5(fname);
  
  //std::cout << "Making HDF5file, timestep="<<nsteps<<" filename="<<fname<<std::endl;
  
  H5::H5File* returnFile = new H5::H5File(fname, H5F_ACC_TRUNC);
  
  return returnFile;
}

/***************************
 Diagnostics
***************************/

void diagn_stability_2D( const double dens[], Particle pa[], double diag_Te[], double diag_ve[], double diag_ne[],
			 double sign, int np, double u0, int nr, int nz, int NR, int NZ, double Omega_pe,
			 double dr, double dz, int steps, int& check, const char *dat_err ) {

  double ne_max = 0., ne_max_old = 0.;
  double ne_av = 0., ne_av_old = 0.; // average over all k's
  int ne_max_j(-1), ne_max_k(-1), ne_av_j(-1);

  // ne/Te ratios
  double noT_max = 0., noT_max_old = 0.;
  double noT_av = 0., noT_av_old = 0.; // average over all k's
  int noT_max_j(-1), noT_max_k(-1), noT_av_j(-1);

  double fne;
  fne = SQU(1./Omega_pe);

  double vz, vr, vt, vz2, vr2, vt2;

  FILE *file;

  // Check ne, ne/Te
  for (int j=0; j<NR; j++) {
    for (int k=0; k<NZ; k++) {
      ne_av_old += sign*dens[j*NZ+k];
      ne_max_old = sign*dens[j*NZ+k];
      
      if ( ne_max < ne_max_old ) {
	ne_max = ne_max_old;
	ne_max_j = j;
	ne_max_k = k;
      }
    }
    if ( ne_av < ne_av_old ) {
      ne_av = ne_av_old;
      ne_av_j = j;
    }
    ne_av_old = 0.;
  }
  ne_av  *= fne/NZ; // in units of n_ref
  ne_max *= fne;    // in units of n_ref
  
  for (int j=0; j<NR; j++) {
    for (int k=0; k<NZ; k++) {
      diag_Te[j*NZ+k]=0.;
      diag_ve[j*NZ+k]=0.;
      diag_ne[j*NZ+k]=0.;
    }
  }

  for (int n=0; n<np; n++) {
    int j = (int)pa[n].p.r;
    int k = (int)pa[n].p.z;
    
    vz = pa[n].p.vz;
    vr = pa[n].p.vr;
    vt = pa[n].p.vt;
    diag_ve[j*NZ+k] += vz + vr + vt;
    
    vz2 = SQU(vz);
    vr2 = SQU(vr);
    vt2 = SQU(vt);
    diag_Te[j*NZ+k] += vz2 + vr2 + vt2;
    
    diag_ne[j*NZ+k]++;
  }
  
  for (int j=0; j<nr; j++) {
    for (int k=0; k<nz; k++) {
      diag_Te[j*NZ+k] /= diag_ne[j*NZ+k];
      //diag_Te[j*NZ+k] -= SQU(diag_ve[j*NZ+k]/diag_ne[j*NZ+k]); // don't extract main velocity
      diag_Te[j*NZ+k] /= SQU(u0);
      
      diag_ne[j*NZ+k] = sign*dens[j*NZ+k];
      
      noT_av_old += diag_ne[j*NZ+k]/diag_Te[j*NZ+k];
      noT_max_old = diag_ne[j*NZ+k]/diag_Te[j*NZ+k];
      
      if ( noT_max < noT_max_old ) {
	noT_max = noT_max_old;
	noT_max_j = j;
	noT_max_k = k;
      }
    }
    if ( noT_av < noT_av_old ) {
      noT_av = noT_av_old;
      noT_av_j = j;
    }
    noT_av_old = 0.;
  }
  noT_av  *= fne/nz; // in units of n_ref/T_ref
  noT_max *= fne;    // in units of n_ref/T_ref
  
  // Write to error file if needed
  file = fopen(dat_err, "a");
  // If error => put check == 1 => stop.
  if ( ne_av > 5. ) {
    fprintf(file, "***ERROR*** GLOBAL(%d steps): <ne>/n_ref= %.5e along j= %d \n", steps, ne_av, ne_av_j);
    check = 1;
  }
  if ( noT_av > 2. ) {
    fprintf(file, "***ERROR*** GLOBAL(%d steps): <ne/Te>/(n_ref/T_ref)= %.5e along j= %d \n", steps, noT_av, noT_av_j);
    check = 1;
  }
  
  // WARNINGS
  if ( ne_av > 1. ) {
    fprintf(file, "***WARNING*** GLOBAL(%d steps): <ne>/n_ref= %.5e along j= %d \n", steps, ne_av, ne_av_j);
  }
  if ( ne_max > 5. ) {
    fprintf(file, "***WARNING*** LOCAL(%d steps): [ne]max/n_ref= %.5e at j= %d k= %d \n", steps, ne_max, ne_max_j, ne_max_k);
  }
  if ( noT_av > 1. ) {
    fprintf(file, "***WARNING*** GLOBAL(%d steps): <ne/Te>/(n_ref/T_ref)= %.5e along j= %d \n", steps, noT_av, noT_av_j);
  }
  if ( noT_max > 2. ) {
    fprintf(file, "***WARNING*** LOCAL(%d steps): [ne/Te]max/(n_ref/T_ref)= %.5e at j= %d k= %d \n", steps, noT_max, noT_max_j, noT_max_k);
  }
  
  fclose(file);
}

void diagn_av_stability( const double dens[], double diag_Te[], double diag_ne[], int n_av,
			 double sign, double u0, int nr, int nz, int NR, int NZ, double Omega_pe,
			 double dr, double dz, int steps, int *check, const char *dat_err ) {

  double ne_max = 0., ne_max_old = 0.;
  double ne_av = 0., ne_av_old = 0.; // average over all k's
  int ne_max_j(-1), ne_max_k(-1), ne_av_j(-1);

  // ne/Te ratios
  double noT_max = 0., noT_max_old = 0.;
  double noT_av = 0., noT_av_old = 0.; // average over all k's
  int noT_max_j(-1), noT_max_k(-1), noT_av_j(-1);

  if (n_av < 0) n_av = 1;
  double fne = SQU(1./Omega_pe)/n_av;

  FILE *file;

  // Check ne, ne/Te
  for (int j=0; j<NR; j++) {
    for (int k=0; k<NZ; k++) {
      ne_av_old += sign*dens[j*NZ+k];
      ne_max_old = sign*dens[j*NZ+k];
      
      if ( ne_max < ne_max_old ) {
	ne_max = ne_max_old;
	ne_max_j = j;
	ne_max_k = k;
      }
    }
    if ( ne_av < ne_av_old ) {
      ne_av = ne_av_old;
      ne_av_j = j;
    }
    ne_av_old = 0.;
  }
  ne_av *= fne/NZ; // in units of n_ref
  ne_max *= fne; // in units of n_ref
  
  for (int j=0; j<nr; j++) {
    for (int k=0; k<nz; k++) {
      diag_Te[j*NZ+k] /= diag_ne[j*NZ+k]; // don't extract main velocity
      diag_Te[j*NZ+k] /= SQU(u0);
      
      if ( diag_Te[j*NZ+k] > 1e-20 ) {
	noT_av_old += sign*dens[j*NZ+k]/diag_Te[j*NZ+k];
	noT_max_old = sign*dens[j*NZ+k]/diag_Te[j*NZ+k];
      }
      
      if ( noT_max < noT_max_old ) {
	noT_max = noT_max_old;
	noT_max_j = j;
	noT_max_k = k;
      }
    }
    if ( noT_av < noT_av_old ) {
      noT_av = noT_av_old;
      noT_av_j = j;
    }
    noT_av_old = 0.;
  }
  noT_av  *= fne/nz; // in units of n_ref/T_ref
  noT_max *= fne;    // in units of n_ref/T_ref

  // Write to error file if needed
  file = fopen(dat_err, "a");

  // IF ERROR => put check == 1 => stop.
  if ( ne_av > 5. ) {
    fprintf(file, "***ERROR*** GLOBAL(%d steps): <ne>/n_ref= %.5e along j= %d \n", steps, ne_av, ne_av_j);
    *check = 1;
  }
  if ( noT_av > 2. ) {
    fprintf(file, "***ERROR*** GLOBAL(%d steps): <ne/Te>/(n_ref/T_ref)= %.5e along j= %d \n", steps, noT_av, noT_av_j);
    *check = 1;
  }

  // WARNINGS
  if ( ne_av > 1. ) {
    fprintf(file, "***WARNING*** GLOBAL(%d steps): <ne>/n_ref= %.5e along j= %d \n", steps, ne_av, ne_av_j);
  }
  if ( ne_max > 5. ) {
    fprintf(file, "***WARNING*** LOCAL(%d steps): [ne]max/n_ref= %.5e at j= %d k= %d \n", steps, ne_max, ne_max_j, ne_max_k);
  }
  if ( noT_av > 1. ) {
    fprintf(file, "***WARNING*** GLOBAL(%d steps): <ne/Te>/(n_ref/T_ref)= %.5e along j= %d \n", steps, noT_av, noT_av_j);
  }
  if ( noT_max > 2. ) {
    fprintf(file, "***WARNING*** LOCAL(%d steps): [ne/Te]max/(n_ref/T_ref)= %.5e at j= %d k= %d \n", steps, noT_max, noT_max_j, noT_max_k);
  }
  
  fclose(file);
  
}
 

