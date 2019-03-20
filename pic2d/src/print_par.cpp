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

  print_par.cpp:
  Prints read-in parameters to output file

***********************************************************************/

#include  "pic.h"
#include  "dim.h"

#include  <time.h>
#include  <stdio.h>
#include  <math.h>

#include "H5Cpp.h"

#define   XTRN  extern

#include  "var.h"
#include  "outp.h"
#include  "arrays1.h"



void print_parameters_2D( void ) {

  printf( "\n" );
  printf( "Input parameters initialised: \n" );
  printf( "\n" );
  printf( "_______________________________________________________________\n" );
  printf( " \n");
  printf( "Scaling and main parameters: \n" );
  printf( " - n_ref:               %g (in 1/cm^3)\n", n_ref );
  printf( " - T_ref:               %g (in eV)\n", T_ref );
  printf( " - Ndb:                 %g\n", Ndb );
  printf( " - nr:                  %d\n", nr );
  printf( " - nz:                  %d\n", nz );
  printf( " - dz:                  %g (in L_Db)\n", dz );
  printf( " - dt:                  %g (in O_pe^-1)\n", Omega_pe );
  printf( " - Ti_over_Te:          %g\n", Ti_over_Te );
  printf( " - mi_over_me:          %g\n", mi_over_me );
  //printf( " - seed:                %llu ( Note: may be changed by initrand() )\n", seed );
  printf( " - RNGbaseSeed          %lu\n", RNGbaseSeed );
  printf( " - numParaThreads       %i\n", numParaThreads );

  printf( " \n" );
  printf( "Timesteps: \n" );
  printf( " - dt_ion:              %d (in dt)\n", dt_ion );
  printf( " - ncoll_el:            %d (in dt)\n", ncoll_el );
  printf( " - ncoll_ion:           %d (in dt)\n", ncoll_ion );
  printf( " - dt_diagn:            %d (in dt)\n", dt_diagn );
  printf( " - nstepsmin:           %d\n", nstepsmin );
  printf( " - nstepsmax:           %d\n", nstepsmax );
  printf( " - nav_start:           %d (in dt)\n", nav_start );
  printf( " - nav_time:            %d (in dt)\n", nav_time );
  printf( " - nav_dt:              %d (in dt)\n", nav_dt );

  printf( " \n" );
  printf( "Fields, particles and boundary conditions: \n" );
  printf( " - Bz_ext:              %g\n", Bz_ext );
  printf( " - Bt_ext:              %g\n", Bt_ext );
  printf( " - e2inj_step:          %d (in dt)\n", e2inj_step );
  printf( " - n2inj_step:          %d (in dt)\n", n2inj_step );
  printf( " - i2inj_step:          %d (in dt)\n", i2inj_step );

  printf("\n");
  printf("Particle boundary parameters:\n");
  printf(" - ArcName: %s\n", pbounds->getName());
  pbounds->print_par();

  printf("\n");
  printf("Circuit parameters:\n");
  printf(" - CircuitName: %s\n", circuit->getName());
  circuit->print_par();

  printf("\n");
  printf("Initial particle distribution:\n");
  if (iParts == NULL) {
    printf(" -- No initial particle distribution initialized --\n");
  }
  else {
    printf(" - Initial particle distribution name: %s\n", iParts->getName());
    iParts->print_par();
  }

  printf( " \n" );
  printf( "Velocities: \n" );
  printf( " - v_te:                %g (in dt/dz)\n", v_te );
  printf( " - v_ti:                %g (in dt_ion/dz)\n", v_ti );
  printf( " - cs:                  %g (in dt/dz)\n", cs );
  printf( " - vi_0:                %g (in dt_ion/dz)\n", vi_0 );

  printf( " \n" );
  printf( "Options (0=yes 1=no): \n" );
  printf( " - Output coordinates:  %d \n", OUT_COORD );
  printf( " - Output vdf:          %d \n", OUT_VDF );
  printf( " - Magnetic push:       %d \n", MAGNETIC );
  printf( " - Continuing old run:  %d \n", CONTINUATION );
  printf( " - Binary output files: %d \n", BINARY_OUTPUT );
  printf( " - Enable collissions:  %c \n", DOCOLL ? 'y' : 'n');
  printf( " - Enable debugging:    %c \n", DODEBUG ? 'y' : 'n');


  printf( "\n" );
  printf( " - Field boundary condition BC = %i \n", BC );
  if (BC == 0)
    printf( "   (Phi=0 at r=nr aka infinity)\n");
  else if (BC == 1)
    printf( "   (DPhi/Dr=0 in r=nr aka infty)\n");
  else if (BC == 2)
    printf( "   (Phi const at r-boundaries, DPhi/Dz=0 at electrodes)\n");
  else if (BC == 3)
    printf( "   (Periodic B.C.)\n");
  else if (BC == 4)
    printf( "   (DPhi/Dr=0 in r=nr aka infty, alternative implementation)\n");

  printf( "\n" );
  printf( "_______________________________________________________________\n" );
  printf( " \n");
  fflush(stdout);

}

void outputfile_addParameterMetadata(H5::H5File* outputFile, const int nsteps) {
  // Adds a new group with metadata useful for scaling

  H5::Group group_metadata       = outputFile->createGroup("/METADATA");
  H5::Group group_metadata_input = outputFile->createGroup("/METADATA/INPUTFILE");
  H5::Group group_metadata_calc  = outputFile->createGroup("/METADATA/CALCULATED");
  H5::Group group_metadata_dynamic = outputFile->createGroup("/METADATA/DYNAMIC");

  H5::DataSpace dataspace_scalar(H5S_SCALAR);

  // -- INPUTFILE METADATA --
  // Attribute: Reference density [cm^-3]
  H5::Attribute attribute_n_ref = group_metadata_input.createAttribute("n_ref", H5::PredType::NATIVE_DOUBLE, dataspace_scalar);
  attribute_n_ref.write(H5::PredType::NATIVE_DOUBLE, &n_ref);

  // Attribute: Reference temperature [eV]
  H5::Attribute attribute_T_ref = group_metadata_input.createAttribute("T_ref", H5::PredType::NATIVE_DOUBLE, dataspace_scalar);
  attribute_T_ref.write(H5::PredType::NATIVE_DOUBLE, &T_ref);

  // Attribute: Particles/Debye cube [-]
  H5::Attribute attribute_Ndb = group_metadata_input.createAttribute("Ndb", H5::PredType::NATIVE_DOUBLE, dataspace_scalar);
  attribute_Ndb.write(H5::PredType::NATIVE_DOUBLE, &Ndb);

  // Attribute: Number of cells (r) [-]
  H5::Attribute attribute_nr = group_metadata_input.createAttribute("nr", H5::PredType::NATIVE_INT, dataspace_scalar);
  attribute_nr.write(H5::PredType::NATIVE_INT, &nr);

  // Attribute: Number of cells (z) [-]
  H5::Attribute attribute_nz = group_metadata_input.createAttribute("nz", H5::PredType::NATIVE_INT, dataspace_scalar);
  attribute_nz.write(H5::PredType::NATIVE_INT, &nz);

  // Attribute: grid size dz [Debyes]
  H5::Attribute attribute_dz = group_metadata_input.createAttribute("dz", H5::PredType::NATIVE_DOUBLE, dataspace_scalar);
  attribute_dz.write(H5::PredType::NATIVE_DOUBLE, &dz);

  // Attribute: Timestep dt [Omega_pe^-1]
  H5::Attribute attribute_Omega_pe = group_metadata_input.createAttribute("Omega_pe", H5::PredType::NATIVE_DOUBLE, dataspace_scalar);
  attribute_Omega_pe.write(H5::PredType::NATIVE_DOUBLE, &Omega_pe);

  // -- CALCULATED METADATA --

  // Attribute: Particles/superparticle ratio [-]
  H5::Attribute attribute_N_sp = group_metadata_calc.createAttribute("N_sp", H5::PredType::NATIVE_DOUBLE, dataspace_scalar);
  attribute_N_sp.write(H5::PredType::NATIVE_DOUBLE, &N_sp);

  // Attribute: Debye length [cm]
  H5::Attribute attribute_Ldb = group_metadata_calc.createAttribute("Ldb", H5::PredType::NATIVE_DOUBLE, dataspace_scalar);
  attribute_Ldb.write(H5::PredType::NATIVE_DOUBLE, &Ldb);

  // Attribute: Plasma frequency [s^-1]
  double O_pe = 56414.6*sqrt(n_ref);
  H5::Attribute attribute_O_pe = group_metadata_calc.createAttribute("O_pe", H5::PredType::NATIVE_DOUBLE, dataspace_scalar);
  attribute_O_pe.write(H5::PredType::NATIVE_DOUBLE, &O_pe);

  // Attribute: Time step [s]
  double dT = Omega_pe/O_pe;
  H5::Attribute attribute_dT = group_metadata_calc.createAttribute("dT", H5::PredType::NATIVE_DOUBLE, dataspace_scalar);
  attribute_dT.write(H5::PredType::NATIVE_DOUBLE, &dT);

  // Attribute: Grid size [cm]
  double dZ = Ldb*dz;
  H5::Attribute attribute_dZ = group_metadata_calc.createAttribute("dZ", H5::PredType::NATIVE_DOUBLE, dataspace_scalar);
  attribute_dZ.write(H5::PredType::NATIVE_DOUBLE, &dZ);

  // Attribute: Total grid size R [cm]
  double R = dZ*nr;
  H5::Attribute attribute_R = group_metadata_calc.createAttribute("R", H5::PredType::NATIVE_DOUBLE, dataspace_scalar);
  attribute_R.write(H5::PredType::NATIVE_DOUBLE, &R);

  // Attribute: Total grid size Z [cm]
  double Z = dZ*nz;
  H5::Attribute attribute_Z = group_metadata_calc.createAttribute("Z", H5::PredType::NATIVE_DOUBLE, dataspace_scalar);
  attribute_Z.write(H5::PredType::NATIVE_DOUBLE, &Z);

  // -- DYNAMIC METADATA --

  // Attribute: Current time in the simulation [ns]
  double simTime = nsteps*Omega_pe*1e9/(56414.6*sqrt(n_ref));
  H5::Attribute attribute_simTime = group_metadata_dynamic.createAttribute("simTime", H5::PredType::NATIVE_DOUBLE, dataspace_scalar);
  attribute_simTime.write(H5::PredType::NATIVE_DOUBLE, &simTime);

  // Attribute: Step number in the simulation
  H5::Attribute attribute_nsteps = group_metadata_dynamic.createAttribute("Nsteps", H5::PredType::NATIVE_INT, dataspace_scalar);
  attribute_nsteps.write(H5::PredType::NATIVE_INT, &nsteps);

}
