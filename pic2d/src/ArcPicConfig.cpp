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

  ArcPicConfig.cpp:
    The purpose of the ArcPicConfig class is to have one consistent place to store
    the formerly global configuration variables controlling things like scaling,
    number of time steps, flags, etc.

***********************************************************************/

#include "ArcPicConfig.h"

#include  "pic.h"
#include  "dim.h"

#include  <time.h>
#include  <stdio.h>
#include  <math.h>

#include <sys/stat.h>
#include <iostream>

#include "H5Cpp.h"

#define   XTRN  extern
#include  "var.h"
#include  "outp.h"
#include  "arrays1.h"
#include "circuit.h"
#include "mydef.h"
#undef XTRN


void ArcPicConfig::print_parameters_2D( void ) {

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
  printf( " - Output coordinates:     %c \n", OUT_COORD     ? 'y' : 'n' );
  printf( " - Output electric field:  %c \n", OUT_EFIELD    ? 'y' : 'n' );
  printf( " - Output vdf:             %c \n", OUT_VDF       ? 'y' : 'n' );
  printf( " - Magnetic push:          %c \n", MAGNETIC      ? 'y' : 'n' );
  printf( " - Continuing old run:     %d \n", CONTINUATION );
  printf( " - Binary output files:    %c \n", BINARY_OUTPUT ? 'y' : 'n' );
  printf( " - Enable collissions:     %c \n", DOCOLL        ? 'y' : 'n' );
  printf( " - Enable debugging:       %c \n", DODEBUG       ? 'y' : 'n' );

  printf("\n");
  printf( " - RNGbaseSeed          %lu\n", RNGbaseSeed );
  printf( " - numParaThreads       %i\n", numParaThreads );

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

void ArcPicConfig::outputfile_addParameterMetadata(H5::H5File* outputFile, const int nsteps) {
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


void ArcPicConfig::calc_parameters_2D( void ) {

  Zmin = 0.;
  Zmax = (double)nz;
  Rmin = 0.;
  Rmax = (double)nr;

  NZ = nz + 1;
  NR = nr + 1;

  // Verboncoeur volume factor
  // Charges never placed on cathode/anode, so don't need special corr factor here
  for (int j=0; j<NR; j++) {
    if      (j==0)   Vcell[j] = PI/3;
    else if (j==nr)  Vcell[j] = PI*(nr-1.0/3.0);
    else             Vcell[j] = 2*PI*j;
  }

  // Dimensionless charge
  qe = -SQU(Omega_pe/dz)/(Ndb*dz);
  qi = -qe;        // for Z=1

  // Lengths normalised to dz, timescales to dt
  lambda_De  = 1./dz;

  me_over_mi = 1./mi_over_me;

  v_te = Omega_pe*lambda_De;
  v_ti = v_te*sqrt(Ti_over_Te*me_over_mi);

  cs = sqrt(me_over_mi)*v_te;

  vi_0 = cs;

  // NORMALISE QUANTITIES
  vi_0 *= dt_ion;
  v_ti *= dt_ion;

  // Timesteps
  nstepsmax = (int)(nstepsmax/Omega_pe+0.1);

  nav_dt    = (int)(dt_out/Omega_pe+0.1);
  nav_start = (int)(av_start/Omega_pe+0.1);
  nav_time  = (int)(av_time/Omega_pe+0.1);

  diagn_start = (int)(av_start/Omega_pe+0.1);

  ncoll_ion  = (int)(ncoll_ion*1./dt_ion+0.5)*dt_ion;

  Ampl = (double) ncoll_ion;

  if (nstepsmax < nav_start + nav_time) {
    std::cerr << std::endl;
    std::cerr << "Error in calc_parameters_2D():" << std::endl;
    std::cerr << " Found nstepsmax < nav_start + nav_time" << std::endl
              << " => Simulation will never reach an output timestep." << std::endl;
    std::cerr << " Got:" << std::endl
              << "  nstepsmax = " << nstepsmax << std::endl
              << "  nav_start = " << nav_start << std::endl
              << "  nav_time  = " << nav_time  << std::endl
              << "  Omega_pe  = " << Omega_pe  << std::endl;
    exit(1);
  }
}




void ArcPicConfig::re_init( void ) {
  // Dimensionless charge
  qe = -SQU(Omega_pe/dz)/(Ndb*dz);
  qi = -qe;        // for Z=1

  // Timesteps
  nav_dt = (int)(dt_out/Omega_pe+0.1);
  nav_start += nav_dt;

  nstepsmin = nsteps + 1;

  // External circuit elements
  circuit->re_init();

  //Particle boundary conditions
  pbounds->re_init(nr, Zmin, Zmax, Rmax);

}




void ArcPicConfig::init_reactions( void ) {
  double T_e = T_ref; //was ref_temp; //2500.  ;       //eV
  double n_e = n_ref; //was ref_dens; //1.0e18; //was 1.0e10
  int Npoints = 1000;   // N of points in fit.

  double E;
  double tmp;
  double Ecoeff, CScoeff;
  double Emin, Emax, Eth, Efc;

  //double lambda_b, eps_b, c_b, Eo;

  FILE *outputfile;

  //double a[9];
  double edata[80], Sdata[80];
  int ind;

  Ldb =    7.43e2*sqrt(T_e/n_e);
  CScoeff = n_e*Ldb/(Ndb*PI*SQU(dz));   //remaining 1/(2j+1) is in collisions.cpp //was n_e*Ldb/Ndb in 1D

  N_sp = n_e*Ldb*Ldb*Ldb/Ndb;

  //Check that folder exists, else we get a segfault
  struct stat st;
  if ( stat("xsct", &st) ) {
    printf("init_reactions(): Folder \"xsct\" did not exsist. Creating it now...\n");
    if ( mkdir("xsct", S_IRWXU) ) {
      printf("init_reactions(): Could not create folder \"xsct\". Abort!\n");
      exit(1);
    }
  }

  /*********************************************************************************
   * e- on Cu  elastic  from J.Phys. B Atom Mol Phys 10 (1977) 3323 Trajmar, et al *
   *********************************************************************************/

  edata[0]   =   6.0;        Sdata[0]    =   49.500  ;
  edata[1]   =  10.0;        Sdata[1]    =   83.300  ;
  edata[2]   =  20.0;        Sdata[2]    =   49.400  ;
  edata[3]   =  60.0;        Sdata[3]    =   23.900  ;
  edata[4]   =  100.0;       Sdata[4]    =   16.000  ;

  Emin =  0.0;          //eV
  Eth =   0.0;          //eV
  Emax =  100.0;        //eV
  Efc =   0.;           //eV

  Ecoeff = (Omega_pe/dz)*(Omega_pe/dz)*(1.+M_ions[2])/M_ions[2]/(0.5*T_e);

  React_Cu_el.Emin = Emin*Ecoeff;
  React_Cu_el.Eth =  Eth*Ecoeff;
  React_Cu_el.Emax = Emax*Ecoeff;
  React_Cu_el.Wfc =  sqrt(Efc*Ecoeff);
  React_Cu_el.N = Npoints;
  React_Cu_el.Estep =  (React_Cu_el.Emax -  React_Cu_el.Emin)/Npoints;

  outputfile = fopen("xsct/Cu_el.dat", "w");

  ind = 0;

  for (int i = 0; i <= React_Cu_el.N; i++ ) {
    E = i*React_Cu_el.Estep + React_Cu_el.Emin;

    tmp = 0.;

    if      (E/Ecoeff <= edata[0])   tmp = Sdata[0];
    else if (E/Ecoeff >= edata[4])   tmp = Sdata[4];
    else {
      while (E/Ecoeff > edata[ind] ) ind++;
      tmp = Sdata[ind - 1] + (Sdata[ind]- Sdata[ind-1])*(E/Ecoeff - edata[ind-1])/(edata[ind]- edata[ind-1]);
    }

    React_Cu_el.CS[i] = CScoeff*tmp*1.e-16;
    fprintf(outputfile, "%f  %15.7e \n", E/Ecoeff, React_Cu_el.CS[i]/CScoeff  );
  }

  fclose(outputfile);

  /*************************************************************************
   * e- + Cu -> Cu+ + 2e  ionization                                       *
   * from  Bolorizadeh, et al J Phys B At Mol Opt Phys 27 (1994) 175       *
   *************************************************************************/

  edata[0]   =  7.8    ;    Sdata[0]    =   0.06  ;
  edata[1]   =  8.0    ;    Sdata[1]    =   0.11  ;
  edata[2]   =  8.2    ;    Sdata[2]    =   0.17  ;
  edata[3]   =  8.4    ;    Sdata[3]    =   0.25  ;
  edata[4]   =  8.6    ;    Sdata[4]    =   0.33  ;
  edata[5]   =  9.0    ;    Sdata[5]    =   0.4800;
  edata[6]   =  9.4    ;    Sdata[6]    =   0.67  ;
  edata[7]   =  10.1   ;    Sdata[7]    =   0.93  ;
  edata[8]   =  11.0   ;    Sdata[8]    =   1.28  ;
  edata[9]   =  12.0   ;    Sdata[9]    =   1.66  ;
  edata[10]  =  15.0   ;    Sdata[10]   =   2.28  ;
  edata[11]  =  20.0   ;    Sdata[11]   =   2.8   ;
  edata[12]  =  30.0   ;    Sdata[12]   =   3.21  ;
  edata[13]  =  40.000 ;    Sdata[13]   =   3.08  ;
  edata[14]  =  50.0   ;    Sdata[14]   =   2.9   ;
  edata[15]  =  100.0  ;    Sdata[15]   =   2.55  ;
  edata[16]  =  200.0  ;    Sdata[16]   =   2.05  ;
  edata[17]  =  290.0  ;    Sdata[17]   =   1.63  ;
  edata[18]  =  465.0  ;    Sdata[18]   =   1.25  ;
  edata[19]  =  710.   ;    Sdata[19]   =   0.97  ;
  edata[20]  =  1000.  ;    Sdata[20]   =   0.68  ;
  edata[21]  =  1900.  ;    Sdata[21]   =   0.45  ;
  edata[22]  =  2100.  ;    Sdata[22]   =   0.40  ;

  Emin =  7.8;          //eV
  Eth =   7.8;          //eV
  Emax =  2100.0;        //eV
  Efc =   0.;           //eV

  Ecoeff = (Omega_pe/dz)*(Omega_pe/dz)*(1.+M_ions[2])/M_ions[2]/(0.5*T_e);

  React_Cu_ion.Emin = Emin*Ecoeff;
  React_Cu_ion.Eth =  Eth*Ecoeff;
  React_Cu_ion.Emax = Emax*Ecoeff;
  React_Cu_ion.Wfc =  sqrt(Efc*Ecoeff);
  React_Cu_ion.N = Npoints;
  React_Cu_ion.Estep =  (React_Cu_ion.Emax -  React_Cu_ion.Emin)/Npoints;

  outputfile = fopen("xsct/Cu_ion.dat", "w");

  ind = 0;
  for (int i = 0; i <= React_Cu_ion.N; i++ ) {
    E = i*React_Cu_ion.Estep + React_Cu_ion.Emin;

    tmp = 0.;

    if      (E/Ecoeff <= edata[0])    tmp = Sdata[0];
    else if (E/Ecoeff >= edata[22])   tmp = Sdata[22];
    else {
      while (E/Ecoeff > edata[ind] ) ind++;
      tmp = Sdata[ind - 1] + (Sdata[ind]- Sdata[ind-1])*(E/Ecoeff - edata[ind-1])/(edata[ind]- edata[ind-1]);
    }

    React_Cu_ion.CS[i] = CScoeff*tmp*1.e-16;
    fprintf(outputfile, "%f  %15.7e \n", E/Ecoeff, React_Cu_ion.CS[i]/CScoeff  );
  }

  fclose(outputfile);

  /********************************************************************
   * Cu+  +  Cu elastic                                               *
   * my fit based on Auberton, J. Phys. D, 36 (2003) 1798             *
   ********************************************************************/

  Emin = 0.01;      //eV
  Eth = 0.0;       //eV
  Emax = 1000.;    //eV
  Efc = 0.0;

  Ecoeff = (Omega_pe/dz)*(Omega_pe/dz)*(M_ions[1]+M_ions[2])/M_ions[1]/M_ions[2]/(0.5*T_e);
  Ecoeff *= dt_ion*dt_ion;

  React_Cup_Cu_el.Emin = Emin*Ecoeff;
  React_Cup_Cu_el.Eth =  Eth*Ecoeff;
  React_Cup_Cu_el.Emax = Emax*Ecoeff;
  React_Cup_Cu_el.Wfc =  sqrt(Efc*Ecoeff);

  React_Cup_Cu_el.N = Npoints;
  React_Cup_Cu_el.Estep =  (React_Cup_Cu_el.Emax -  React_Cup_Cu_el.Emin)/React_Cup_Cu_el.N;

  // Now - calculate coefficients 2 go 2 Dim.less units !!
  outputfile = fopen("xsct/Cup_Cu_el.dat", "w");

  for (int i = 0; i <= React_Cup_Cu_el.N; i++ ) {
    E = i*React_Cup_Cu_el.Estep + React_Cup_Cu_el.Emin;

    React_Cup_Cu_el.CS[i] = CScoeff*2.513E-15*pow((E/Ecoeff),-0.1623041);
    fprintf(outputfile, "%f %15.7e \n", E/Ecoeff, React_Cup_Cu_el.CS[i]/CScoeff  );
  }

  fclose(outputfile);

  /*******************************************************************************
   * Cu  +  Cu elastic                                                           *
   * my fit based on CX, Rcu = 1.4A for hard sphere model and ArAr CS asymptotic *
   *******************************************************************************/

  Emin = 0.01;      //eV
  Eth = 0.0;       //eV
  Emax = 1000.;    //eV
  Efc = 0.0;

  Ecoeff = (Omega_pe/dz)*(Omega_pe/dz)*(M_ions[2]+M_ions[2])/M_ions[2]/M_ions[2]/(0.5*T_e);
  Ecoeff *= dt_ion*dt_ion;

  React_Cu_Cu.Emin = Emin*Ecoeff;
  React_Cu_Cu.Eth =  Eth*Ecoeff;
  React_Cu_Cu.Emax = Emax*Ecoeff;
  React_Cu_Cu.Wfc =  sqrt(Efc*Ecoeff);

  React_Cu_Cu.N = Npoints;
  React_Cu_Cu.Estep =  (React_Cu_Cu.Emax -  React_Cu_Cu.Emin)/React_Cu_Cu.N;

  // Now - calculate coefficients 2 go 2 Dim.less units !!
  outputfile = fopen("xsct/Cu_Cu.dat", "w");


  for (int i = 0; i <= React_Cu_Cu.N; i++ ) {
    E = i*React_Cu_Cu.Estep + React_Cu_Cu.Emin;

    React_Cu_Cu.CS[i] = CScoeff*5.524E-16*pow((E/Ecoeff),-0.1623041);
    fprintf(outputfile, "%f %15.7e \n", E/Ecoeff, React_Cu_Cu.CS[i]/CScoeff  );
  }

  fclose(outputfile);

}

