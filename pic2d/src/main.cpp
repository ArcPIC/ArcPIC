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

  main.cpp:
  Main program control and time stepping loop.
  Writing of LOCKFILE, out/timeIndex.dat, mainstats.dat

***********************************************************************/

#include  <stdio.h>
#include <iostream>
#include <iomanip>
#include  <math.h>
#include  <sys/stat.h>
#include  <time.h>
#include  <string.h>
#include  <stdlib.h>

#include  <slu_ddefs.h>

#include <omp.h>

#include "H5Cpp.h"

#include "ArcPicConfig.h"
#include "ParticleSpecies.h"

#define   XTRN
#include  "pic.h"
#include  "var.h"
#include  "dim.h"
#include  "arrays1.h"
#include  "outp.h"
#include  "mydef.h"

#include  "random.h"

#include  "init.h"
#include  "my_mem.h"
#include  "my_time.h"
#include  "phi.h"
#include  "efield.h"
#include  "push.h"
#include  "e_ion.h"
#include  "moms.h"
#include  "aver.h"
#include  "filenames.h"
#include  "outputz.h"
#include  "engy.h"
#include  "colls.h"
#include  "input.h"
#include  "vdf.h"
#include  "backup.h"
#include  "checkbounds.h"
#include  "weightPotential.h"


void print_time( double t );

int main () {

  //Safety: avoid unintended double-starts
  struct stat st;
  if (not stat("LOCKFILE", &st) ) {
    printf("main(): LOCKFILE exists, refusing to start.\n");
    exit(1);
  }
  FILE* lockfile = fopen("LOCKFILE", "w");
  fprintf(lockfile, "This is a lockfile - ArcPic will refuse to start as long as it's present.\n");
  fprintf(lockfile, "Goal: Avoid unintensional restarts which corrupts older data.\n");
  fclose(lockfile);

  // Diagnostics
  int check_stab = 0; // PIC stability ok: 0=yes, 1=no
  int check_dens = 0; // Density well resolved: 0=yes, 1=no
  double f_dens;      // Factor for density check
  if ( snprintf( ferr, LEN_FILENAME, "%s", "error.log") >= LEN_FILENAME ) {
    printf("Error detected when making filename for ferr, generated filename = '%s'\n", ferr);
    exit(1);
  }

  //Particle arrays
  ParticleSpecies* electrons;
  std::vector<ParticleSpecies*> ionSpecies;
  std::vector<ParticleSpecies*> neutralSpecies;

  // SuperLU parameters
  SuperMatrix L_slu, U_slu;
  int *perm_c_slu;
  int *perm_r_slu;
  double *rhs_slu;

  //Time index file
  // (written every output step, used to correlate step number with simulation time
  // for the sake of plotting etc)
  FILE* timeIndex;
  //Main loop stats file (written every time step,
  // contains vital overall statistics such as particle counts)
  FILE* mainStats;

  // Initialise variables
  picConfig = ArcPicConfig();

  double dvt;
  n_aver = 0, n_aver_ion = 0, n_aver_diagn = 0;

  //TODO: This should be backed up
  double induced_cathode_prev = 0.0;

  // Collision sanity checks
  Vec3d mcheck = {0.0,0.0,0.0};
  double echeck = 0.0;

  //Read input.txt
  input();

  int NG = (nr+1)*(nz+1);
  //Sanity check of injection steps
  if (n2inj_step % dt_ion != 0 or n2inj_step == 0) {
    printf("ERROR!! n2inj_step=%i doesn't occur on dt_ion=%i\n", n2inj_step, dt_ion);
    printf("Aborting!\n");
    exit(1);
  }
  if (i2inj_step % dt_ion != 0 or i2inj_step == 0) {
    printf("ERROR!! i2inj_step=%i doesn't occur on dt_ion=%i\n", i2inj_step, dt_ion);
    printf("Aborting!\n");
    exit(1);
  }

  //Initialize multithreading
  if ( numParaThreads > omp_get_num_procs() ||  numParaThreads > omp_get_thread_limit() ) {
      printf("Error in main.cpp during initialization: Can't have numParaThreads=%i > num_procs=%i or thread_limit=%i\n",
           numParaThreads, omp_get_num_procs(), omp_get_thread_limit());
    exit(1);
  }
  else if (numParaThreads < 1) {
    printf("Error in main.cpp during initialization: Can't have numParaThreads=%i < 1\n", numParaThreads);
    exit(1);
  }
  omp_set_num_threads(numParaThreads);

  //File names for continuation
  if ( snprintf( fbckp, LEN_FILENAME, "out/bckp_%03i.dat", CONTINUATION ) >= LEN_FILENAME ) {
    printf("Error detected when making filename for fbckp, generated filename = '%s'\n", fbckp);
    exit(1);
  }
  if ( CONTINUATION > 0 ) {
    if ( snprintf( frestart, LEN_FILENAME, "out/bckp_%03i.dat", CONTINUATION-1 ) >= LEN_FILENAME ) {
      printf("Error detected when making filename for frestart, generated filename = '%s'\n", frestart);
      exit(1);
    }
  }

  picConfig.dr = picConfig.dz;

  double time_start = omp_get_wtime(); //Real walltime [s], relative to some (thread-dependent) starting point

  time_t now_is = my_time();
  printf("*** Starting up 2D Arc-PIC code, date:  *** %s\n", ctime(&now_is) );

  //********INITIALISING...********//

  // STARTING NEW RUN
  if ( CONTINUATION == 0 ) {
    nstepsmin = 0;
    printf("-- Starting new run. (%d steps) -- \n", nstepsmin);

    allocate_arrays( nr, nz, &perm_c_slu, &perm_r_slu, &rhs_slu );

    //Initialize scaling
    picConfig.calc_parameters_2D();

    M_ions[0] =  picConfig.mi_over_me;                 //  H+
    M_ions[1] =  63.546/1.0073*picConfig.mi_over_me-1; //  Cu+
    M_ions[2] =  63.546/1.0073*picConfig.mi_over_me;   //  Cu

    for (int sort=0; sort < NSpecies; sort++) {
      cs_ions[sort] =  sqrt(M_ions[0]/M_ions[sort])*picConfig.vi_0;
      vt_ions[sort] =  sqrt(M_ions[0]/M_ions[sort])*picConfig.v_ti;
    }

    //Create particle arrays
    electrons              = new ParticleSpecies(nr,nz, "e-", 1.0, picConfig.qe);                    //e-
    ionSpecies.push_back    (new ParticleSpecies(nr,nz, "H+", 63.546/1.0073*picConfig.mi_over_me-1, picConfig.qi));//Cu+
    ionSpecies.push_back    (new ParticleSpecies(nr,nz, "Cu+", 63.546/1.0073*picConfig.mi_over_me-1, picConfig.qi));//Cu+
    neutralSpecies.push_back(new ParticleSpecies(nr,nz, "Cu", 63.546/1.0073*picConfig.mi_over_me  , 0.0));//Cu

    rng_initialize(RNGbaseSeed, numParaThreads);
    printf("\n");
    rng_printStatus();
    printf("\n");

    // External circuit elements
    circuit->init();
    //Particle boundary conditions
    pbounds->init(nr, picConfig.Zmin, picConfig.Zmax, picConfig.Rmax);
    //Initial particle distribution
    if (iParts != NULL) {
      iParts->init();
    }

    // Add initial conditions, outputting
    if ( BC == 2 || BC == 3 ) potential_factorise_BC23( nr, nz, NR, NZ, picConfig.dr, picConfig.dz, &L_slu, &U_slu, &perm_c_slu, &perm_r_slu );
    else                        potential_factorise_2D( nr, nz, NR, NZ, picConfig.dr, picConfig.dz, &L_slu, &U_slu, &perm_c_slu, &perm_r_slu );

    // Initialize field and sputtering arrays
    for (int i=0; i<Lastion*(nr+1)*(nz+1); i++) {
      E_ion_r[i]=0.;
      E_ion_z[i]=0.;
    }
    for (int i=0; i<(nr+1)*(nz+1); i++) {
      E_grid_r[i]=0.;
      E_grid_z[i]=0.;
      phi[i]=0.;
    }

    //Initialize time index file
    timeIndex = fopen("out/timeIndex.dat", "w");
    fprintf(timeIndex, "# StepNum SimTime[ns]\n");
    fflush(timeIndex);

    mainStats = fopen("mainStats.dat", "w");
    fprintf(mainStats, "# StepNum SimTime[ns] StepTime[s] n_electrons n_ions n_neutrals [superparticles]\n");
    fflush(mainStats);

    // INITIAL PARTICLE DISTRIBUTION
    if (iParts != NULL) {

      iParts->inject_e(electrons);
      for (auto neutral : neutralSpecies) {
        iParts->inject_n(neutral);
      }
      for (auto ion : ionSpecies) {
        iParts->inject_i(ion);
      }

      if ( OUT_COORD ) {
        // Output coordinates: position and velocity
        if ( BINARY_OUTPUT ) {
          H5::H5File* h5OutFile_initDist = createH5File_timestep( 0, "initDist" );
          H5::Group group_coords = h5OutFile_initDist->createGroup("/COORDS");

          out_coords_2D_h5( electrons, 1, picConfig.Omega_pe, picConfig.dz, "ELECTRONS", group_coords );
          out_coords_2D_h5( ionSpecies[1], dt_ion, picConfig.Omega_pe, picConfig.dz, "IONS", group_coords );
          out_coords_2D_h5( neutralSpecies[0], dt_ion, picConfig.Omega_pe, picConfig.dz, "NEUTRALS", group_coords );
          h5OutFile_initDist->close();
          delete h5OutFile_initDist;
          h5OutFile_initDist = NULL;
        }
        else {
          file_names_2D( 0 );
          out_coords_2D( electrons, 1, picConfig.Omega_pe, picConfig.dz, fr_e );
          out_coords_2D( ionSpecies[1], dt_ion, picConfig.Omega_pe, picConfig.dz, fr_i );
          out_coords_2D( neutralSpecies[0], dt_ion, picConfig.Omega_pe, picConfig.dz, fr_n );
        }
      }

    }

    // INITIAL ENERGY
    kin_pot_en( electrons, ionSpecies[1], neutralSpecies[0],
                &En_e, En_i+1, En_i+Lastion, &En_p, &En_tot,
                1./ionSpecies[1]->mass, 1./neutralSpecies[0]->mass,
                phi, NR, NZ, picConfig.Omega_pe, picConfig.dz );
    printf( "...... Initial energy .................... \n" );

    std::cout   << std::setw(8) << std::setfill('.') << std::left << electrons->name
                << std::right   << std::setfill(' ')
                << " : np = " << std::setw(10) << electrons->GetN()
                << ", En = "  << std::setw(10) << En_e << std::endl;
    for (auto ion : ionSpecies) {
      size_t sort=1; // TODO
      std::cout   << std::setw(8) << std::setfill('.') << std::left << ion->name
                  << std::right   << std::setfill(' ')
                  << " : np = " << std::setw(10) << ion->GetN()
                  << ", En = "  << std::setw(10) << En_i[sort] << std::endl;
    }
    for (auto neutral : neutralSpecies) {
      size_t sort=2; // TODO
      std::cout   << std::setw(8) << std::setfill('.') << std::left << neutral->name
                  << std::right   << std::setfill(' ')
                  << " : np = " << std::setw(10) << neutral->GetN()
                  << ", En = "  << std::setw(10) << En_i[sort] << std::endl;
    }

    printf( "Pot. en..:               En_pot = %-9.5f \n",  En_p );
    printf( "...... Total: En_tot = %-9.6f ........ \n", En_tot);
    printf( "\n");

    H5::H5File* h5OutFile_0 = NULL;
    if ( BINARY_OUTPUT ) {
      std::cout << "Creating HDF5 output file for initial data " << std::flush;
      h5OutFile_0 = createH5File_timestep( 0 );
      std::cout << h5OutFile_0->getFileName() << " ..." << std::flush << std::endl;
    }
    else {
      file_names_2D( 0 );
    }

    // TO TEST
    electrons->UpdateDensityMap( Vcell );
    ionSpecies[1]->UpdateDensityMap( Vcell );
    if ( BINARY_OUTPUT ) {
      H5::Group group_dens = h5OutFile_0->createGroup("/DENSITY");
      out_dens_2D_h5( electrons->densMap,     1, -1., nr, nz, NZ, picConfig.Omega_pe, picConfig.dr, picConfig.dz, "ELECTRONS", group_dens );
      out_dens_2D_h5( ionSpecies[1]->densMap, 1,  1., nr, nz, NZ, picConfig.Omega_pe, picConfig.dr, picConfig.dz, "IONS",      group_dens );
    }
    else {
      out_dens_2D( electrons->densMap,     1, -1., nr, nz, NZ, picConfig.Omega_pe, picConfig.dr, picConfig.dz, fn_e );
      out_dens_2D( ionSpecies[1]->densMap, 1,  1., nr, nz, NZ, picConfig.Omega_pe, picConfig.dr, picConfig.dz, fn_i );
    }
    //TODO: Maybe not needed?
    electrons->ZeroDensityMap();
    for (auto ion : ionSpecies) {
      ion->ZeroDensityMap();
    }
    for (auto neutral : neutralSpecies) {
      neutral->ZeroDensityMap();
    }

    // I. GET POTENTIAL
    H5::Group group_emfield_0;
    if ( BINARY_OUTPUT ) {
      group_emfield_0 = h5OutFile_0->createGroup("/EMFIELD");
      out_phi_2D_h5(phi, 1, nr, nz, NZ, picConfig.Omega_pe, picConfig.dr, picConfig.dz, "POTENTIAL0", group_emfield_0);
    }
    else{
      file_names_2D( 0 );
      out_phi_2D(phi, 1, nr, nz, NZ, picConfig.Omega_pe, picConfig.dr, picConfig.dz, fphi);
    }

    if ( BC == 2 || BC == 3 ) potential_backsolve_BC23( nr, nz, NR, NZ, picConfig.dz, circuit->getU0(), circuit->getUNz(),
                                                        phi, L_slu, U_slu, perm_c_slu, perm_r_slu,
                                                        electrons->densMap, ionSpecies[1]->densMap, &rhs_slu );
    else                        potential_backsolve_2D( nr, nz, NR, NZ, picConfig.dz, circuit->getU0(), circuit->getUNz(),
                                                        phi, L_slu, U_slu, perm_c_slu, perm_r_slu,
                                                        electrons->densMap, ionSpecies[1]->densMap, &rhs_slu );

    printf("\n");

    // II. CALCULATE FIELD
    electric_field_2D( phi, E_grid_r, E_grid_z, E_ion_r, E_ion_z, nr, nz, NR, NZ );

    printf("\n");

    if ( BINARY_OUTPUT ) {
      out_phi_2D_h5(phi, 1, nr, nz, NZ, picConfig.Omega_pe, picConfig.dr, picConfig.dz, "POTENTIAL", group_emfield_0);
      if ( OUT_EFIELD ) {
        out_efield_2D_h5( E_av_z, E_av_r, 1, nr, nz, NZ, picConfig.Omega_pe, picConfig.dr, picConfig.dz, "EFIELD", group_emfield_0 );
      }
    }
    else{
      file_names_2D( 1 );
      out_phi_2D(phi, 1, nr, nz, NZ, picConfig.Omega_pe, picConfig.dr, picConfig.dz, fphi);
    }

    //Initialize reactions etc. just before writing metadata to h5OutFile_0
    picConfig.init_reactions();
    picConfig.print_parameters_2D();

    if ( BINARY_OUTPUT ) {
      //Add metadata
      picConfig.outputfile_addParameterMetadata(h5OutFile_0, nsteps);

      h5OutFile_0->close();
      delete h5OutFile_0;
      h5OutFile_0 = NULL;
    }

  }
  else {
    //Temporary
#warning backup/restart disabled!
    printf("Restart facility not fully implemented!\n");
    exit(1);

    // CONTINUING PREVIOUS RUN
    allocate_arrays( nr, nz, &perm_c_slu, &perm_r_slu, &rhs_slu );

    read_data( frestart );
    re_init();

    timeIndex = fopen("out/timeIndex.dat", "a");
    mainStats = fopen("mainStats.dat", "a");

    // Check the data read in:
    checkbounds_2D( electrons, picConfig.Rmin, picConfig.Rmax, picConfig.Zmin, picConfig.Zmax );
    for (auto ion : ionSpecies) {
      checkbounds_2D( ion, picConfig.Rmin, picConfig.Rmax, picConfig.Zmin, picConfig.Zmax );
    }
    for (auto neutral : neutralSpecies) {
      checkbounds_2D( neutral, picConfig.Rmin, picConfig.Rmax, picConfig.Zmin, picConfig.Zmax );
    }

    // Compute L and U matrices
    if ( BC == 2 || BC == 3 ) potential_factorise_BC23( nr, nz, NR, NZ, picConfig.dr, picConfig.dz, &L_slu, &U_slu, &perm_c_slu, &perm_r_slu );
    else                        potential_factorise_2D( nr, nz, NR, NZ, picConfig.dr, picConfig.dz, &L_slu, &U_slu, &perm_c_slu, &perm_r_slu );

    printf("-- Continuing from old run. (%d steps) -- \n", nstepsmin);

    picConfig.init_reactions();
    picConfig.print_parameters_2D();

    // Continue with loop - electron push
  }


  printf( "*** Beginning main loop *** \n" );
  print_time( omp_get_wtime()-time_start);
  printf("\n");

  double step_time = omp_get_wtime(); // time for each step, written to mainStats.dat

  //********BEGIN LOOP********//
  for (nsteps=nstepsmin; nsteps<=nstepsmax; nsteps++) {

    // III. MOVE PARTICLES
    if ( MAGNETIC ) {
      push_magnetic_2D( electrons, E_grid_r, E_grid_z, picConfig.Bz_ext, picConfig.Bt_ext, 1., NZ );
    }
    else  {
      push_2D         ( electrons, E_grid_r, E_grid_z,                     NZ );
    }

    pbounds->remove_e( electrons );

    if (DODEBUG) {
      if(not checkbounds_2D(electrons, picConfig.Rmin, picConfig.Rmax, picConfig.Zmin, picConfig.Zmax)) {
        std::cerr << "Error detected by checkbounds_2D() for electrons (right after remove_e() ), nsteps = " << nsteps << "." << std::endl;
        if (electrons->GetN() < 100){
          electrons->PrintParticles();
        }
        else {
          std::cerr << " (to many particles to print)" << std::endl;
        }
        exit(1);
      }
    }

    if ( nsteps/e2inj_step*e2inj_step  == nsteps ) {
      pbounds->inject_e(electrons, E_grid_z);
    }

    if (DODEBUG) {
      if(not checkbounds_2D(electrons, picConfig.Rmin, picConfig.Rmax, picConfig.Zmin, picConfig.Zmax)) {
        std::cerr << "Error detected by checkbounds_2D() for electrons, nsteps = " << nsteps << "." << std::endl;
        if (electrons->GetN() < 100){
          electrons->PrintParticles();
        }
        else {
          std::cerr << " (to many particles to print)" << std::endl;
        }
        exit(1);
      }
    }

    // IV. CALCULATE DENSITIES
    electrons->UpdateDensityMap( Vcell );

    if ( nsteps/dt_ion*dt_ion == nsteps ) { //BEGIN IF ION STEP

      scalEion_2D( E_ion_r, E_ion_z, NR, NZ, dt_ion, M_ions );

      // Push and remove

      for (auto ion : ionSpecies) {
        unsigned int sort = 1; // TODO: FIXME!
        if ( MAGNETIC ) {
          push_magnetic_2D( ion, E_ion_r + sort*NG, E_ion_z + sort*NG, picConfig.Bz_ext, picConfig.Bt_ext, -1.*dt_ion/ion->mass, NZ );
        }
        else {
          push_2D( ion, E_ion_r + sort*NG, E_ion_z + sort*NG, NZ );
        }
        pbounds->remove_i(ion, sort);
      }

      for (auto neutral : neutralSpecies) {
        push_neutrals_2D( neutral );
        pbounds->remove_n( neutral );
      }

      // Inject

      if ( nsteps/i2inj_step*i2inj_step  == nsteps ) {
        unsigned int sort = 1; // TODO: FIXME!
        for (auto ion : ionSpecies) {
          pbounds->inject_i(ion,E_grid_z,sort);
        }
      }

      if ( nsteps/n2inj_step*n2inj_step  == nsteps ) {
        for (auto neutral : neutralSpecies) {
          pbounds->inject_n(neutral, E_grid_z);
        }
      }

      //Compute density maps

      for (auto ion : ionSpecies) {
        if (DODEBUG) {
          if(not checkbounds_2D(ion, picConfig.Rmin, picConfig.Rmax, picConfig.Zmin, picConfig.Zmax)) {
            std::cerr << "Error detected by checkbounds_2D() for ion '" << ion->name << "', nsteps = " << nsteps << "." << std::endl;
            if (ion->GetN() < 100){
              ion->PrintParticles();
            }
            else {
              std::cerr << " (to many particles to print)" << std::endl;
            }
            exit(1);
          }
        }

        ion->UpdateDensityMap( Vcell );
      }

      for (auto neutral : neutralSpecies) {
        if (DODEBUG) {
          if(not checkbounds_2D(neutral, picConfig.Rmin, picConfig.Rmax, picConfig.Zmin, picConfig.Zmax)) {
            std::cerr << "Error detected by checkbounds_2D() for neutral '"  << neutral->name << "', nsteps = " << nsteps << "." << std::endl;
            if (neutral->GetN() < 100){
              neutral->PrintParticles();
            }
            else {
              std::cerr << " (to many particles to print)" << std::endl;
            }
            exit(1);
          }
        }

        // neutral density only used for outputting
        if( nsteps >= nav_start ) {
          neutral->UpdateDensityMap(Vcell);
        }
      }

      if( nsteps >= nav_start ) {
        //TODO: All species
        aver_moments_2D( mom_ion + NG, n_aver_ion, ionSpecies[1], nr, nz, NZ );
        aver_moments_2D( mom_ion + Lastion*NG, n_aver_ion, neutralSpecies[0], nr, nz, NZ ); // should be _SN()
        average_2D( ionSpecies[1]->densMap, n_i_av + NG, NR, NZ, n_aver_ion ); // NEW: 25.8.2010
        average_2D( neutralSpecies[0]->densMap, n_i_av + Lastion*NG, NR, NZ, n_aver_ion ); // NEW: 25.8.2010

        n_aver_ion++;
        //printf("Averaging ion moments \n");
      }

      initEion_2D( E_ion_r, E_ion_z, NR, NZ );

    } //END IF ION STEP

    // V. MCC COLLISIONS (IF NOT DISABLED)
    if (DOCOLL == true) {
      // Coulomb collisions
      if (ncoll_el>0 && nsteps/ncoll_el*ncoll_el == nsteps) {
        // e-e Coulomb collisions
        electrons->Order2D();
        coll_el_knm_2D(electrons, nr, nz, NZ, 0, &mcheck, &echeck, ncoll_el);

        // i-i Coulomb collisions
        ionSpecies[1]->Order2D(); //TODO: All species!
        coll_el_knm_2D(ionSpecies[1], nr, nz, NZ, 1, &mcheck, &echeck, ncoll_el);
      }

      // Other collisions
      if(ncoll_ion > 0 && nsteps/ncoll_ion*ncoll_ion == nsteps ) {

        //TODO: This is often not needed!
        electrons->Order2D();
        //TODO: This is often not needed!
        ionSpecies[1]->Order2D(); //TODO: All species!
        neutralSpecies[0]->Order2D(); //TODO: All species!

        //elastic Cu+ Cu collisions
        coll_ion_neutral_noSP_2D( neutralSpecies[0],
                                  ionSpecies[1],
                                  nr, nz, NZ, React_Cup_Cu_el, &mcheck, &echeck );

        //elastic Cu Cu collisions
        coll_n_n_2D( neutralSpecies[0], nr, nz, NZ, React_Cu_Cu, &mcheck, &echeck );

        //elastic el Cu collisions
        coll_el_all_fake_2D( neutralSpecies[0],
                             electrons, nr, nz, NZ, React_Cu_el );

        // e + Cu = Cu+ + 2e
        coll_el_neutrals_2D( neutralSpecies[0],
                             electrons, ionSpecies[1],
                             nr, nz, NZ, React_Cu_ion, &mcheck, &echeck );
      }
    }

    // VI. UPDATE POTENTIAL
    //printf("induced cathode charge elec = %g, ions = %g, deltaQ = %g ",
    //       inducedCharge_cathode(elec, nr_e, e_order), inducedCharge_cathode(ions+NPART, nr_i[1], i_order+NG),
    //       pbounds->getDeltaQ());

    //TODO: Generalize

    double induced_cathode = inducedCharge_cathode(ionSpecies[1]) - inducedCharge_cathode(electrons);

    //positive for negative particles leaving surface and positive coming towards it
    double induced_cathode_delta = induced_cathode - induced_cathode_prev;
    //printf("induced current = %g wall current = %g total current = %g\n",
    //       induced_cathode_delta, pbounds->getDeltaQ(), pbounds->getDeltaQ() + induced_cathode_delta);

    induced_cathode_prev = induced_cathode;

    //circuit->timestep(pbounds->getDeltaQ(), nsteps);
    circuit->timestep(pbounds->getDeltaQ()+induced_cathode_delta, nsteps);

    //Notify the pbounds if this is a output timestep
    pbounds->timestep(nsteps, nsteps >= nav_start && nsteps == nav_start+nav_time);



    if ( BC == 2 || BC == 3 ) potential_backsolve_BC23( nr, nz, NR, NZ, picConfig.dz, circuit->getU0(), circuit->getUNz(),
                                                        phi, L_slu, U_slu, perm_c_slu, perm_r_slu,
                                                        electrons->densMap, ionSpecies[1]->densMap, &rhs_slu );
    else                        potential_backsolve_2D( nr, nz, NR, NZ, picConfig.dz, circuit->getU0(), circuit->getUNz(),
                                                        phi, L_slu, U_slu, perm_c_slu, perm_r_slu,
                                                        electrons->densMap, ionSpecies[1]->densMap, &rhs_slu );

    // VII. CALCULATE FIELD
    electric_field_2D( phi, E_grid_r, E_grid_z, E_ion_r, E_ion_z, nr, nz, NR, NZ );

    // VIII. DIAGNOSTICS AND OUTPUT
    step_time = omp_get_wtime()-step_time;
    fprintf(mainStats, "%08i %010e %f %zu %zu %zu \n",
            nsteps, nsteps*picConfig.Omega_pe*1e9/(56414.6*sqrt(picConfig.n_ref)), step_time,
            electrons->GetN(), ionSpecies[1]->GetN(), neutralSpecies[0]->GetN());
    fflush(mainStats);
    step_time = omp_get_wtime();

    if (nsteps >= diagn_start) {
      aver_diagn_2D( electrons->densMap, diagn_dens, electrons, diagn_Te, diagn_ne, n_aver_diagn, nr, nz, NR, NZ );
      n_aver_diagn++;

      if (nsteps == diagn_start+nav_time) {
        diagn_av_stability( diagn_dens, diagn_Te, diagn_ne, n_aver_diagn, -1., picConfig.v_te, nr, nz,
                            NR, NZ, picConfig.Omega_pe, picConfig.dr, picConfig.dz, nsteps, &check_stab, ferr );
        if (check_stab == 1) {
          //printf("*** ERROR: stability violated. Stopping *** \n");
          //goto The_END; //STOP!!
          printf("*** ERROR: stability violated. *** \n");
          check_stab = 0;
          fflush(stdout);
        }

        n_aver_diagn = 0;
        diagn_start += dt_diagn;
      }
    }

    if (nsteps >= nav_start) {
      if ( OUT_VDF ) {
        //TODO: Generalize for all species

        dvt = 6.66666667/picConfig.v_te;
        vel_dst_along_2D( vdf_ez, vdf_er, vdf_eabs, electrons,         nr, nz, dvt, n_aver );

        dvt = 6.66666667/cs_ions[1];
        vel_dst_along_2D( vdf_iz, vdf_ir, vdf_iabs, ionSpecies[1],     nr, nz, dvt, n_aver );

        dvt = 6.66666667/cs_ions[Lastion];
        vel_dst_along_2D( vdf_nz, vdf_nr, vdf_nabs, neutralSpecies[0], nr, nz, dvt, n_aver );
      }

      aver_moments_2D( mom_el, n_aver, electrons, nr, nz, NZ );
      average_2D( phi,      phi_av, NR, NZ, n_aver );
      average_2D( electrons->densMap,      n_e_av, NR, NZ, n_aver ); // NEW: 25.8.2010
      average_2D( E_grid_z, E_av_z, NR, NZ, n_aver ); // E-field output
      average_2D( E_grid_r, E_av_r, NR, NZ, n_aver ); // E-field output

      n_aver++;
      //printf("Averaging electron moments \n");

      if (nsteps == nav_start) {
        printf("Start averaging .......................... \n");
        printf("Omega_pe*t= %-6.0f, t=%f ns (%d steps) \n",
        nsteps*picConfig.Omega_pe, nsteps*picConfig.Omega_pe*1e9/(56414.6*sqrt(picConfig.n_ref)), nsteps);
      }
      if( nsteps == nav_start+nav_time ) {
        printf("*** Outputting now. (%d steps) *** \n", nsteps);
        print_time( omp_get_wtime()-time_start);
        printf("No. electrons: %zu, ions: %zu, neutrals: %zu \n",
               electrons->GetN(),ionSpecies[1]->GetN(),neutralSpecies[0]->GetN());
        printf("\n");

        fprintf(timeIndex, "%08i %010e\n", nsteps, nsteps*picConfig.Omega_pe*1e9/(56414.6*sqrt(picConfig.n_ref)) );
        fflush(timeIndex);

        if ( BINARY_OUTPUT ) {
          std::cout << "Writing HDF5 output file " << std::flush;
          H5::H5File* h5OutFile = createH5File_timestep( nsteps );
          std::cout << h5OutFile->getFileName() << " ..." << std::flush;

        //Add metadata
        picConfig.outputfile_addParameterMetadata(h5OutFile, nsteps);

        //Density
        H5::Group group_dens = h5OutFile->createGroup("/DENSITY");
        out_dens_2D_h5( n_e_av,              n_aver,    -1., nr, nz, NZ, picConfig.Omega_pe, picConfig.dr, picConfig.dz, "ELECTRONS", group_dens );
        out_dens_2D_h5( n_i_av + NG,         n_aver_ion, 1., nr, nz, NZ, picConfig.Omega_pe, picConfig.dr, picConfig.dz, "IONS",      group_dens );
        out_dens_2D_h5( n_i_av + Lastion*NG, n_aver_ion, 1., nr, nz, NZ, picConfig.Omega_pe, picConfig.dr, picConfig.dz, "NEUTRALS",  group_dens );

        //Electromagnetic fields
        H5::Group group_emfield = h5OutFile->createGroup("/EMFIELD");
        // Electrostatic potential
        out_phi_2D_h5( phi_av, n_aver, nr, nz, NZ, picConfig.Omega_pe, picConfig.dr, picConfig.dz, "POTENTIAL", group_emfield );
        if ( OUT_EFIELD ) {
          out_efield_2D_h5( E_av_z, E_av_r, n_aver, nr, nz, NZ, picConfig.Omega_pe, picConfig.dr, picConfig.dz, "EFIELD", group_emfield );
        }

        //Position & velocity
        if ( OUT_COORD ) {
          H5::Group group_coords = h5OutFile->createGroup("/COORDS");
          out_coords_2D_h5( electrons, 1, picConfig.Omega_pe, picConfig.dz, "ELECTRONS", group_coords );
          out_coords_2D_h5( ionSpecies[1], dt_ion, picConfig.Omega_pe, picConfig.dz, "IONS", group_coords );
          out_coords_2D_h5( neutralSpecies[0], dt_ion, picConfig.Omega_pe, picConfig.dz, "NEUTRALS", group_coords );
        }

        //Velocity moments
        // Average velocity
        H5::Group group_velavg = h5OutFile->createGroup("/VELOCITY_AVERAGE");
        out_vels_2D_h5( mom_el, nr, nz, NZ, picConfig.cs*sqrt(M_ions[0]/M_ions[1]), picConfig.dr, picConfig.dz, "uez", "uer", "uet", group_velavg);
        out_vels_2D_h5( mom_ion + NG, nr, nz, NZ, picConfig.cs*dt_ion*sqrt(M_ions[0]/M_ions[1]), picConfig.dr, picConfig.dz, "uiz", "uir", "uit", group_velavg);
        // Temperature
        H5::Group group_temp = h5OutFile->createGroup("/VELOCITY_TEMP");
        out_temps_2D_h5(       mom_el,        picConfig.cs*sqrt(M_ions[0]/M_ions[1]), picConfig.me_over_mi*M_ions[0]/M_ions[1],
                         nr, nz, NZ, picConfig.dr, picConfig.dz, "Tez", "Ter", "Tet", group_temp);
        out_temps_2D_h5( mom_ion + NG, picConfig.cs*dt_ion*sqrt(M_ions[0]/M_ions[1]), 1.0,
                         nr, nz, NZ, picConfig.dr, picConfig.dz, "Tiz", "Tir", "Tit", group_temp);

          if ( OUT_VDF ) {
            H5::Group group_vdf = h5OutFile->createGroup("/VELOCITY_DIST");
            out_fv_along_2D_h5( vdf_ez, vdf_er, vdf_eabs, nr, nz, "vdfez", "vdfer", "vdfeabs", group_vdf );
            out_fv_along_2D_h5( vdf_iz, vdf_ir, vdf_iabs, nr, nz, "vdfCupz", "vdfCupr", "vdfCupabs", group_vdf );
            out_fv_along_2D_h5( vdf_nz, vdf_nr, vdf_nabs, nr, nz, "vdfCuz", "vdfCur", "vdfCuabs", group_vdf );
          }

          h5OutFile->close();
          delete h5OutFile;
          h5OutFile=NULL;
          std::cout << " done." << std::endl << std::flush;
        }
        else{
          file_names_2D( nsteps );
          out_dens_2D( n_e_av,              n_aver,    -1., nr, nz, NZ, picConfig.Omega_pe, picConfig.dr, picConfig.dz, fn_e );
          out_dens_2D( n_i_av + NG,         n_aver_ion, 1., nr, nz, NZ, picConfig.Omega_pe, picConfig.dr, picConfig.dz, fn_i );
          out_dens_2D( n_i_av + Lastion*NG, n_aver_ion, 1., nr, nz, NZ, picConfig.Omega_pe, picConfig.dr, picConfig.dz, fn_n );

          out_phi_2D( phi_av, n_aver, nr, nz, NZ, picConfig.Omega_pe, picConfig.dr, picConfig.dz, fphi );

          out_vels_2D ( mom_el, nr, nz, NZ, picConfig.cs*sqrt(M_ions[0]/M_ions[1]), picConfig.dr, picConfig.dz, fv_ez, fv_er, fv_et );
          out_temps_2D( mom_el, picConfig.cs*sqrt(M_ions[0]/M_ions[1]), picConfig.me_over_mi*M_ions[0]/M_ions[1],
                        nr, nz, NZ, picConfig.dr, picConfig.dz, fT_ez, fT_er, fT_et);

          out_vels_2D( mom_ion + NG, nr, nz, NZ, picConfig.cs*dt_ion*sqrt(M_ions[0]/M_ions[1]), picConfig.dr, picConfig.dz, fv_iz, fv_ir, fv_it);
          out_temps_2D( mom_ion + NG, picConfig.cs*dt_ion*sqrt(M_ions[0]/M_ions[1]), 1., nr, nz, NZ, picConfig.dr, picConfig.dz, fT_iz, fT_ir, fT_it);

          if ( OUT_VDF ) {
            out_fv_along_2D( vdf_ez, vdf_er, vdf_eabs, nr, nz, fvdf_ez, fvdf_er, fvdf_eabs );
            out_fv_along_2D( vdf_iz, vdf_ir, vdf_iabs, nr, nz, fvdf_iz, fvdf_ir, fvdf_iabs );
            out_fv_along_2D( vdf_nz, vdf_nr, vdf_nabs, nr, nz, fvdf_nz, fvdf_nr, fvdf_nabs );
          }

          if ( OUT_EFIELD ) {
            out_efield_2D( E_av_z, E_av_r, n_aver, nr, nz, NZ, picConfig.Omega_pe, picConfig.dr, picConfig.dz, fEz, fEr );
          }

          if ( OUT_COORD ) {
            // Output coordinates: position and velocity
            out_coords_2D( electrons, 1, picConfig.Omega_pe, picConfig.dz, fr_e );
            out_coords_2D( ionSpecies[1], dt_ion, picConfig.Omega_pe, picConfig.dz, fr_i );
            out_coords_2D( neutralSpecies[0], dt_ion, picConfig.Omega_pe, picConfig.dz, fr_n );
          }
        }

        kin_pot_en( electrons, ionSpecies[1], neutralSpecies[0],
                    &En_e, En_i+1, En_i+Lastion, &En_p, &En_tot,
                    1./M_ions[1], 1./M_ions[2], phi, NR, NZ, picConfig.Omega_pe, picConfig.dz );
        printf( "...... Energy balance .................... \n" );

        std::cout   << std::setw(8) << std::setfill('.') << std::left << electrons->name
                    << std::right   << std::setfill(' ')
                    << " : np = " << std::setw(10) << electrons->GetN()
                    << ", En = "  << std::setw(10) << En_e << std::endl;
        for (auto ion : ionSpecies) {
          size_t sort=1; // TODO
          std::cout   << std::setw(8) << std::setfill('.') << std::left << ion->name
                      << std::right   << std::setfill(' ')
                      << " : np = " << std::setw(10) << ion->GetN()
                      << ", En = "  << std::setw(10) << En_i[sort] << std::endl;
        }
        for (auto neutral : neutralSpecies) {
          size_t sort=2; // TODO
          std::cout   << std::setw(8) << std::setfill('.') << std::left << neutral->name
                      << std::right   << std::setfill(' ')
                      << " : np = " << std::setw(10) << neutral->GetN()
                      << ", En = "  << std::setw(10) << En_i[sort] << std::endl;
        }

        printf( "Pot. en..:               En_pot = %-9.5f \n",  En_p );
        printf( "...... Total: En_tot = %-9.6f ........ \n", En_tot );
        printf( "\n" );

        fflush( stdout );

        // Check densities
        f_dens = SQU(1./picConfig.Omega_pe)/n_aver;
        for (int k=0; k<=nz; k++ ) {
          if ( (n_i_av[NG+k] > 5./f_dens) || (n_e_av[k] < -5./f_dens) ) {
            printf("ERROR: UNDERRESOLVED n_i_av = %.4e > %.4e or n_e_av = %.4e < %.4e for k = %d \n",
                   n_i_av[NG+k],5./f_dens,n_e_av[k],-5./f_dens,k);
            fflush(stdout);
            check_dens = 1;
            //goto The_END;
          }
        }

        n_aver = 0;
        n_aver_ion = 0;

        // Backup data (if we are not under-resolved)
#warning backup/restart disabled!
        //if ( check_dens == 0 ) save_data( fbckp );

        nav_start += nav_dt;
      }
    }

  }
  //********END LOOP********//
  printf("Done!.......................... \n");
  fflush(stdout);

  // The_END: save_data( fbckp );
 The_END: delete_arrays( &perm_c_slu, &perm_r_slu, &rhs_slu );

  delete circuit;
  delete pbounds;

  delete   electrons;

  for (auto ion : ionSpecies) {
    delete ion;
    ion = NULL;
  }
  ionSpecies.clear();

  for (auto neutral : neutralSpecies) {
    delete neutral;
    neutral = NULL;
  }
  neutralSpecies.clear();

  fclose(timeIndex);

  fclose(mainStats);

  printf( "*** Total runtime *** \n" );
  print_time( omp_get_wtime()-time_start);
  printf("\n");
  fflush(stdout);

  return 0;
}

void  print_time( double t ) {
  // t = time in seconds
  double  tt, hour, min, sec;

  t   = (floor)(t+0.5);
  tt  = modf(t/3600., &hour);
  tt  = modf(tt*60.,  &min);
  sec = t-hour*3600.-min*60.;

  printf("Computation time: %6.0f s (%2.2i:%2.2i:%2.2i) \n",
         t, (int)hour, (int)min, (int)sec );
}
