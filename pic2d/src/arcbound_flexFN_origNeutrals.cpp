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

  arcbounds_flexFN_origNeutrals.cpp:
  Particle boundary conditions,
  Fowler-Nordheim where beta and alpha are functions of the r-coordinate,
  with sputtering and evaporation similar to Helga's thesis.

***********************************************************************/

#include "arcbound_flexFN_origNeutrals.h"
#include "mydef.h"

#include "random.h"

#define   XTRN extern
#include "var.h"
#include "arrays1.h"
#undef XTRN

#include <cmath>

using namespace std;

// ******** Implementation of FlexFN_twoComp_origNeutrals ******************
FlexFN_twoComp_origNeutrals::FlexFN_twoComp_origNeutrals(std::vector<char*>& options) {
  if (options.size() != 9) {
    cout << "Error in FlexFN_twoComp_origNeutrals(): Expected 9 options, got "
	 << options.size() << endl;
    exit(1);
  }

  sscanf(options[0], "%*[^:]%*[:] %lg", &(this->alpha1));
  sscanf(options[1], "%*[^:]%*[:] %lg", &(this->alpha2));
  sscanf(options[2], "%*[^:]%*[:] %lg", &(this->beta1));
  sscanf(options[3], "%*[^:]%*[:] %lg", &(this->beta2));
  sscanf(options[4], "%*[^:]%*[:] %u",  &(this->idx1));
  sscanf(options[5], "%*[^:]%*[:] %lg", &(this->SEY));
  sscanf(options[6], "%*[^:]%*[:] %lg", &(this->r_Cu_e));

  char foo;
  sscanf(options[7],  "%*[^:]%*[:] %c",  &(foo));
  if      (foo == 'y') this->doHeatspike = true;
  else if (foo == 'n') this->doHeatspike = false;
  else {
    printf("Error in FlexFN_twoComp_origNeutrals::FlexFN_twoComp_origNeutrals(): doHeatspike has to be either 'y' or 'n' \n");
    exit(1);
  }

  sscanf(options[8], "%*[^:]%*[:] %u",  &(this->file_timestep));
}
void FlexFN_twoComp_origNeutrals::print_par() const {
  printf( " - alpha1               %g \n", alpha1);
  printf( " - alpha2               %g \n", alpha2);
  printf( " - beta1                %g \n", beta1);
  printf( " - beta2                %g \n", beta2);
  printf( " - idx1                 %u \n", idx1);
  printf( " - SEY:                 %g \n", SEY);
  printf( " - r_Cu_e:              %g \n", r_Cu_e);
  printf( " - doHeatspike          %c \n", doHeatspike ? 'y' : 'n');
  printf( " - file_timestep        %u \n", file_timestep );
}

void FlexFN_twoComp_origNeutrals::init(unsigned int nr, double zmin, double zmax, double rmax) {
  ArcBounds::init(nr,zmin,zmax,rmax);

  FN_alpha    = new double*[2];
  FN_alpha[0] = new double[nr+1];
  FN_alpha[1] = new double[nr+1];

  FN_beta     = new double*[2];
  FN_beta[0]  = new double[nr+1];
  FN_beta[1]  = new double[nr+1];

  for (unsigned int i = 0; i < nr+1; i++) {
    if (i <= idx1) {
      FN_alpha[0][i] = alpha1;
      FN_beta[0][i]  = beta1;
    }
    else {
      FN_alpha[0][i] = 0.0;
      FN_beta[0][i]  = 0.0;
    }
    FN_alpha[1][i] = alpha2;
    FN_beta[1][i]  = beta2;
  }

  FN_current = new double[nr+1];

  FN_current_sum = new double[nr+1];
  for (unsigned int i = 0; i < nr+1; i++) FN_current_sum[i] = 0.0;


}
void FlexFN_twoComp_origNeutrals::re_init(unsigned int nr, double zmin, double zmax, double rmax) {
  FlexFN_twoComp_origNeutrals::init(nr, zmin, zmax, rmax);
}
FlexFN_twoComp_origNeutrals::~FlexFN_twoComp_origNeutrals() {
  for (size_t i = 0; i < 2; i++) {
    delete[] FN_alpha[i];
    delete[] FN_beta[i];
  }
  delete[] FN_alpha;
  delete[] FN_beta;

  delete[] FN_current;

  delete[] FN_current_sum;
}


void FlexFN_twoComp_origNeutrals::inject_e(ParticleSpecies* pa, double const Ez[]) {
  //Fowler-nordheim two-component injection
  for (unsigned int i = 0; i < nr+1; i++) FN_current[i] = 0;
  for (unsigned int i = 0; i < 2; i++) {
    calcFN_current(Ez,FN_alpha[i],FN_beta[i],FN_current);
  }
  injectFN(pa, FN_current, Ez);
  for (unsigned int i = 0; i < nr+1; i++) FN_current_sum[i] += FN_current[i];

  //SEY at cathode
  double r1, r2; // temp variables from Gaussian RNG
  double v_inj_e = v_te*0.01; //electron injection velocity
  // SEY: inject with incident r-coordinates, empty SEY variables!
  size_t totY = 0;
  for (size_t k=0; k < sput_cathode_SEY.size(); k++ ) {

    pa->ExpandBy( sput_cathode_SEY[k].Y );
    for (int i=0; i < sput_cathode_SEY[k].Y; i++ ) {

      // Gaussian scheme
      do { r1 = sqrt(-2.*log(RAND+1.e-20)); }
      while( r1 > 5. );
      r2 = RAND * TWOPI;
      pa->vr.push_back( r1*cos(r2)*v_inj_e );
      pa->vt.push_back( r1*sin(r2)*v_inj_e );

      // Gaussian scheme
      do { r1 = sqrt(-2.*log(RAND+1.e-20)); }
      while( r1 > 5. );
      pa->vz.push_back( r1*v_inj_e );

      pa->z.push_back( zmin + pa->vz.back() );
      pa->r.push_back( sput_cathode_SEY[k].r + pa->vr.back() );
      if ( pa->r.back() < 0 )      pa->r.back() = 1.e-20;
      else if ( pa->r.back() > nr) pa->r.back() = 2*nr - pa->r.back();

      pa->m.push_back( 1 );

      current_e[0] += 2*( int(pa->r.back()) ) + 1;

      totY++;
    }
  }

  injected_e[0] += totY;
  sput_cathode_SEY.clear();
}

void FlexFN_twoComp_origNeutrals::inject_n(ParticleSpecies* pa, double const Ez[]) {

  double v_inj_i = vt_ions[2];

  // Cathode, neutral evaporation (per-cell)
  double r1, r2;
  for (size_t i = 0; i <= nr; i++) { //Loop over meshpoints/injection cells
    //Calculate number to inject
    double tmp = r_Cu_e * FN_current_sum[i];
    size_t n2inject_evap = size_t ( tmp );
    tmp -= n2inject_evap;
    if ( RAND <= tmp ) n2inject_evap++;
    FN_current_sum[i] = 0.0;

    pa->ExpandBy(n2inject_evap);

    //Inject the particles: Uniform distribution over the cell
    for ( size_t k=0; k<n2inject_evap; k++ ) {
      // Velocity:
      // Gaussian scheme
      do { r1 = sqrt(-2.*log(RAND+1.e-20)); } while( r1 > 5. );
      r2 = RAND * TWOPI;
      pa->vr.push_back( r1*cos(r2)*v_inj_i ); // do not suppress
      pa->vt.push_back( r1*sin(r2)*v_inj_i ); // do not suppress
      // Gaussian scheme
      do { r1 = sqrt(-2.*log(RAND+1.e-20)); } while( r1 > 5. );
      pa->vz.push_back( r1*v_inj_i );

      // uniform random position on annulus
      double a = i > 0  ? (i-0.5) : 0.0;
      double b = i < nr ? (i+0.5) : nr;
      pa->r.push_back( sqrt((b*b-a*a)*RAND+a*a) );
      pa->z.push_back( zmin );

      //Fractional timestep push (euler method, no acceleration)
      double Rp = RAND; //1-R; how (fractionally) far into the timestep
                        // are we at z=0?
      pa->r.back()  += Rp*n2inj_step*pa->vr.back();
      pa->z.back()  += Rp*n2inj_step*pa->vz.back();

      if ( pa->r.back() < 0 )      pa->r.back() = -pa->r.back(); //Reflect on axis
      else if ( pa->r.back() > nr) pa->r.back() = 2*nr - pa->r.back();

      pa->m.push_back( 1 ); //SN;
    }

    injected_n[0] += n2inject_evap;
  }

  // Cathode, sputtering
  inject_sput(pa, sput_cathode, true, v_inj_i);
  // Anode, sputtering only
  inject_sput(pa, sput_anode, false, v_inj_i);
}
void FlexFN_twoComp_origNeutrals::inject_sput(ParticleSpecies* pa, std::vector<Sput> &sput, bool isCathode, double v_inj) {

  double r1, r2;

  size_t m = sput.size();
  size_t tot = 0;
  for ( size_t k=0; k<m; k++ ) {

    pa->ExpandBy( sput[k].Y );
    for ( int i=0; i<sput[k].Y; i++ ) {

      // Gaussian scheme
      do { r1 = sqrt(-2.*log(RAND+1.e-20)); } while( r1 > 5. );
      r2 = RAND * TWOPI;
      pa->vr.push_back( r1*cos(r2)*v_inj );
      pa->vt.push_back( r1*sin(r2)*v_inj );

      // Gaussian scheme
      do { r1 = sqrt(-2.*log(RAND+1.e-20)); } while( r1 > 5. );
      pa->vz.push_back( isCathode ? r1*v_inj : -r1*v_inj );

      pa->z.push_back( (isCathode ? zmin : zmax) + pa->vz.back() );
      pa->r.push_back( sput[k].r + pa->vr.back() );

      if ( pa->r.back() < 0 )      pa->r.back() = 1.e-20;
      else if ( pa->r.back() > nr) pa->r.back() = 2*nr - pa->r.back();

      pa->m.push_back( 1 ); //SN;

      // Sum up total yield for outputting
      tot++;
    }
  }

  injected_n[isCathode ? 0 : 1] += tot;
  sput.clear();
}
void FlexFN_twoComp_origNeutrals::remove_i(ParticleSpecies* pa, unsigned int sort) {
  if(sort != 1){
    if (pa->GetN() != 0) {
      printf("Error detected in FlexFN_twoComp_origNeutrals::remove_i(): %zu particles of sort=%u\n", pa->GetN(),sort);
      exit(1);
    }
    return;
  }

  size_t n_lost = 0;

  // check threshold of sputtering yield, only for ions, only at the cathode
  // threshold = 1e7 A/cm^2
  double threshold = 1.e7*(Ndb*dt_ion*Omega_pe) /
    (6.7193e-12*n_ref*sqrt(T_ref));
  double current [nr];
  for (unsigned int i=0; i<nr; i++ ) current[i] = 0.;

  //Store sputtering particles here, use different formulaes
  // to determine how much is actually sputtered depending on threshold
  static vector<Sput> sput_cathode_temp;

  //Total cathode yield (for sanity check)
  int Ycat_sum = 0;

  // SEY integer and fractional part
  int SEY_i = int ( SEY );
  double SEY_f = SEY - SEY_i;

  for (size_t n=0; n<pa->GetN(); n++ ) {
    //"Infinity"
    if ( pa->r[n] >= rmax ) {
      removed_i[sort][2]++;
      removedIons.push_back(pa->GetOneParticle(n));
      n_lost++;
      continue;
    }
    //Cathode
    else if ( pa->z[n] < zmin ) {
      removed_i[sort][0]++;
      current_i[sort][0] += 2*(int(pa->r[n]))+1;
      current_cathode[ int(pa->r[n]) ] += 1;
      removedIons.push_back(pa->GetOneParticle(n));

      Sput newSput = calc_sput(pa->GetOneParticle(n), cs_ions[sort], current);
      if ( newSput.Y != 0 ) {
        sput_cathode_temp.push_back(newSput);
        Ycat_sum += newSput.Y;
      }

      // SEY = 0.5 = constant, only from ions hitting the cathode
      // SEY with registering r-coordinates
      // Set reasonable threshold for incident ion energy (e.g. 100 eV)
      if ( (SQU(pa->vz[n]) + SQU(pa->vr[n]) + SQU(pa->vt[n]))*T_ref/(2*SQU(cs_ions[sort])) > 100 ) {
        if ( RAND <= SEY_f ) {
          Sput foo;
          foo.r = pa->r[n];
          foo.Y = SEY_i+1;
          sput_cathode_SEY.push_back(foo);
        }
        else if ( SEY_i > 0 ) {
          Sput foo;
          foo.r = pa->r[n];
          foo.Y = SEY_i;
          sput_cathode_SEY.push_back(foo);
        }
      }
      n_lost++;
      continue;
    }
    //Anode
    else if ( pa->z[n] >= zmax ) {
      removed_i[sort][1]++;
      current_i[sort][1] += 2*(int(pa->r[n]))+1;
      current_anode[ int(pa->r[n]) ] -= 1;
      removedIons.push_back(pa->GetOneParticle(n));

      Sput newSput = calc_sput(pa->GetOneParticle(n), cs_ions[sort], NULL);
      if ( newSput.Y != 0 ) sput_anode.push_back(newSput);

      n_lost++;
      continue;
    }
    //Implicit else: keep this particle
    if (n_lost > 0) {
      //Only need to move particles if something has been lost
      pa->z[n-n_lost] = pa->z[n];
      pa->r[n-n_lost] = pa->r[n];
      pa->vz[n-n_lost] = pa->vz[n];
      pa->vr[n-n_lost] = pa->vr[n];
      pa->vt[n-n_lost] = pa->vt[n];
      pa->m[n-n_lost] = pa->m[n];
    }
  }
  //Delete the final n_lost particles
  if (n_lost > 0 ){
    size_t np = pa->GetN();
    pa->z.resize(np-n_lost);
    pa->r.resize(np-n_lost);
    pa->vz.resize(np-n_lost);
    pa->vr.resize(np-n_lost);
    pa->vt.resize(np-n_lost);
    pa->m.resize(np-n_lost);
  }

  // Enhanced yield?
  bool check_enh = false;
  double r, sigma(1e6*nr); //Crash if sigma not initialized
  for (unsigned int i=0; i<nr; i++ ) {
    // Rescale current with area of cell / area Ldb^2
    current[i] /= PI * (2*i + 1) * SQU(dr);

    if ( current[i] >= threshold ) {
      check_enh = true;
      sigma = i+1;
    }
  }

  if ( check_enh == false  or doHeatspike == false) {
    // Yamamura-Tawara fitting for Cu -> Cu
    for (size_t j=0; j<sput_cathode_temp.size(); j++ ) {
      sput_cathode.push_back(sput_cathode_temp[j]);
    }
  }
  else if ( check_enh == true ) {
    // Enhanced yield from MD
    // Generate Gaussian distribution for r-coord's
#warning "Strange model - always generate 1000 Cu's?"
    for (int j=0; j<1000; j++ ) { // 1000/SN
      do { r = sigma * sqrt(-2.*log(RAND+1.e-20)); } while ( r >= nr );
      Sput foo;
      foo.r = r;
      foo.Y = 1;
      sput_cathode.push_back(foo);
    }

    // Print control message
    if ( Ycat_sum > 1000 ) { // 1000/SN
      printf("*** UNDERESTIMATED *** Yamamura sputtering yield Ycat_sum ( %d ) > enhanced sputtering yield ( %d )! \n", Ycat_sum, 1000 );
      fflush( stdout );
    }
  }
  sput_cathode_temp.clear();

  if (check_enh == true and doHeatspike == false) {
    printf("Skipped heatspike sputtering as it is switched off, yield=%d\n", Ycat_sum);
    fflush(stdout);
  }
}
void FlexFN_twoComp_origNeutrals::remove_n(ParticleSpecies* pa){
  size_t n_lost = 0;

  for (size_t n=0; n<pa->GetN(); n++ ) {
    //"Infinity"
    if ( pa->r[n] >= rmax ) {
      removedNeutrals.push_back(pa->GetOneParticle(n));
      removed_n[2]++;
      n_lost++;
      continue;
    }
    //Cathode
    else if ( pa->z[n] < zmin ) {
      removed_n[0]++;
      current_n[0] += 2*(int(pa->r[n]))+1;
      removedNeutrals.push_back(pa->GetOneParticle(n));

      Sput newSput = calc_sput(pa->GetOneParticle(n), cs_ions[NSpecies-1], NULL);
      if ( newSput.Y != 0 ) sput_cathode.push_back(newSput);

      n_lost++;
      continue;
    }
    //Anode
    else if ( pa->z[n] >= zmax ) {
      removed_n[1]++;
      current_n[1] += 2*(int(pa->r[n]))+1;
      removedNeutrals.push_back(pa->GetOneParticle(n));

      Sput newSput = calc_sput(pa->GetOneParticle(n), cs_ions[NSpecies-1], NULL);
      if ( newSput.Y != 0 ) sput_anode.push_back(newSput);

      n_lost++;
      continue;
    }
    //Implicit else: keep this particle
    if (n_lost > 0) {
      //Only need to move particles if something has been lost
      pa->z[n-n_lost] = pa->z[n];
      pa->r[n-n_lost] = pa->r[n];
      pa->vz[n-n_lost] = pa->vz[n];
      pa->vr[n-n_lost] = pa->vr[n];
      pa->vt[n-n_lost] = pa->vt[n];
      pa->m[n-n_lost] = pa->m[n];
    }
  }
  //Delete the final n_lost particles
  if (n_lost > 0 ){
    size_t np = pa->GetN();
    pa->z.resize(np-n_lost);
    pa->r.resize(np-n_lost);
    pa->vz.resize(np-n_lost);
    pa->vr.resize(np-n_lost);
    pa->vt.resize(np-n_lost);
    pa->m.resize(np-n_lost);
  }
}

Sput FlexFN_twoComp_origNeutrals::calc_sput(const Particle pa,
					     const double cs,
					     double* current_enhancedY) {
  double nrg = SQU(pa.p.vz) + SQU(pa.p.vr) + SQU(pa.p.vt);
  nrg *= T_ref/(2*SQU(cs));
  if ( nrg > 23.383 ) { //in eV
    // Register fluxes going through each cell for enhanced Y
    if ( current_enhancedY != NULL ) {
      current_enhancedY[size_t(pa.p.r)] += 1.; // to be rescaled!
    }

    // Reduced energy
    double eps = (4.45394e-06)*nrg;

    // Nuclear stopping cross section
    double sn = 8205*3.441*sqrt(eps)*log(eps + 2.718) /
      (1 + 6.355*sqrt(eps) + eps*(6.882*sqrt(eps) - 1.708));

    // Sputtering yield
    double p = 0.042*(0.2525/3.49) * (sn/(1 + 0.010124*0.1573*pow(eps,0.3))) *
      pow(1 - sqrt(23.383/nrg),2.5);

    // Fractional part of p handled probabilistically
    double Y = int ( p );
    p -= Y;
    if ( RAND <= p ) Y++;

    // Register r-coordinates of bombarding particles
    if ( Y > 0 ) {
      Sput ret;
      ret.r = pa.p.r;
      ret.Y = Y;
      return ret;
    }
  }

  //No sputtering
  Sput ret;
  ret.r = 0.0;
  ret.Y = 0;
  return ret;
}
