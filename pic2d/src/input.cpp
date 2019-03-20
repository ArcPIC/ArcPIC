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

  input.cpp:
  Reads in user-given parameters from input.txt

***********************************************************************/

#include  "pic.h"

#define XTRN extern
#include  "var.h"
#include  "circuit.h"
#include  "arcbounds.h"
#undef XTRN

#include "filenames.h"
#include "input.h"

#include  <stdio.h>

#include <cstdlib>
#include <sys/stat.h>
#include <iostream>
//#include <string.h>
#include <cstdlib>
#include <cstring>

#include <assert.h>
#include <limits>

using namespace std;

void input( void ) {

  FILE *in_file;

  //Check that file exists, else we get a segfault
  struct stat st;
  if ( stat("input.txt", &st) ) {
    printf("input(): Could not find file \"input.txt\". Aborting!\n");
    exit(1);
  }
  //Also check that output folder exist, create if neccessary
  if ( stat("out", &st) ) {
    printf("input(): Folder \"out\" did not exist. Creating it now...\n");
    if ( mkdir("out", S_IRWXU) ) {
      printf("input(): Could not create folder \"out\". Abort!\n");
      exit(1);
    }
  }

  in_file = fopen("input.txt","r");

  //SCALING PARAMETERS
  n_ref = parseDouble(in_file,"n_ref");
  T_ref = parseDouble(in_file, "T_ref");
  Ndb   = parseDouble(in_file, "Ndb");


  fscanf(in_file,"%*[^:]%*[:]");
  fscanf(in_file,"%d %d", &nr, &nz);

  dz       = parseDouble(in_file, "dz");
  Omega_pe = parseDouble(in_file, "Omega_pe");

  // TIMESTEPS
  dt_ion    = parseInt(in_file, "dt_ion");
  ncoll_el  = parseInt(in_file, "ncoll_el");
  ncoll_ion = parseInt(in_file, "ncoll_ion");
  dt_diagn  = parseInt(in_file, "dt_diagn");
  nstepsmax = parseInt(in_file, "nstepsmax");


  dt_out   = parseDouble(in_file, "dt_out");
  av_start = parseDouble(in_file, "av_start");
  av_time  = parseDouble(in_file, "av_time");

  // FIELDS, PARTILCES AND BOUNDARY CONDITIONS
  Bz_ext = parseDouble(in_file, "Bz_ext");
  Bt_ext = parseDouble(in_file, "Bt_ext");

  //Injection timesteps
  fscanf(in_file,"%*[^:]%*[:]");
  fscanf(in_file,"%d", &e2inj_step);

  fscanf(in_file,"%*[^:]%*[:]");
  fscanf(in_file,"%d", &n2inj_step);

  fscanf(in_file,"%*[^:]%*[:]");
  fscanf(in_file,"%d", &i2inj_step);

  //New-style particle boundary conditions
  pbounds = ArcBounds::LoadArcBounds(in_file);

  // Circuit parameters
  circuit = Circuit::LoadCircuit(in_file); //Global var in var.h

  //Initial particles
  iParts = InitialParticles::LoadInitialParticles(in_file);

  // OPTIONS
  OUT_COORD     = parseYN(in_file, "OUT_COORD");
  OUT_EFIELD    = parseYN(in_file, "OUT_EFIELD");
  OUT_VDF       = parseYN(in_file, "OUT_VDF");

  MAGNETIC      = parseYN(in_file, "MAGNETIC");

  CONTINUATION  = parseInt(in_file,"CONTINUATION");

  BINARY_OUTPUT = parseYN(in_file, "BINARY_OUTPUT");
  DOCOLL        = parseYN(in_file, "DOCOLL");
  DODEBUG       = parseYN(in_file, "DODEBUG");

  // MISCELLANEOUS
  Ti_over_Te = parseDouble(in_file, "Ti_over_Te");
  mi_over_me = parseDouble(in_file, "mi_over_me");

  fscanf(in_file,"%*[^:]%*[:]");
  fscanf(in_file,"%lu", &RNGbaseSeed);

  BC             = parseInt(in_file, "BC");
  numParaThreads = parseInt(in_file,"numParaThreads");

  fclose (in_file);

}

bool parseYN(FILE* in_file, std::string errorVariable) {
  char foo;
  fscanf(in_file,  "%*[^:]%*[:] %c",  &(foo));
  if      (foo == 'y') return true;
  else if (foo == 'n') return false;
  else {
    cerr << "Error in parsing(): '" << errorVariable << "' has to be either 'y' or 'n'" << endl;
    cerr << "Got: '" << foo << "'" << endl;
    exit(1);
  }
}

double parseDouble(FILE* in_file, std::string errorVariable) {
  double retVal = std::numeric_limits<double>::quiet_NaN();

  char   buffer[LINE_MAXLEN+1];
  memset(buffer, '\0', LINE_MAXLEN+1);

  fscanf(in_file, "%*[^:]%*[:] %sLINE_MAXLEN", buffer);// '%sLINE_MAXLEN':
                                                       // Read from string, maximum NAME_MAXLEN characters
                                                       // (macro variable LINE_MAXLEN is inserted, making it into e.g. '%s300')
                                                       // Note that fscanf will read the given number of characters,
                                                       // and then append a '\0' after that! That's why the allocation is LINE_MAXLEN+1.
  if (buffer[LINE_MAXLEN-1] != '\0') {
    cerr << "Error in parseDouble() when reading '" << errorVariable << "': Possible truncation!" << endl;
    exit(1);
  }

  std::string strbuf = string(buffer);
  size_t idx = 0;

  try {
    retVal = std::stod(strbuf, &idx);
  }
  catch (const std::invalid_argument& ia) {
    cerr << "Invalid argument in parseDouble() when reading variable '" << errorVariable << "'." << endl
         << "Got: '" << buffer << "'" << endl
         << "Expected a floating point number! (exponential notation is accepted)" << endl;
    exit(1);
  }
  catch (const std::out_of_range& ia) {
    cerr << "Invalid argument in parseDouble() when reading variable '" << errorVariable << "'." << endl
         << "Got: '" << buffer << "'" << endl
         << "Value is out of range for double!" << endl;
    exit(1);
  }

  if (idx != strbuf.length()) {
    cerr << "Invalid argument in parseDouble() when reading variable '" << errorVariable << "'." << endl
         << "Got: '" << buffer << "'" << endl
         << "Junk at the end of the argument: '" << strbuf.substr(idx) << "'" << endl;
    exit(1);
  }

  return retVal;
}

double parseInt(FILE* in_file, std::string errorVariable) {
  int retVal = 0;

  char   buffer[LINE_MAXLEN+1];
  memset(buffer, '\0', LINE_MAXLEN+1);

  fscanf(in_file, "%*[^:]%*[:] %sLINE_MAXLEN", buffer);// '%sLINE_MAXLEN':
                                                       // Read from string, maximum NAME_MAXLEN characters
                                                       // (macro variable LINE_MAXLEN is inserted, making it into e.g. '%s300')
                                                       // Note that fscanf will read the given number of characters,
                                                       // and then append a '\0' after that! That's why the allocation is LINE_MAXLEN+1.
  if (buffer[LINE_MAXLEN-1] != '\0') {
    cerr << "Error in parseInt() when reading '" << errorVariable << "': Possible truncation!" << endl;
    exit(1);
  }

  std::string strbuf = string(buffer);
  size_t idx = 0;

  try {
    retVal = std::stoi(strbuf, &idx);
  }
  catch (const std::invalid_argument& ia) {
    cout << "Invalid argument in parseInt() when reading variable '" << errorVariable << "'." << endl
         << "Got: '" << buffer << "'" << endl
         << "Expected an integer!" << endl;
    exit(1);
  }
    catch (const std::out_of_range& ia) {
    cerr << "Invalid argument in parseInt() when reading variable '" << errorVariable << "'." << endl
         << "Got: '" << buffer << "'" << endl
         << "Value is out of range for double!" << endl;
    exit(1);
  }

  if (idx != strbuf.length()) {
    cerr << "Invalid argument in parseInt() when reading variable '" << errorVariable << "'." << endl
         << "Got: '" << buffer << "'" << endl
         << "Junk at the end of the argument: '" << strbuf.substr(idx) << "'" << endl;
    exit(1);
  }

  return retVal;
}

char* readInputSection (FILE* in_file, vector<char*>& options_ret, bool acceptNone) {
  assert (LINE_MAXLEN > NAME_MAXLEN+4); // Should handle '*** NAME'

  char wantName [NAME_MAXLEN+1];
  memset(wantName, '\0', NAME_MAXLEN+1);

  char* retName = NULL;
  vector<char*> options;

  fscanf(in_file, "%*[^:]%*[:]");
  fscanf(in_file, "%sNAME_MAXLEN", wantName);   // '%sNAME_MAXLEN':
                                                // Read from string, maximum NAME_MAXLEN characters
                                                // (macro variable NAME_MAXLEN is inserted, making it into e.g. '%s64')
                                                // Note that fscanf will read the given number of characters,
                                                // and then append a '\0' after that! That's why the allocation is NAME_MAXLEN+1.
  if (wantName[NAME_MAXLEN-1] != '\0') {
    cerr << "Error in readInputSection(): Possible truncation in wantName" << endl;
    exit(1);
  }

  //Move the file pointer past the '\n'
  // foo should be much longer that NAME_MAXLEN (foo also used below)
  char foo[LINE_MAXLEN];
  memset(foo, '\0', LINE_MAXLEN);
  fgets(foo,LINE_MAXLEN,in_file); // fgets will read at max LINE_MAXLEN-1 characters, and then append a '\0'.
                                  // Thus the allocation of LINE_MAXLEN bytes and checking at position end-1 (= LINE_MAXLEN-2).
  if (foo[LINE_MAXLEN-2] != '\0') {
    cerr << "Error in readInputSection(): Possible buffer overflow in foo (looking for EOL in main header)."
         << " Check input.txt, could be lots of trailing whitespace!" << endl;
    exit(1);
  }

  cout << "wantName='" << wantName << "'\n";

  //Loop over config sections
  const unsigned int safetyMax_sections = 100;
  unsigned int safety_section = 0;
  while (true) {
    //Loop sanity check
    if ( safety_section >= safetyMax_sections ) {
      cerr << "Error in readInputSection: Infinite loop in parser, or more than"
           << "safetyMax_sections =" << safetyMax_sections << " config sections." << endl;
      exit(1);
    }

    //Find the right section
    char thisName [NAME_MAXLEN+1];
    memset(thisName, '\0', NAME_MAXLEN+1);

    memset(foo, '\0', LINE_MAXLEN);
    fgets(foo, LINE_MAXLEN,in_file);
    if (foo[LINE_MAXLEN-2] != '\0') {
      cerr << "Error in readInputSection(): Possible buffer overflow in foo (looking for EOL in section scan)."
           << " Check input.txt, could be lots of trailing whitespace!" << endl;
      exit(1);
    }
    sscanf(foo, "*** %sNAME_MAXLEN", thisName);
    if (thisName[NAME_MAXLEN-1]!= '\0') {
      cerr << "Error in readInputSection(): Possible truncation in thisName" << endl;
      exit(1);
    }

    cout << "thisName='" << thisName << "'\n";

    assert(NAME_MAXLEN >= 3);
    if ( ! strncmp(thisName, "END", 3) ) break;

    //Loop over options for thisName
    const unsigned int safetyMax_options = 300;
    unsigned int safety_option = 0;
    while (true) {
      //Loop sanity check
      if ( safety_option >= safetyMax_options ) {
        cerr << "Error in readInputSection: Infinite loop in parser, or more than"
             << "safetyMax_options =" << safetyMax_options << " options while parsing config section"
             << " '" << thisName << "'." << endl;
        exit(1);
      }

      //Get the next line
      char* line = new char[LINE_MAXLEN];
      memset(line, '\0', LINE_MAXLEN);
      fgets(line, LINE_MAXLEN, in_file);
      if ( line[LINE_MAXLEN-2] != '\0') {
        cerr << "Error in readInputSection(): Possible truncation when reading lines of config section"
             << " '" << thisName << "'." <<endl
             << " Check input.txt, could be lots of trailing whitespace!!" << endl;
        exit(1);
      }
      //printf("Got line '%s'\n", line); //DEBUG

      //Check for "///"
      assert(LINE_MAXLEN >= 3);
      if ( ! strncmp(line, "///", 3) ) break;

      //Collect the line
      options.push_back(line);

      safety_option++;
    }

    //Got the right section!
    if ( ! strncmp( thisName, wantName, NAME_MAXLEN ) ) {
      cout << "In readInputSection(), got name = '" << wantName << "'" << endl;
      cout << "Got options:" << endl;
      for (size_t i = 0; i < options.size(); i++) {
        cout << "\t[" << i << "]:" << options[i]; //endl not needed, already have \n at end of lines
      }

      if (retName != NULL) {
        cerr << "Error in readInputSection: Found the wanted name twice in input.txt!" << endl;
        exit(1);
      }

      //Don't return immediatly, need to parse the rest of the sections
      retName = new char[NAME_MAXLEN];
      strncpy(retName, wantName, NAME_MAXLEN);      // We have previously checked that wantName has a \0 in NAME_MAXLEN-1.
      for (size_t i = 0; i < options.size(); i++) {
        char* line_ret = new char[LINE_MAXLEN];
        strncpy(line_ret, options[i], LINE_MAXLEN); // fgets makes sure that options[i] has a \0 in LINE_MAXLEN-1 (and we check that there is no overflow)
        options_ret.push_back(line_ret);
      }
    }

    //Wrong section, clear the options array before next iteration
    for (size_t i = 0; i < options.size(); i++) {
      delete[] options[i];
    }
    options.clear();

    safety_section++;
  }

  if (retName != NULL) {
    //Implicitly also return options_ret
    return retName;
  }
  else if ( ! strncmp( wantName, "None", NAME_MAXLEN ) and acceptNone ) {
    retName = new char[NAME_MAXLEN];
    strncpy(retName, wantName, NAME_MAXLEN);
    return retName;
  }
  //Implicit else
  cerr << "Error in readInputSection(): Couldn't find section corresponding to"
       << "'" << wantName << "'" << endl;
  exit(1);
  return NULL; //avoid compiler warning
}
