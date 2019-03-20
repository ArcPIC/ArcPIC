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

  input.h:
  Header file for input.cpp

***********************************************************************/

#ifndef INPUT_H
#define INPUT_H

#include <cstdio>
#include <vector>

//Toplevel input.txt parser
void input( void );

/* This reads a section from the input file on the form
 *
 * (something) : <name_to_use> <short comment bla bla>
 * *** <name1>
 * <Option 1> : <value> <comment bla bla>
 * <Option 2> : <value> <comment bla bla>
 * <Option 3> : <value> <comment bla bla>
 * ///
 * *** <name2> <short comment bla bla>
 * <Option 1> : <value> <comment bla bla>
 * (etc)
 * ///
 * *** END
 *
 * Arguments:
 * - in_file: pointer to the input.txt file to read from
 * - options_ret: Return-by-argument (see below)
 * - acceptNone: Set to true to accept the wantName "None"
 *               corresponding to no input section.
 *               "None" is the returned as retName.
 *
 * Return values:
 * - char* containing the name of selected section
 * - vector<char*> with the option values
 *
 * It is the callers responsibility that char*'s are deleted.
 */
char* readInputSection (FILE* in_file, std::vector<char*>& options_ret, bool acceptNone=false);

//Parse a y/n flag, return a bool
// exit(1) with error message if not 'y' or 'n'.
bool parseYN(FILE* in_file, std::string errorVariable);

//Parse a double-value flag
// exit(1) with error message if not it cannot be successfully converted.
double parseDouble(FILE* in_file, std::string errorVariable);

//Parse a int-value flag
// exit(1) with error message if not it cannot be successfully converted.
double parseInt(FILE* in_file, std::string errorVariable);

//Buffer sizes used in the parsing functions
#define NAME_MAXLEN 64  //Maximum allowed length of Name
#define LINE_MAXLEN 300 //max length of line

#endif
