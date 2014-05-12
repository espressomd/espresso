/*
  Copyright (C) 2010,2012,2013 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
    Max-Planck-Institute for Polymer Research, Theory Group
  
  This file is part of ESPResSo.
  
  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/
#ifndef BINARYFILE_TCL_H
#define BINARYFILE_TCL_H
/** \file binary_file_tcl.hpp
    This file defines a binary file format for the particle data.

The functions in this file are no longer supported by the Espresso
team. Use them at your own risk!

    It is the header file for \ref binary_file_tcl.cpp "binary_file_tcl.c" and provides
    functions to read and write binary particle data to Tcl channels.
    It also defines the structures and constants necessary to interprete/
    generate particle data files in this format. THE FILE FORMAT IS HARDWARE
    DEPENDENT, SINCE RAW DATA TYPES ARE WRITTEN!!!

    <p>
    The file format conists of the following:
    <ol>
    <li> \ref MDHeader with \ref MDHeader::magic set to \ref MDMAGIC ("MD01") without trailing 0.
    <li> then \ref MDHeader::n_rows chars out of
    \ref POSX,\ref POSY,\ref POSZ,\ref VX,\ref VY,\ref VZ,\ref FX,\ref FY,\ref FZ,
    \ref Q,\ref TYPE which determine the data that is contained in the following and
    correspond to the respective fields of \ref Particle.
    <li> now follows the particle data, always starting with one integer giving the
    particles identity (this row is NOT defined before!), then the data for the fields
    as defined before. The data types
    for \ref POSX,\ref POSY,\ref POSZ,\ref VX,\ref VY,\ref VZ,\ref FX,\ref FY,\ref FZ,
    \ref Q are double, for \ref TYPE integer.
    <li> The file is finished with a single integer -1 (replacing the particle identity
    of a subsequent particle).
    </ol>
*/
#include "parser.hpp"

/** This string is to be put in the \ref MDHeader::magic field of \ref MDHeader
    to allow unique identification of binary packed MD data.
*/
#define MDMAGIC "MD01"

/** The header of the binary file format. Magic is used to identify the file and should have
    a value of \ref MDMAGIC ("MD01") without trailing 0. */
struct MDHeader {
  /** Magic Identifier. Must be \ref MDMAGIC ("MD01") without trailing 0. */
  char magic[4];
  /** Number of data rows contained in the following data. */
  int  n_rows;
};

/** \name Field Codes.
    Possible field codes to follow \ref MDHeader. */
/*@{*/
/** Row contains the x component of the position. */
#define POSX  0
/** Row contains the y component of the position. */
#define POSY  1
/** Row contains the z component of the position. */
#define POSZ  2
/** Row contains the x component of the velocity. */
#define VX    3
/** Row contains the y component of the velocity. */
#define VY    4
/** Row contains the z component of the velocity. */
#define VZ    5
/** Row contains the x component of the force. */
#define FX    6
/** Row contains the y component of the force. */
#define FY    7
/** Row contains the z component of the force. */
#define FZ    8
/** Row contains the mass. */
#define MASSES 9
/** Row contains the charge. */
#define Q     10
/** Row contains the type. */
#define TYPE  11
 /** Row contains the x component of the magnetic moment. */
 #define MX    12
 /** Row contains the y component of the magnetic moment. */
 #define MY    13
 /** Row contains the z component of the magnetic moment. */
 #define MZ    14
/** Row contains the solvation free energies (SHANCHEN) */
#define SOLVATION 15

/*@}*/

/**************************************************************
 * functions
 **************************************************************/

/** \name Exported Functions */
/*@{*/
/** Implements the writemd Tcl command. The first argument gives the
    channel to write to, all subsequent arguments give the fields to write.
    The order of the data in the file is given by the argument order.
    The names of the field are the same as above in lowercase. */
int tclcommand_writemd(ClientData data, Tcl_Interp *interp,
	    int argc, char **argv);

/** Implements the readmd Tcl command. Note that for reading in new particles,
    all three position fields are mandatory. It doesn't matter where they are
    put, readmd searches for them. readmd also takes care to initialize the
    \ref particle_node map and the \ref node_grid, if necessary.
*/
int tclcommand_readmd(ClientData data, Tcl_Interp *interp,
	   int argc, char **argv);
/*@}*/

#endif
