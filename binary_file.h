// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2009; all rights reserved unless otherwise stated.
#ifndef BINARYFILE_H
#define BINARYFILE_H
/** \file binary_file.h
    This file defines a binary file format for the particle data.

The functions in this file are no longer supported by the Espresso
team. Use them at your own risk!

    It is the header file for \ref binary_file.c "binary_file.c" and provides
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
#include <tcl.h>

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
int writemd(ClientData data, Tcl_Interp *interp,
	    int argc, char **argv);

/** Implements the readmd Tcl command. Note that for reading in new particles,
    all three position fields are mandatory. It doesn't matter where they are
    put, readmd searches for them. readmd also takes care to initialize the
    \ref particle_node map and the \ref node_grid, if necessary.
*/
int readmd(ClientData data, Tcl_Interp *interp,
	   int argc, char **argv);
/*@}*/

#endif
