// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2006; all rights reserved unless otherwise stated.
#ifndef CONFIG_H
#define CONFIG_H

/** \file config.h
 
   This file contains all preprocessor flags deciding which 
   \ref features to turn on/off.  It is recommended to turn
   everything off which you do not need in order to optimize the
   performance of Espresso for your problem. There are also quite a
   number of features which are turned off by default since they are
   used only rarely.

   To access the information on the compilation status of the code you
   are working with in your Espresso Tcl-script, use the corresponding
   \ref tcl_features "Tcl-commands".

   If you add a new feature to Espresso, you also have to add the
   corresponding lines in the function \ref compilation_callback and
   to add documentation in <tt>doc/text/features.doc</tt>.
 
   <b>Responsible:</b>
   <a href="mailto:arnolda@mpip-mainz.mpg.de">Axel</a>
*/

#include <archconfig.h>

#define PARTIAL_PERIODIC
#define ELECTROSTATICS 
/* #define ROTATION */
#define EXTERNAL_FORCES
#define CONSTRAINTS
/* #define MASS */
/* #define EXCLUSIONS */
/* #define COMFORCE */
/* #define COMFIXED */
/* #define MOLFORCES */
/* #define BOND_CONSTRAINT */

/* #define TABULATED */
#define LENNARD_JONES
/* #define LJ_WARN_WHEN_CLOSE */
/* #define MORSE */
/* #define LJCOS */
/* #define LJCOS2 */
/* #define BUCKINGHAM */
/* #define SOFT_SPHERE */

/* Note: Activate ONLY ONE bonded angle potential out of the following! */
/* #define BOND_ANGLE_HARMONIC */
#define BOND_ANGLE_COSINE
/* #define BOND_ANGLE_COSSQUARE */

/* #define NEMD */
/* #define NPT */ 
/* #define DPD */
/* #define LB */

/************************************************/
/** \name Default Parameter Settings            */
/************************************************/
/*@{*/

/** CELLS: Default value for the maximal number of cells per node. */
#define CELLS_MAX_NUM_CELLS 216

/** P3M: Default for number of interpolation points of the charge
    assignment function. */
#define P3M_N_INTERPOL 32768

/** P3M: Default for boundary condition: Epsilon of the surrounding
    medium. */
#define P3M_EPSILON 0.0

/** P3M: Default for offset of first mesh point from the origin (left
    down corner of the simulation box. */
#define P3M_MESHOFF 0.5

/** P3M: Default for the number of Brillouin zones taken into account
    in the calculation of the optimal influence function (aliasing
    sums). */
#define P3M_BRILLOUIN 1

/** P3M: Maximal mesh size that will be checked. The current setting
         limits the memory consumption to below 1GB, which is probably
	 reasonable for a while. */
#define P3M_MAX_MESH 128

/** Precision for capture of round off errors. */
#define ROUND_ERROR_PREC 1.0e-14

/** Tiny angle cutoff for sinus calculations */
#define TINY_SIN_VALUE 1e-10
/** Tiny angle cutoff for cosine calculations */
#define TINY_COS_VALUE 0.9999999999

/** maximal number of iterations in the RATTLE algorithm before it bails out. */
#define SHAKE_MAX_ITERATIONS 1000

/*@}*/

/********************************************/
/* \name exported functions of config.c     */
/********************************************/
/*@{*/
#include <tcl.h>

/** callback for version status. */
int version_callback(Tcl_Interp *interp);

/** callback for compilation status. */
int compilation_callback(Tcl_Interp *interp);
/*@}*/

#endif

#if defined LB
#define USE_TEMPORARY
#endif
