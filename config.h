// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2004; all rights reserved unless otherwise stated.

/** \file config.h 
 *
 *  This file contains all Precompiler Flags deciding which features
 *  of Espresso to turn on/off. It is recommended to turn everything
 *  off which you do not need in order to optimize the performance of
 *  Espresso for your problem. 
 *
 *  There are also quite a number of features which are turned off by
 *  default since they are used only rarely.
 *
 *  You can get information on the compilation status of the code you
 *  are working with by using the tcl command \ref tcl_code_info
 *  "code_info" within your tcl_script. It is highly recommended to
 *  store this information with your simulation data in order to
 *  maintain the reproducibility of your results.
 *
 *  If you add a new compile flag you also have to add the
 *  corresponding lines in the function \ref compilation_callback.
 *
 *  \warning The documentation of this file only contains the
 *  precompiler flags that you have turned on. See the source code of
 *  this file for more options!
 *
 *  <b>Responsible:</b>
 *  <a href="mailto:arnolda@mpip-mainz.mpg.de">Axel</a>
*/



/** if defined, the code will be slower, but with the \ref #periodic
    array you can choose which coordinates are bound to p.b.c and
    which are not. If not defined, all coordinates are bound to
    p.b.c. 

    Has effect on: \ref per_callback, \ref #fields, and functions in
    \ref domain_decomposition.c, \ref grid.c, \ref interaction_data.c,
    \ref layered.c, \ref statistics_chain.c
*/
#define PARTIAL_PERIODIC

/** if defined, you will get a warning when particles approach nearer than
    0.9 sigma, because then it's likely the integration will blow up.
*/
/* #define LJ_WARN_WHEN_CLOSE */

#define ELECTROSTATICS

/** Compiler flag to enable describing and processing particle orientations.

This will allow to use such particle properties as quart, omega, and torque. */
/* #define ROTATION */

/** Compiler flag to enable external forces. E.g. apply a fixed external force
    to a particle or fix a particle in space. */
#define EXTERNAL_FORCES

/** Compiler Flag to enable constraints, eg walls, spheres. 
    See \ref constraint.h and \ref interaction_data.h */
#define CONSTRAINTS

/** Compiler Flag to enable COMFORCE potential */
/* #define COMFORCE */

/** Compiler Flag to enable COMFIXED potential */
/* #define COMFIXED */

/************************************************/
/** \name available short--ranged potentials
    For optimization it might be useful to switch
    off the ones you don't need */
/*@{*/

/** to use tabulated potential*/
/* #define TABULATED */

/** Lennard-Jones */
#define LENNARD_JONES

/** Lennard-Jones with cosine tail */
/* #define LJCOS */

/*@}*/

/************************************************/
/** \name Switches for bonded interactions      */
/************************************************/
/*@{*/

/* NOTE: Turn on one and only one of the following switches!!! */

/** Harmonic bond angle potential:      V = 1/2 k (phi - phi0)^2 */
#define BOND_ANGLE_HARMONIC
/** Cosine bond angle potential:        V = k (1+cos(phi-phi0)) */
/* #define BOND_ANGLE_COSINE */
/** Cosine square bond angle potential: V = 1/2 k (cos(phi)-cos(phi0))^2 */
/* #define BOND_ANGLE_COSSQUARE */

/*@}*/

/***********************************************************/
/** \name Simulation methods, Integrators and Thermostats  */
/***********************************************************/
/*@{*/

/** NEMD (Non Eqilibrium Molecular Dynamics).
    This is used to perform shear simulations */
/* #define NEMD */

/** Allows to use (N,p,T)-ensembles during integration as well */
/* #define NPT */

/** DPD Thermostat (Dissipative Particle Dynamics) 
    Flag needed only because DPD acts like a short range potential
*/
/* #define DPD */

/*@}*/

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

/** Precision for capture of round off errors. */
#define ROUND_ERROR_PREC 1.0e-14

/** Tiny angle cutoff for sinus calculations */
#define TINY_SIN_VALUE 1e-10
/** Tiny angle cutoff for cosine calculations */
#define TINY_COS_VALUE 0.9999999999

/*@}*/
#ifndef CONFIG_H
#define CONFIG_H

#include <tcl.h>
/** callback for version status. */
MDINLINE int version_callback(Tcl_Interp *interp)
{
  Tcl_AppendResult(interp, "ESPRESSO: v1.6.3b (Icheb), Last Change: 12.06.2004", (char *) NULL);
  return (TCL_OK);
}

/** callback for compilation status. */
MDINLINE int compilation_callback(Tcl_Interp *interp)
{
  Tcl_AppendResult(interp, "{ Compilation status ", (char *) NULL);
#ifdef PARTIAL_PERIODIC
  Tcl_AppendResult(interp, "{ PARTIAL_PERIODIC } ", (char *) NULL);
#endif
#ifdef LJ_WARN_WHEN_CLOSE
  Tcl_AppendResult(interp, "{ LJ_WARN_WHEN_CLOSE } ", (char *) NULL);
#endif
#ifdef ELECTROSTATICS
  Tcl_AppendResult(interp, "{ ELECTROSTATICS } ", (char *) NULL);
#endif
#ifdef ROTATION
  Tcl_AppendResult(interp, "{ ROTATION } ", (char *) NULL);
#endif
#ifdef EXTERNAL_FORCES
  Tcl_AppendResult(interp, "{ EXTERNAL_FORCES } ", (char *) NULL);
#endif
#ifdef CONSTRAINTS
  Tcl_AppendResult(interp, "{ CONSTRAINTS } ", (char *) NULL);
#endif
#ifdef COMFORCE
  Tcl_AppendResult(interp, "{ COMFORCE } ", (char *) NULL);
#endif
#ifdef COMFIXED
  Tcl_AppendResult(interp, "{ COMFIXED } ", (char *) NULL);
#endif
#ifdef TABULATED
  Tcl_AppendResult(interp, "{ TABULATED } ", (char *) NULL);
#endif
#ifdef LENNARD_JONES
  Tcl_AppendResult(interp, "{ LENNARD_JONES } ", (char *) NULL);
#endif
#ifdef LJCOS
  Tcl_AppendResult(interp, "{ LJCOS } ", (char *) NULL);
#endif
#ifdef BOND_ANGLE_HARMONIC
  Tcl_AppendResult(interp, "{ BOND_ANGLE_HARMONIC } ", (char *) NULL);
#endif
#ifdef BOND_ANGLE_COSINE
  Tcl_AppendResult(interp, "{ BOND_ANGLE_COSINE } ", (char *) NULL);
#endif
#ifdef BOND_ANGLE_COSSQUARE
  Tcl_AppendResult(interp, "{ BOND_ANGLE_COSSQUARE } ", (char *) NULL);
#endif
#ifdef NEMD
  Tcl_AppendResult(interp, "{ NEMD } ", (char *) NULL);
#endif
#ifdef NPT
  Tcl_AppendResult(interp, "{ NPT } ", (char *) NULL);
#endif
#ifdef DPD
  Tcl_AppendResult(interp, "{ DPD } ", (char *) NULL);
#endif
  Tcl_AppendResult(interp, "}", (char *) NULL);
  return (TCL_OK);
}

#endif
