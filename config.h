/** \file config.h 
 *
 *  <b>Responsible:</b>
 *  <a href="mailto:arnolda@mpip-mainz.mpg.de">Axel</a>
*/

/** if defined, the code will be slower, but with the \ref #periodic
    array you can choose which coordinates are bound to p.b.c and
    which are not. If not defined, all coordinates are bound to
    p.b.c. 

    Has effect on: \ref per_callback, \ref find_node, \ref #fields, 
    \ref #cells_init and \ref #sort_particles_into_cells.
*/
#define PARTIAL_PERIODIC 

/** if defined, you will get a warning when particles approach nearer than
    0.9 sigma, because then it's likely the integration will blow up.
*/
/* #define LJ_WARN_WHEN_CLOSE */

#define ELECTROSTATICS

/** Compiler flag to enable describing and processing particle orientations.
This will allow to use such particle properties as quart, lambda, and torque. */
/** #define DIPOLAR_INTERACTION */

/** Compiler flag to enable external forces. E.g. apply a fixed external force
    to a particle or fix a particle in space. */
#define EXTERNAL_FORCES

/** Compiler Flag tp enable constraints, eg walls, spheres. 
    See \ref constraint.h and \ref interaction_data.h */
#define CONSTRAINTS

/************************************************/
/** \name Default Parameter Settings            */
/************************************************/
/*@{*/

/** CELLS: Default value for the maximal number of cells per node. */
#define CELLS_MAX_NUM_CELLS 512

/** P3M: Default for number of interpolation points of the charge
    assignment function. */
#define P3M_N_INTERPOL 32768

/** P3M: Default for baundary condition: Epsilon of the surrounding
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

/*@}*/
#ifndef CONFIG_H
#define CONFIG_H

#include <tcl.h>
/** callback for version status. */
MDINLINE int version_callback(Tcl_Interp *interp)
{
  Tcl_AppendResult(interp, "ESPRESSO/tcl_md: Version: 0.99b RC01, Last Change: 16.04.2003", (char *) NULL);
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
#ifdef EXTERNAL_FORCES
  Tcl_AppendResult(interp, "{ EXTERNAL_FORCES } ", (char *) NULL);
#endif
#ifdef CONSTRAINTS
  Tcl_AppendResult(interp, "{ CONSTRAINTS } ", (char *) NULL);
#endif
#ifdef DIPOLAR_INTERACTION
  Tcl_AppendResult(interp, "{ DIPOLAR_INTERACTION} ", (char *) NULL);
#endif
  Tcl_AppendResult(interp, "}", (char *) NULL);
  return (TCL_OK);
}

#endif
