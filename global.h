#ifndef GLOBAL_H
#define GLOBAL_H
/** \file global.h
    This file contains the code for access to globally defined
    variables using the script command setmd.
*/

#include <tcl.h>

/**********************************************
 * description of global variables
 * add any variable that should be handled
 * automatically in global.c. This includes
 * distribution to other nodes and
 * read/user-defined access from Tcl.
 **********************************************/

/** Field is of type integer in \ref Datafield. */
#define TYPE_INT    0
/** Field is of type double in \ref Datafield. */
#define TYPE_DOUBLE 1

/** Maximal size of an array in \ref Datafield. */
#define MAX_DIMENSION 64

/** \type int (SetCallback)(Tcl_Interp *interp, void *data)
    Type for the write callback procedure of \ref Datafield */
typedef int (SetCallback)(Tcl_Interp *interp, void *data);

/** Type describing variables that are accessible via Tcl. */
typedef struct {
  /** Physical address of the variable. */
  void        *data;
  /** Type of the variable, either \ref TYPE_INT or \ref TYPE_DOUBLE.*/
  int          type;
  /** Dimension of the variable. */
  int          dimension;
  /** Name used in the Tcl script. */
  const char  *name;
  /** changeproc is called with an array of type \ref type and
      dimension \ref dimension every time a setmd is done in Tcl script code.
      The procedure is assumed to actually set the value and return TCL_OK or
      not change the value and return TCL_ERROR and an appropriate error message
      in the interpreters result stack. If you use the \ref ro_callback, every time
      an error will be issued, i. e. the value cannot be changed from Tcl.
  */
  SetCallback *changeproc;
} Datafield;

extern const Datafield fields[];

/** index of \ref nprocs in \ref fields */
#define FIELD_NPROCS 0
/** index of \ref processor_grid in \ref fields */
#define FIELD_PGRID  1
/** index of \ref local_box_l in \ref fields */
#define FIELD_LBOXL  2
/** index of \ref box_l in \ref fields */
#define FIELD_BOXL   3
/** index of \ref n_total_particles in \ref fields */
#define FIELD_NTOTAL 4
/** index of \ref n_particle_types in \ref fields */
#define FIELD_NPTYPE 5
/** index of \ref n_interaction_types in \ref fields */
#define FIELD_NITYPE 6
/** index of \ref time_step in \ref fields */
#define FIELD_TSTEP  7
/** index of \ref max_cut in \ref fields */
#define FIELD_MCUT   8
/** index of \ref skin in \ref fields */
#define FIELD_SKIN   9
/** index of \ref max_range in \ref fields */
#define FIELD_RANGE 10
/** index of \ref friction_gamma in \ref fields */
#define FIELD_GAMMA 11

/**********************************************
 * misc procedures
 **********************************************/

/** tcl procedure for datafield access */
int setmd(ClientData data, Tcl_Interp *interp,
	  int argc, char **argv);
#endif
