#ifndef GLOBAL_H
#define GLOBAL_H
/** \file global.h
    This file contains the code for access to globally defined
    variables using the script command setmd. \ref add_vars "Here"
    you can find details on how to add new variables in the interpreter's
    space.

    <b>Responsible:</b>
    <a href="mailto:arnolda@mpip-mainz.mpg.de">Axel</a>
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
  /** Dimension of the variable. Limited to \ref MAX_DIMENSION */
  int          dimension;
  /** Name used in the Tcl script. */
  const char  *name;
  /** changeproc is called with an array of type \ref Datafield#type
      and dimension \ref Datafield#dimension every time a setmd is
      done in Tcl script code.  The procedure is assumed to actually
      set the value and return TCL_OK or not change the value and
      return TCL_ERROR and an appropriate error message in the
      interpreters result stack. 
  */
  SetCallback *changeproc;
  /** Minimal number of characters needed for identification. */
  int min_length;
} Datafield;

extern const Datafield fields[];

/** \name Field Enumeration
    These numbers identify the variables given in
    \ref #fields for use with \ref mpi_bcast_parameter.
*/
/*@{*/
/** index of \ref box_l in \ref #fields */
#define FIELD_BOXL        0  
/** index of \ref #cell_grid in  \ref #fields */
#define FIELD_CELLGRID    1
/** index of \ref #cell_size in  \ref #fields */
#define FIELD_CELLSIZE    2
/** index of \ref friction_gamma in \ref #fields */
#define FIELD_GAMMA       3
/** index of \ref local_box_l in \ref #fields */
#define FIELD_LBOXL       4
/** index of \ref max_cut in \ref #fields */
#define FIELD_MCUT        5
/** index of \ref transfer_rate  in \ref #fields */
#define FIELD_MAXNUMCELLS 6
/** index of \ref max_seen_particle in \ref #fields */
#define FIELD_MAXPART     7
/** index of \ref max_range in \ref #fields */
#define FIELD_MAXRANGE    8
/** index of \ref max_skin in  \ref #fields */
#define FIELD_MAXSKIN     9
/** index of \ref n_nodes in \ref #fields */
#define FIELD_NNODES      10
/** index of \ref n_total_particles in  \ref #fields */
#define FIELD_NPART       11
/** index of \ref n_particle_types in \ref #fields */
#define FIELD_NPARTTYPE   12
/** index of \ref node_grid in \ref #fields */
#define FIELD_NODEGRID    13
/** index of \ref #periodic in \ref #fields */
#define FIELD_PERIODIC    14
/** index of \ref #skin in \ref #fields */
#define FIELD_SKIN        15
/** index of \ref #temperature in \ref #fields */
#define FIELD_TEMPERATURE 16
/** index of \ref sim_time in  \ref #fields */
#define FIELD_SIMTIME     17
/** index of \ref time_step in \ref #fields */
#define FIELD_TIMESTEP    18
/** index of \ref timing_samples in  \ref #fields */
#define FIELD_TIMINGSAMP   19
/** index of \ref transfer_rate  in \ref #fields */
#define FIELD_TRANSFERRATE 20
/** index of \ref rebuild_verletlist in \ref #fields */
#define FIELD_VERLETFLAG   21
/** index of \ref verlet_reuse in  \ref #fields */
#define FIELD_VERLETREUSE  22
/*@}*/

/**********************************************
 * misc procedures
 **********************************************/

/// Implements the Tcl command \ref tcl_setmd. It allows to modify simulation parameters
int setmd(ClientData data, Tcl_Interp *interp,
	  int argc, char **argv);

/** Implements the Tcl command \ref #tcl_code_info.  It provides information on the
    Version, Compilation status and the debug status of the used
    code. */
int code_info(ClientData data, Tcl_Interp *interp,
	 int argc, char **argv);

#endif
