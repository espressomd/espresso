#ifndef GLOBAL_H
#define GLOBAL_H
/** \file global.h
    This file contains the code for access to globally defined
    variables using the script command setmd.

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

/** \name Field Enumeration
    These numbers identify the variables given in
    \ref fields for use with e. g. \ref mpi_bcast_parameter.
*/
/*@{*/
/** index of \ref n_nodes in \ref fields */
#define FIELD_NNODES 0
/** index of \ref node_grid in \ref fields */
#define FIELD_NGRID  1
/** index of \ref local_box_l in \ref fields */
#define FIELD_LBOXL  2
/** index of \ref box_l in \ref fields */
#define FIELD_BOXL   3
/** index of \ref max_seen_particle in \ref fields */
#define FIELD_MAXPART 4
/** index of \ref n_particle_types in \ref fields */
#define FIELD_NITYPE 5
/** index of \ref sim_time in  \ref fields */
#define FIELD_SIM_TIME 6
/** index of \ref time_step in \ref fields */
#define FIELD_TIME_STEP  7
/** index of \ref max_cut in \ref fields */
#define FIELD_MCUT   8
/** index of \ref skin in \ref fields */
#define FIELD_SKIN   9
/** index of \ref max_range in \ref fields */
#define FIELD_RANGE  10
/** index of \ref friction_gamma in \ref fields */
#define FIELD_GAMMA 11
/** index of \ref rebuild_verletlist in \ref fields */
#define FIELD_VERLET 12
/** index of \ref p3m_struct::bjerrum in \ref fields */
#define FIELD_BJERRUM 13
/** index of \ref p3m_struct::alpha  in \ref fields */
#define FIELD_P3M_ALPHA 14
/** index of \ref p3m_struct::r_cut in \ref fields */
#define FIELD_P3M_RCUT 15
/** index of \ref p3m_struct::mesh in \ref fields */
#define FIELD_P3M_MESH 16
/** index of \ref p3m_struct::cao  in \ref fields */
#define FIELD_P3M_CAO 17
/** index of \ref p3m_struct::epsilon  in \ref fields */
#define FIELD_P3M_EPSILON 18
/** index of \ref p3m_struct::mesh_off  in \ref fields */
#define FIELD_P3M_MESH_OFF 19
/** index of \ref transfer_rate  in \ref fields */
#define FIELD_TRANSFER_RATE 20
/** index of \ref transfer_rate  in \ref fields */
#define FIELD_MAXNUMCELLS 21
/** index of \ref periodic in \ref fields */
#define FIELD_PERIODIC 22
/** index of \ref temperature in \ref fields */
#define FIELD_TEMPERATURE 23
/** index of \ref lj_force_cap in \ref fields */
#define FIELD_LJFORCECAP 24
/** index of \ref start_time in  \ref fields */
#define FIELD_START_TIME 25
/** index of \ref n_total_particles in  \ref fields */
#define FIELD_N_PART 26
/*@}*/

/**********************************************
 * misc procedures
 **********************************************/

/** Implements the Tcl command setmd. It allows to modify simulation parameters.
    These parameters must be declared in \ref fields. If you write your own
    callback, REMEMBER THE OTHER NODES!!! The value will only be set by you
    and only on the master node. Using \ref mpi_bcast_parameter may be useful. */
int setmd(ClientData data, Tcl_Interp *interp,
	  int argc, char **argv);
#endif
