#ifndef GLOBAL_H
#define GLOBAL_H
/** \file global.h
    This file contains all globally available variables and related stuff
    like setmd.
*/

#include <tcl.h>
#include "interaction_data.h"
#include "particle_data.h"

/**********************************************
 * global variables
 * put everything here that should be available
 * to more than one module, especially those
 * variables that should be treated as Tcl
 * datafields.
 * Please mark in which fiel the variable is
 * defined.
 **********************************************/

/*****************************************
 * mpi related stuff from communication.c
 *****************************************/

extern int this_node;
extern int nprocs;

/****************************************
 * topology from grid.c
 ****************************************/

extern int processor_grid[3];
extern int pe_pos[3];
extern int neighbors[6];
extern double boundary[6];

/*****************************************
 * box dimensions from global.c
 *****************************************/
extern double box_l[3];
extern double local_box_l[3];
extern double my_left[3];
extern double my_right[3];

/*****************************************
 * thermostat data from thermostat.c
 *****************************************/

extern double friction_gamma;

/****************************************
 * particle data from particle_data.c
 ****************************************/

/** size of local particle array. */
extern int   max_particles;
/** number of particles belonging to that node. */
extern int     n_particles;
/** number of ghost particle belonging to that node. */
extern int     n_ghosts;
/** local particle array. */
extern Particle *particles;

/* total number of particles in the system. */
extern int n_total_particles;
/* used only on master node: particle->node mapping */
extern int  *particle_node;

/** Mapping between particle identity and local index. 
 *    You find the local index of particle i 
 *    at position i of this field. 
 *    A particle that is not in the processors domain 
 *    (including its ghostshell) is marked with -1.
 */
extern int *local_index;

/****************************************
 * nonbonded interactions from interaction_data.c
 ****************************************/

/** number of particle types. */
extern int n_particle_types;

/* size is n_particle_types^2 */
extern IA_parameters *ia_params;

/** number of interaction types. */
extern int n_interaction_types;

/****************************************
 * bonded interactions from global.c
 ****************************************/

/** possible values for bonded_ia_type */
#define BONDED_IA_DIHEDRAL 0
#define BONDED_IA_ANGLE    1

typedef struct {
  int bonded_ia_type;
  union {
    struct {
      int dummy;
    } dihedral;
    struct {
      int dummy;
    } angle;
  } body;
} Bonded_ia_parameters;

extern int n_bonded_ia; 
extern Bonded_ia_parameters *bonded_ia_params;

/****************************************
 * integration from integrator.c
 ****************************************/

extern double time_step;
extern double max_cut;
extern double skin;
extern double max_range;
extern double max_range2;

extern int calc_forces_first;

/****************************************
 * Verlet list from verlet.c
 ****************************************/
extern int   n_verletList;
extern int max_verletList;
extern int    *verletList;

extern int rebuild_verletlist;

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
