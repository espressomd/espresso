#ifndef GLOBAL_H
#define GLOBAL_H
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

/****************************************
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

/****************************************
 * box dimensions from global.c
 *****************************************/
extern double box_l[3];
/* size of local box. */
extern double local_box_l[3];
/* left corner of local box. */
extern double my_left[3];
/* right corner of local box. */
extern double my_right[3];

/****************************************
 * thermostat data from thermostat.c
 *****************************************/

/** friction coefficient */
extern double friction_gamma;

/****************************************
 * particle data from global.c
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
 * bonded interactions from particle_data.c
 ****************************************/

/** possible values for bonded_ia_type */
#define BONDED_IA_DIHEDRAL 0
#define BONDED_IA_ANGLE    1

typedef union {
  int bonded_ia_type;
  struct {
    int dummy;
  } dihedral;
  struct {
    int dummy;
  } angle;
} Bonded_ia_parameters;

/** field defining the bonded ia types */
extern int n_bonded_ia;
extern Bonded_ia_parameters *bonded_ia_params;

/** reallocate particles bonds */
void realloc_bonds(int index, int size);

/****************************************
 * integration from integrator.c
 ****************************************/

/** time step for integration */
extern double time_step;
/** maximal interaction cutoff. */
extern double max_cut;
/** verlet list skin. */
extern double skin;
/** maximal interaction range (max_cut + skin). */
extern double max_range;
/** maximal interaction range squared. */
extern double max_range2;

/** Flag for integrator.
    Wether to calculate the forces before the first step. */
extern int calc_forces_first;

/****************************************
 * Verlet list from verlet.c
 ****************************************/
extern int   n_verletList;
extern int max_verletList;
extern int    *verletList;

/** Flag for rebuilding the verlet list. */
extern int rebuild_verletlist;

/**********************************************
 * description of global variables
 * add any variable that should be handled
 * automatically in global.c. This includes
 * distribution to other nodes and
 * read/user-defined access from Tcl.
 **********************************************/

/** possible field types */
#define TYPE_INT    0
#define TYPE_DOUBLE 1

/** maximal dimension of a writable datafield */
#define MAX_DIMENSION 64

/* set callback procedure */
typedef int (SetCallback)(Tcl_Interp *interp, void *data);

/* variable descriptor */
typedef struct {
  void        *data;      /* physical address */
  int          type;      /* int or double */
  int          dimension; /* field dimension */
  const char  *name;      /* name assigned in Tcl */
  SetCallback *changeproc;/* procedure called if value should be
			     changed. Maybe ro_callback for
			     non-writeable variables */
} Datafield;

extern const Datafield fields[];

/* identifiers of variables in fields */
#define FIELD_NPROCS 0
#define FIELD_PGRID  1
#define FIELD_LBOXL  2
#define FIELD_BOXL   3
#define FIELD_NTOTAL 4
#define FIELD_NPTYPE 5
#define FIELD_NITYPE 6
#define FIELD_TSTEP  7
#define FIELD_MCUT   8
#define FIELD_SKIN   9
#define FIELD_RANGE 10
#define FIELD_GAMMA 11

/**********************************************
 * misc procedures
 **********************************************/

/** tcl procedure for datafield access */
int setmd(ClientData data, Tcl_Interp *interp,
	  int argc, char **argv);
#endif
