#ifndef GLOBAL_H
#define GLOBAL_H
#include <tcl.h>

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
extern int neighbors[6];

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
 * particle data from global.c
 ****************************************/

/* total number of particles in the system. */
extern int n_total_particles;

/** Field to hold particle information
 *  of local particles. */
typedef struct {
  int    identity;
  int    type;

  double p[3];
  double q;

  double v[3];
  double f[3];

  int   n_bonds;
  int max_bonds;
  int    *bonds;
} Particle;

extern int     n_particles;
extern int   max_particles;
extern Particle *particles;

/** first unused particle entry */
extern int   min_free_particle;

extern int     n_ghosts;
extern int   max_ghosts;
extern Particle *ghosts;
/** first unused ghost entry */
extern int   min_free_ghost;

/** Mapping between particle identity and local index. 
 *    You find the local index of particle i 
 *    at position i of this field. 
 *    A particle that is not in the processors domain 
 *    (including its ghostshell) is marked with -1.
 */
extern int *local_index;

/****************************************
 * nonbonded interactions from global.c
 ****************************************/

/** number of particle types. */
extern int n_particle_types;

/** field containing the interaction parameters for
 *  nonbonded interactions. Access via
 * get_ia_param(i, j), i,j < n_particle_types */
typedef struct {
  double LJ_cutoff;
  double LJ_shift;
  double LJ_offset;

  /* don't which else, since electrostatic is different...
     but put rest here too. */
} IA_parameters;

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
 * description of variables from global.c
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
  void        *data;
  int          type;
  int          dimension;
  const char  *name;
  SetCallback *changeproc;
} Datafield;

extern const Datafield fields[];

/**********************************************
 * procedures for access to global variables
 **********************************************/

/** initialize data fields */
void init_data();

/** call if topology (grid, box dim, ...) changed */
void changed_topology();

/*******************************
 * particle storage
 *******************************/

/** allocate storage for local particles and ghosts.
    Given size is rounded up to multiples of
    PART_INCREMENT */
void reallocate_particles(int size);

/** search for a specific particle, returns field index */
int got_particle(int part);

/** add a particle, returns new field index */
int add_particle(int part);

/** free a particle */
void free_particle(int index);

/*******************************
 * bonded ia information
 *******************************/

/** reallocate particles bonds */
void realloc_bonds(int index, int size);

/*******************************
 * nonbonded ia access
 *******************************/

IA_parameters *get_ia_param(int i, int j);

#endif
