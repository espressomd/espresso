#ifndef GLOBAL_H
#define GLOBAL_H
/**********************************************
 * global variables
 **********************************************/

/* mpi related stuff from communication.c */
extern int this_node;
extern int nprocs;

/* topology from communication.c */
extern int processor_grid[3];
extern int neighbors[6];

/** box dimensions (primary simulation box). */
extern double box_l[3];
/** size of local box. */
extern double local_box_l[3];
/** left corner of local box. */
extern double my_left[3];
/** right corner of local box. */
extern double my_right[3];

/* particle data */
/** total number of particles in the system. */
extern int n_total_particles;

typedef struct {
  int    identity;
  int    type;

  double p[3];
  double q;

  double v[3];
  double f[3];

  int n_bond;
  int *bonds;
  int *bond_type;
} Particle;

extern int     n_particles;
extern int     n_ghosts;
extern int   max_particles;
/** Field to hold particle information of local particles. */
extern Particle *particles;

/** Mapping between particle identity and local index. 
 *    You find the local index of particle i 
 *    at position i of this field. 
 *    A particle that is not in the processors domain 
 *    (including its ghostshell) is marked with -1.
 */
extern int *local_index;

/** number of particle types. */
extern int n_particle_types;
/** number of interaction types. */
extern int n_interaction_types;

/* Integration */
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

/* ghost communication */
extern int double_sided;

/* Verlet list */
extern int   n_verletList;
extern int max_verletList;
extern int    *verletList;

/** Flag for rebuilding the verlet list. */
extern int rebuild_verletlist;
/** Flag for integrator. 
    Wether to calculate the forces befor the first step. */
extern int calc_forces_first;

/** initialize data fields */
void init_data();

/** search for a specific particle, returns field index */
int got_particle(int part);

/** add a particle, returns new field index */
int add_particle(int part);

#endif
