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
/* Integration */
extern double time_step;
extern double max_cut;
extern double skin;
extern double max_range;
extern double max_range2;

/* Verlet list */
extern int   n_verletList;
extern int max_verletList;
extern int    *verletList;

extern int rebuild_verletlist;

/** initialize data fields */
void init_data();

/** search for a specific particle, returns field index */
int got_particle(int part);

/** add a particle, returns new field index */
int add_particle(int part);

#endif
