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

/* box dimensions */
extern double box_l[3];
extern double my_left[3];
extern double my_right[3];

/* particle data */
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
extern Particle *particles;

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
