#ifndef PARTICLE_DATA_H
#define PARTICLE_DATA_H
#include <tcl.h>

/************************************************
 * data types
 ************************************************/

/** Field to hold particle information
 *  of local particles. */
typedef struct {
  int    identity;
  int    type;

  /* periodically folded position */
  double p[3];
  double p_old[3];
  /* index of the simulation box image where the particle really sits */
  int    i[3];
  double q;

  double v[3];
  double f[3];

  int   n_bonds;
  int max_bonds;
  int    *bonds;
} Particle;

/************************************************
 * functions
 ************************************************/

/** tcl procedure for particle access */
int part(ClientData data, Tcl_Interp *interp,
	 int argc, char **argv);

/** allocate storage for local particles and ghosts.
    Given size is rounded up to multiples of
    PART_INCREMENT */
void realloc_particles(int size);

/** search for a specific particle, returns field index */
int got_particle(int part);

/** add a particle, returns field index */
int add_particle(int part);

/** allocate space for a particle, returns field index */
int alloc_particle();

/** fold particle coordinates to primary simulation box */
void fold_particle(double pos[3],int image_box[3]);

/** unfold particle coordinates to physical position */
void unfold_particle(double pos[3],int image_box[3]);

/** free a particle */
void free_particle(int index);

/** add particle to particle->node map */
void map_particle_node(int part, int node);

/** rebuild particle->node map from scratch */
void build_particle_node();

/** update n_total_particles on slave nodes
 *  and invalidate particle_node */
void particle_finalize_data();

#endif
