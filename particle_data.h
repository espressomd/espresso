#ifndef PARTICLE_DATA_H
#define PARTICLE_DATA_H
/** \file particle_data.h
    For more information on particle_data,
    see \ref particle_data.c "particle_data.c"
 */

#include <tcl.h>
/************************************************
 * data types
 ************************************************/

/** Struct holding all particle information
 *  of the particles. */
typedef struct {
  /** unique identifier for the particle */
  int    identity;
  int    type;

  /** periodically folded position */
  double p[3];
  /** position in the last time step for Verlet list*/
  double p_old[3];
  /** index of the simulation box image where the particle really sits */
  int    i[3];
  /** charge */
  double q;

  /** velocity */
  double v[3];
  /** force */
  double f[3];

  int   n_bonds;
  int max_bonds;
  int    *bonds;
} Particle;

/************************************************
 * functions
 ************************************************/

/** implementation of the tcl command part */
int part(ClientData data, Tcl_Interp *interp,
	 int argc, char **argv);

/** allocate storage for local particles and ghosts.
    \param size the size to provide at least. It is rounded
    up to multiples of \ref PART_INCREMENT. */
void realloc_particles(int size);

/** search for a specific particle.
    \param part the identity of the particle to search
    \return its field index or -1 if particle is not on this node */
int got_particle(int part);

/** add a particle and initialize it.
    \param part the identity of the particle to add
    \return the new field index */
int add_particle(int part);

/** allocate space for a particle.
    \return the new field index */
int alloc_particle();

/** fold particle coordinates to primary simulation box.
    \param pos the position...
    \param image_box and the box

    Both pos and image_box are I/O,
    i. e. a previously folded position will be folded correctly.
*/
void fold_particle(double pos[3],int image_box[3]);

/** unfold particle coordinates to physical position.
    \param pos the position...
    \param image_box and the box

    Both pos and image_box are I/O, i.e. image_box will be (0,0,0)
    afterwards.
*/
void unfold_particle(double pos[3],int image_box[3]);

/** add particle to \ref particle_node.
    This procedure is only used on the master node in script mode.
    \param part the identity of the particle
    \param node the node it is stored on
*/
void map_particle_node(int part, int node);

/** rebuild \ref particle_node from scratch.
    After a simulation step \ref particle_node has to be rebuild
    since the particles might have gone to a different node.
*/
void build_particle_node();

/** update \ref n_total_particles on slave nodes and
    invalidate \ref particle_node. This has to be done
    at the beginning of the integration */
void particle_finalize_data();

#endif
