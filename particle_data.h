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
  /** particle type, used for non bonded interactions. */
  int    type;

  /** periodically folded position */
  double p[3];
  /** position in the last time step befor last Verlet list update */
  double p_old[3];
  /** index of the simulation box image where the particle really sits */
  int    i[3];
  /** charge */
  double q;

  /** velocity */
  double v[3];
  /** force */
  double f[3];

  /** size of field \ref bonds. */
  int   n_bonds;
  /** allocated size of field \ref bonds. */
  int max_bonds;
  /** field containing the bond information for the particle.
   * 
   * For each bond, where the particle is the main particle, bonds
   * contains the type of the bonded interaction and the identity of
   * the other participating particles.  
   */
  int    *bonds;
} Particle;

/************************************************
 * exported variables
 ************************************************/

/** Size of local particle array. */
extern int   max_particles;
/** Number of particles belonging to that node. */
extern int     n_particles;
/** Number of ghost particle belonging to that node. */
extern int     n_ghosts;
/** Local particle array. 
 *
 *  The figure shows how local particles and ghosts are stored in the
 *  particle array \anchor particle_array.  
 *
 *  \image html particles.gif  "local particle array"
 *  \image latex particles.eps "local particle array" width=8cm
*/
extern Particle *particles;

/** Highest particle number seen so far. If you leave out some
    particle numbers, this number might be higher than the
    true number of particles. On the other hand, if you start
    your particle numbers at 0, the total number of particles
    is larger by 1.
*/
extern int max_seen_particle;

/** Capacity of the particle_node array. */
extern int  n_particle_node;
/** Used only on master node: particle->node mapping. */
extern int  *particle_node;

/** Mapping between particle identity and local index. 
    You find the local index of particle i at position
    i of this field. 
    A particle that is not in the nodes domain 
    (including its ghostshell) is marked with -1.
*/
extern int *local_index;

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
/** (re)allocate storage for particle bonds.
    Bonds are located in the structure \ref Particle .
    \param part     local index of particle where the bonds are located. 
    \param new_size New size of the particles[].bonds field.
 */
void realloc_part_bonds(int part, int new_size);
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
