#ifndef PARTICLE_DATA_H
#define PARTICLE_DATA_H
/** \file particle_data.h

    <b>Responsible:</b>
    <a href="mailto:arnolda@mpip-mainz.mpg.de">Axel</a>
 *
 *  \todo Put ghost structure as a sub structure into the particle structure, so everybody knows whats going into the ghost communication.
 *
    For more information on particle_data,
    see \ref particle_data.c "particle_data.c"
 */

#include <tcl.h>
#include "config.h"
#include "utils.h"

/************************************************
 * data types
 ************************************************/

typedef struct {
  /** unique identifier for the particle. */
  int    identity;
  /** particle type, used for non bonded interactions. */
  int    type;  

#ifdef ELECTROSTATICS
  /** charge. */
  double q;
#endif

  /** periodically folded position. */
  double p[3];
} ReducedParticle;

/** Struct holding all particle information
 *  of the particles. */
typedef struct {
  ReducedParticle r;

  /** position in the last time step befor last Verlet list update. */
  double p_old[3];
  /** index of the simulation box image where the particle really sits. */
  int    i[3];

  /** force. */
  double f[3];
  /** velocity. */
  double v[3];

  /** bonded interactions list. */
  IntList bl;
} Particle;

/** List of particles. The particle array is resized using a sophisticated
    (we hope) algorithm to avoid unnecessary resizes.
    Access using \ref realloc_particles, \ref got_particle,...
*/
typedef struct {
  /** The particles payload */
  Particle *part;
  /** Number of particles contained */
  int n;
  /** Number of particles that fit in until a resize is needed */
  int max;
} ParticleList;

/** List of reduced particles (e.g. ghost particles). */
typedef struct {
  /** The reduced particles payload */
  ReducedParticle *part;
  /** Number of reduced particles contained */
  int n;
  /** Number of reduced particles that fit in until a resize is needed */
  int max;
} RedParticleList;


/************************************************
 * exported variables
 ************************************************/

/** Highest particle number seen so far. If you leave out some
    particle numbers, this number might be higher than the
    true number of particles. On the other hand, if you start
    your particle numbers at 0, the total number of particles
    is larger by 1.
*/
extern int max_seen_particle;

/** Capacity of the \ref particle_node / \ref local_index. */
extern int  max_particle_node;
/** Used only on master node: particle->node mapping. */
extern int  *particle_node;
/** id->particle mapping on all nodes. */
extern Particle   **local_particles;

/************************************************
 * functions
 ************************************************/

/** implementation of the tcl command part */
int part(ClientData data, Tcl_Interp *interp,
	 int argc, char **argv);

/** initialize a particle list.
 *  Use with care and ONLY for initialization! */
void init_particleList(ParticleList *pList);

/** allocate storage for local particles and ghosts.
    \param plist the list on which to operate
    \param size the size to provide at least. It is rounded
    up to multiples of \ref PART_INCREMENT. */
void realloc_particles(ParticleList *plist, int size);

/** search for a specific particle.
    \param plist the list on which to operate 
    \param id the identity of the particle to search
    \return a pointer to the particle structure or NULL if particle is
    not in this list */
Particle *got_particle(ParticleList *plist, int id);

/** append a particle at the end of a particle List.
    reallocates particles if necessary!
    \param plist List to append the particle to.
    \param part  Particle to append.  
    \return Pointer to new location of the particle. */
Particle *append_particle(ParticleList *plist, Particle *part);

/** remove a particle from one particle List and append it to  another.
    Refill the destList with last particle. 
    reallocates particles if necessary.
    \param destList   List where the particle is appended.
    \param sourceList List where the particle will be removed.
    \param ind        Index of the particle in the sourceList.
    \return Pointer to new location of the particle.
 */
Particle *move_particle(ParticleList *destList, ParticleList *sourceList, int ind);

/** allocate space for a particle.
    \param plist the list on which to operate
    \return the new field index */
Particle *alloc_particle(ParticleList *plist);

/** initialize a reduced particle list (ghosts).
 *  Use with care and ONLY for initialization! */
void init_redParticleList(RedParticleList *pList);

/** allocate storage for reduced particles (ghosts).
    \param plist the list on which to operate
    \param size the size to provide at least. It is rounded
    up to multiples of \ref PART_INCREMENT. */
void realloc_redParticles(RedParticleList *plist, int size);



/** remove bond from particle if possible */
int try_delete_bond(Particle *part, int *bond);

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

/** update \ref max_seen_particle on slave nodes and
    invalidate \ref particle_node. This has to be done
    at the beginning of the integration.
*/
void particle_finalize_data();

/** initialize the \ref local_particles structure. Called from \ref integrate_vv_init */
void local_particles_init();

/** free the \ref local_particles structure. Called from \ref integrate_vv_exit */
void local_particles_exit();

#endif
